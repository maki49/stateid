# 激发态不可约表示开发交接（2026-10-07）

## Resume 入口

在 `/home/fortneu/work/codelab/stateid` resume 时，先阅读本文件、`docs/abacus-conventions.md`、`docs/tutorials/03-nv-minus.md` 与 `tests/test_symmetry.py`。目标是从真实 ABACUS Γ 点 NAO 数据得到可审计的轨道及 TDA 根子空间 irrep；不要只按近简并或 band 编号输出标签。

本文件为独立审计 fork 的交接，未实现核心功能，未修改 io、tests 或 README。审计时主线程正在添加 `io/abacus_transitions.py`、对应测试及教程，请保留那些在途改动。当前版本 pyproject 为 0.2.0。

## 结论

数学核心已具备，真实数据适配链尚不完整。当前 NV 下载数据可以支持“前两根以 ↓126→↓127/128 为主、近简并且横向偏振，与候选 E 态相容”，不能严格证明这些 KS 轨道为 a1/e 或激发态为多电子 ³E。主要瓶颈是同几何 C/S/AO 标签、自动 AO 对称操作、原生 MPI LR 向量恢复与参考态相位。

## 已有 API 与限制

- `io.read_wavefunctions` / `read_gamma_wavefunctions`：单帧 Γ 实数或 k 点复数文本；`OrbitalData.coefficients` 为 `(nao, nband)`，明确 spin，文件 band 标签一基。不支持追加多帧和二进制。
- `AbacusReader.read_overlap`：stateid NPY 或 ABACUS 新旧文本 CSR；可 Fourier 到 S(k)。`read_output` 尚未实现。
- `LRData`：x 为 `(nph, nstate)`，可有 y，`transitions` 为 `(nph,5)` 的 k、hole spin、occ、particle spin、virt；`read_lr_npz` 可读手动导出数据。`read_lr_eigenvectors` 原生 ABACUS 分片仍拒绝解析。
- `symmetry.project_representation`：非正交 AO 下 `D=Q† S T Q`，区分系数操作 T 与矩阵元 M=ST，返回 character、closure、metric error。
- `symmetry.analyze_c3v`：要求 E/C3/C3²/sv1/sv2/sv3 六操作，验证 metric、子空间闭合、C3v 群关系及同类 character 后匹配；已有 `[1,1,1]` A1、`[1,1,-1]` A2、`[2,-1,0]` E。
- `particle_hole_operation`：单自旋、不截断完整 ph 张量积，`chi_ref * kron(D_occ.conj(), D_virt)`；要求一维纯参考态。不能直接适配 openshell 双自旋/任意稀疏 ph 表或 full LR。
- `harmonics.magnetic_order` / `complex_to_real_basis` / `rotation_in_real_basis`：实球谐基变换；不是从 STRU 构造 AO 操作的工具。
- 主线程新增 `read_transition_analysis`：读取阈值裁剪后的实 TDA 组成、偶极、振子强度，不重归一、不推断自旋多重度或 irrep；不可代替完整 X。

## 当前 NV 数据与缺失材料

数据根目录：`/home/fortneu/work/tex/lr-grad/run_data/18_NV_diamond/222/pbe`。

可用：各阶段 INPUT.remote_run/STRU.remote_run、KS eig_occ、逐步 excitation energy、每步四个 amplitude rank 分片、avg/JT 最后一步覆盖式 `trans_analysis_updown_tda.dat`。未发现 WFC_NAO/NAO wavefunction 文本、S(R) CSR 或 AO shell 标签文件。`eig_occ` 只有能级与占据，无法恢复轨道 character。

avg 最后一步 thresholded 组成：root 0 的 ↓126→↓127/128 权重约 0.828112+0.135772；root 1 对应 0.135977+0.828194。遗漏约 3.6% 不能凭裁剪表恢复相位/完整态或严格群闭合。文件只对应最后一步，不能据此验证 GS 几何 step0 成分。

准备最小静态高对称 GS 几何重算（不再 relax），明确输出同次 C↑/C↓、S(R)、准确 STRU、轨道文件与 AO 元数据、全部 X、构建版本与并行信息；先核实实际版本 INPUT 参数，而非在交接里猜输出开关。优先此数据采集，不必先实现通用群库。保存 atom/species 顺序、晶胞单位、同一离子步、文件 hash 与 occupation/window 索引。

## 缺口 1：晶体几何到 AO 操作 T(R)

建议新增 `io/structure.py`、`symmetry/ao.py`，不把解析器塞进 character 匹配函数。

需要：STRU 晶胞和坐标、元素/轨道文件映射；每个 AO 的 atom_index、ell、radial zeta/index、m（ABACUS `(0,+1,-1,...)`）、全局 AO index。C 的行不能仅凭 13 AO/原子猜排列。DZP 可含多径向 shell，旋转只能在同一 atom/ell/zeta block 内混合 m。

可用 spglib 获取分数坐标 `(W,t)` 及原子周期映射；转换主动笛卡尔旋转时明确列/行向量，记录原点、平移和返回晶格矢量。Γ 点返回晶格矢量相位为 1，但仍应记录映射，方便未来 k 点。含镜面时处理 improper rotation：在每个 ell 块使用 inversion parity (-1)^ell；先明确 Wigner D 主动/被动约定。

构造 T 的每列表示 R|phi_nu>，须验证 atom permutation species/shell 一致、T†ST=S、T(R1R2)=T(R1)T(R2)。不要仅以实 character 检验转置或逆：加入不交换操作及直接 s/p/d 函数点值测试。真实结构近似对称时报告每操作最大原子偏差和 symprec，禁止放宽容差后默认为严格对称。

## 缺口 2：真实 LR MPI 分片适配

参考实际计算开发树：`/home/fortneu/lr-grad/abacus-develop/source/source_esolver/esolver_lr_lcao_tddft.cpp`，约 730–750 行 writer，985–1055 行 setup/index。仓库 `docs/abacus-conventions.md` 引用另一路 `/home/fortneu/abacus/abacus-develop`，不能混用版本。

- `Excitation_Amplitude_openshell_stepN_1.dat` 的末尾 1 是 `my_rank+1`，不是 root；每个分片包含全部 roots。
- X state-major。每 root 的 local 长度为 `nk*(local_size_up+local_size_down)`；openshell 单一 X 容器，先 up 所有 k，再 down 所有 k。
- 单 spin/k block 是 `(nvirt,nocc)` 的二维 block-cyclic 分布；本地 virt 行快，offset=`local_occ*local_nvirt+local_virt`；必须用 global2local row/col 恢复，而不是按 rank 顺序拼接。
- 当前完整窗口 ↑ nocc=128/nvirt=691，↓ nocc=126/nvirt=693，共 nph=175766（Γ）；6 roots。band window offset也必须记录，不要推广为所有案例都从 band1 开始。
- 需要每 rank process-grid 坐标、BLACS grid/rank 映射、block size、source process coordinates、k ordering、spin/window 数据。仅“四个分片”不能推断 2×2 grid 或 block size；缺 metadata 应显式失败。

建议 sidecar manifest 或明确导出 NPZ，记录这些 metadata 及 source hash；读取器检测重复/遗漏 global ph、每 rank 长度和全部 root 范数。用串行与非等尺寸 2D 多 rank fixture 对照，覆盖 nocc_up≠nocc_down、nk>1、末尾短 block，和错误把 rank 当 root 的负例。

## 缺口 3：根子空间的物理 irrep

高对称结构同时检查电子参考态是否保留群（几何对称不等于 KS 解对称）。先验证 Cσ†SCσ≈I，再分别得到完整 invariant occupied/virtual window 的 Dσ(R)。截断窗口切开简并壳层必须报不闭合。

TDA ph 操作可先按每 spin 做矩阵无关 contraction，避免构造 175766² dense matrix（约 247 GB 仅实 double，complex 更大）。按 i outer/a inner 将 x reshape 为 `(nocc,nvirt)`，与 `kron(D_occ.conj(),D_virt)` 一致的主动作用为 `D_occ.conj() @ X @ D_virt.T`；先用小矩阵逐项和 dense kron 交叉验证。双 spin blocks 联合，在一个统一参考态相位下组合。

取前两根整体 Q，不逐根标 ex/ey；计算每操作 2×2 `D_root=Q†U_ph Q`，report root orthonormality、泄漏 closure、group relations、同类一致性和 character 残差。满足 `[2,-1,0]` 才称 E 子空间。随机 U(2) 混合后 character/closure 应不变，单 E 分量应失败。近简并能量是辅助证据，不是 irrep 判据。

参考态相位 `chi_ref(R)` 不能默认 1 后输出总多电子标签。对于不变 occupied alpha/beta 空间，一维 Slater reference 的空间相位来自 `det(D_occ_up)*det(D_occ_down)`；验证模为1、群乘法与纯参考态条件。泛化需处理分数占据或非单一 determinant，首版应明确拒绝。

NV 的 a1/e 是 KS 轨道标签；a1*⊗e=E 是激发算符，结合参考态 A2 后 A2⊗E=E 是总空间 irrep。³ 的自旋上标需参考 Ms、S²/自旋污染及 spin-conserving LR 的适用边界，不能由 beta 激发或 E character 单独推出。full LR X/Y 不宜直接用 TDA 欧氏概率/闭合，应另设具有明确 norm 的实现。

JT 几何近似 Cs=C1h 时 E 可关联到 A′+A″，不能在真实 C1 几何上强行恢复 E。未来添加 Cs character table、明确的镜面 operation validator 和 parent-child correlation；先保留 C3v 阶段验证。

## 最小里程碑与验收

1. **采集可复现高对称数据**：C/S/STRU/AO manifest，一次 Γ 静态重算；明确 ROOT energies/geometry provenance。可先手动 AO 元数据，但 schema 需验证。
2. **AO operations builder**：小型非正交 s/p/d 多 atom fixture，周期边界 mapping、proper/improper、active inverse 负例；通过 metric 与完整群关系。
3. **真实轨道 irrep**：↓126 `[1,1,1]`，↓127/128 整体 `[2,-1,0]`；若失败输出诊断而非猜标签。完整 occupied/virtual invariance单独检查。
4. **LR adapter + matrix-free subspace analysis**：MPI fixture 与序列参考一致、phase/window/map严格校验；NV root pair `[2,-1,0]` 或清晰拒绝，随机 gauge 稳定。
5. **物理 evidence report**：JSON + Markdown，分开 geometry/orbital/excitation-operator/total-state/spin证据，包含 residual、tolerances、missing data。偶极对前两根使用整体响应张量 `sum d_s d_s†`，报告相对 NV 轴 longitudinal fraction 和 transverse-plane eigenvalues；避免逐根偶极方向受简并 gauge 影响。

运行现有测试建议 `PYTHONPATH=src python -m unittest discover -s tests`；新增测试必须覆盖真实布局和错误路径，不能只复写实现。第一版无需推广到所有空间群、多 k、SOC或 antiunitary。完成前应检查主线程 changes 是否已落地，不覆盖 `abacus_transitions.py`。

## 可以直接用于下次 resume 的提示

> 阅读 docs/development/excited-state-irrep-handoff-2026-10-07.md 和现有 conventions/tests，先核实并采集同几何 C/S/STRU/AO 标签。实现 Γ 点 AO symmetry operation builder 及有明确 metadata 的 LR adapter，再以 matrix-free TDA 根对子空间投影验证 C3v E。区分轨道 irrep、参考态相位、总空间 irrep与自旋；对缺数据和非闭合空间明确失败。保留现有 transition-analysis parser，不修改 ABACUS 源码，先给出可审计 fixture与误差报告。

## 补充：实际 ABACUS 输出参数与逐步同步审计

审计源码 HEAD `d6f3b3b99`（`/home/fortneu/lr-grad/abacus-develop`），2026-10-07。本节替代前述“需核实输出开关”的待办，未修改源码或算例。

推荐未来 Γ 点 NV 静态/relax 输入添加：

```text
out_wfc_lcao  1
out_hsr       1 12
out_stru      1
out_freq_ion  1
out_app_flag  0
out_mul       1
out_wfc_lr    1
```

`out_mat_hs2 1 12` 是可用 legacy alias，可替代 out_hsr；不要两个同时设。优先 out_hsr 文本 1，第二整数 precision。`out_wfc_lcao 1` 文本波函数精度内部固定 8，并不随 CSR precision 改变。`out_element_info` 非完整 AO mapping，通常无需开；它输出元素径向/赝势信息。

### 逐步行为（不是仅最后离子步）

- LR 每离子步先调用 KS runner，因此 C/S 对应当步 SCF 参考态和当步结构：`source/source_esolver/esolver_lr_lcao_tddft.cpp:714`。
- KS after_scf 调用统一 LCAO 输出：`source/source_esolver/esolver_ks_lcao.cpp:625,646`；门控是 `source/source_io/module_ctrl/ctrl_scf_lcao.cpp:108-125`。out_freq_ion=1 则每步输出；0 也会在每次 KS 收敛（或 scf_nmax）输出 C/S，不能误解为“只输出 relax 最后一步”。正整数 N 为 modulo 步抽样，首步0包括，最后一步未必命中；不用它作为“仅最终”开关。
- out_app_flag=0 波函数写到 OUT/WFC（目录创建 `source/source_io/module_parameter/read_input.cpp:287`），文件名含 g=istep+1，每帧独立文本；例 Γ spin2 `WFC/wfs1g1_nao.txt` 和 `wfs2g1_nao.txt`，第二步 g2。命名 `source/source_base/module_out/filename.cpp:105-110,135`，实际 writer `source/source_io/module_wf/write_wfc_nao.cpp:269-287`。
- out_app_flag=1 波函数同名追加，每帧重新 nbands/nlocal header，无明确离子步标记：`write_wfc_nao.cpp:285`、`:69-85`。stateid 当前单帧 reader 拒绝此多帧形式。不要选它做新数据采集。
- 文本 CSR 在 out_app_flag=0 下每步独立 `srg1_nao.csr`、`srg2_nao.csr`，H 为 hrs1g1/hrs2g1等；命名 `source/source_io/module_hs/hsr_writer.cpp:32-58`。header 是 `--- Ionic Step istep+1 ---`，包含 ucell：`:86-113`。out_app_flag=1 同名 sr_nao.csr 逐帧追加。注意当前 text writer 对 istep>0 使用 append 开文件，即使文件名已有 g；新输出目录每帧正常唯一，但重复在已有目录运行时 g2及以后的文件可能追加重复帧。重跑用新 suffix/干净目录，parser需检查帧。
- Γ only 的 CSR 是已折叠到 R=(0,0,0) 的 SΓ/HΓ，带 representation note：`hsr_writer.cpp:290-293`。可以直接用于 Γ metric，不能当未折叠 S(R) 再做多 k。非 Γ 才保留真正 R 分解。CSR 是全局 gather，不需 MPI rank 重建：`:326-350`。
- out_stru1：relax 的 STRU_NOW 每步覆盖，STRU{istep+1} 由 out_freq_ion 抽样；是能量刚计算的移动前构型，`source/source_relax/relax_driver.cpp:205-228`。STRU_FINAL 在退出时写，`:264-275`。循环先 esolve、stru_out、再 relax move，`:58-64`；若未收敛且退出发生在移动后，FINAL可能是尚未计算的新构型，不能盲配最终 C/S。使用当步 STRU编号及 header geometry-evaluated 信息。scf/nscf只输出FINAL，无STRU1，但波函数/CSR仍g1。
- LR amplitude 与 energy 由 out_wfc_lr1 每次求解写；excited relax文件含 step=istep 零基，末尾 rank+1（不是root）：`esolver_lr_lcao_tddft.cpp:733-748,786,828`。**LR step0 ↔ C/S g1 ↔ STRU1；step13 ↔ g14 ↔ STRU14**。out_freq_ion不控制LR writer，所以抽样N>1会导致部分LR步没有C/S/STRU；推荐1。
- `trans_analysis_updown_tda.dat` 并非逐步保留文件，仍可能覆盖到最后一步；不能代替逐步完整 X。

### AO labels 已有输出：用 out_mul1，不需新增标签参数

`ctrl_scf_lcao.cpp:556-567` 调用 Mulliken `cal_mag`，在 `source/source_io/module_mulliken/cal_mag.h:49-53` rank0写 `OUT/Orbital`。写出逻辑 `source/source_cell/cell_index.cpp:264-284`，包含 species、l、实球谐分量 index m、zeta(Z+1)、sym。顺序为 atom-major、atom内AO顺序；每行对应全局AO，可以以数据行序号确定global AO index。

**文件header的 #io 注释误称 orbital index，实际数值是 zero-based atom iat，对同原子所有AO重复。** m并非有符号磁量子数，是0..2l的实球谐index（`:215-225`）；需通过 stateid magnetic_order(l)映射，而不是直接当负/正m。zeta字段完整保留，可正确识别径向shell。Orbital一般覆盖同名文件，但原子/基组顺序在常规relax中固定，可一次保存并hash。out_dos或out_proj_band也有另一路 Orbital输出，不必额外开启它们获取标签。

### MPI LR manifest 仍是缺口

参数 `nb2d` 真实存在，控制2D block size：`source/source_io/module_parameter/read_inp_estruc.cpp:1142-1149`。0为自动，保存INPUT里的0不能替代运行时真实block size。`Parallel_2D` 能提供grid、coord、descriptor，`source/source_base/parallel_2d.cpp:115-130,149-170`；LR复用KS BLACS context与block size，`esolver_lr_lcao_tddft.cpp:990-1002`。

当前没有发现可输出完整 LR布局sidecar的INPUT开关。常规日志有维度/窗口但未见逐rank BLACS grid coord/全局rank映射；不应宣称上述flags已解决分片重建。若保留原生MPI amplitudes，后续需源码新增manifest（实际block size、grid/coord/rank映射、每spin local dimensions、window/k metadata、列行次序），或有完整验证的独立导出工具。串行小fixture可先消除此依赖；nb2d显式值只是减少未知，仍不足保证跨版本rank映射。

建议立刻用现有flags保存C/S/STRU/Orbital/逐步X，后续stateid实现AO操作与布局adapter。高对称静态GS几何的一次计算是最小可用数据；无需为获取标签而重新整条JT弛豫。将实际输出文件加入下载白名单，特别是WFC目录、srg*.csr、Orbital、STRU数字文件及INPUT.info。

## 新增真实小样例入口：7 原子 NV 单胞（主线程，2026-10-07）

主线程已生成同几何真实 C/S/AO 标签：

- 根目录 `/home/fortneu/work/tex/lr-grad/checks/nv-state-fd/single_state`。
- `OUT.ABACUS/WFC/wfs1g1_nao.txt`、`WFC/wfs2g1_nao.txt`：91 AO / 91 bands，alpha occupied16、beta occupied14。
- `OUT.ABACUS/srg1_nao.csr`、`OUT.ABACUS/Orbital`：同第一离子步矩阵与 AO 标签。
- 几何取根目录 `STRU`；写本节时第一步 `STRU1` 仍待完整力计算结束后产生，应在使用时复核文件与步标签。
- 联合读取验证报告 `/home/fortneu/work/tex/lr-grad/checks/nv-state-fd/metadata-joint-check.json`：max|C†SC-I| alpha=1.894e-8、beta=1.388e-8。

正在计算3个 TDA roots、out_wfc_lr=1、串行；完成后的原生幅度可成为不需要 MPI 拼接的最小真实 irrep 开发 fixture。**7 原子周期模型不是稀释的真实 NV 缺陷模型，不能套用63原子的126/127/128轨道编号或直接迁移其物理归属。** 它用于 IO、AO 操作、对称性投影和FD方法正确性的小系统验证。解析后自行按实际占据/能级/character选候选空间，禁止在代码中硬编码大型NV轨道编号。

波函数文本精度 `source/source_io/module_wf/write_wfc_nao.cpp:80,172` 硬编码 setprecision(8)；不能假定 out_ndigits 会提升 WFC 精度。12位CSR配合8位WFC已得到上述1e-8级正交误差，设irrep/闭合容差时需考虑实际文件精度与几何近似对称性，分别报告几何偏差、metric与群关系误差，不用无依据的机器精度门槛。
