# 08 · 用真实非 Γ 点 NV⁻ 输出验证 Bloch 对称性

本节对应 2026-10-07 在 WSL 本地完成的一次 PBE 单点计算。
它把同次计算的 STRU、Orbital、C(k)、S(R)、H(R) 接到 stateid，
并检查数值闭合、群关系和特征标，而不是由能量简并直接赋标签。

## 1. 计算与可追溯数据

计算目录：

~~~text
$HOME/work/LR-TDDFT-cal/grad/18_NV_diamond/111/pbe/stateid_non_gamma_20261007/
~~~

构型由该 NV 项目已有的 gen_stru.py 生成：1×1×1 常规金刚石晶胞，
删除一个 C，并用 N 替换其 [111] 近邻，得到 6 C + 1 N。
生成时名义晶格常数为 3.5686645447981364 Å；
原始 STRU 的 lat0=1.889726 Bohr 是近似换算，按该 ABACUS 的
BOHR_TO_A=0.5291770，实际晶格长度为 3.5686628871624562 Å。
不作结构弛豫。这是高缺陷浓度的周期模型，用于验证 IO/对称性；
其能带编号和能量不能替代 63 原子 NV⁻ 的 126/127/128。

| 设置 | 本次取值 |
|---|---|
| 程序 | ABACUS v3.11.0-beta10，运行日志 commit d6f3b3b99 |
| 可执行文件 | /home/fortneu/lr-grad/abacus-develop/build/abacus_std_para |
| 电子结构 | esolver_type=ksdft，PBE，LCAO DZP，91 AO / 91 bands |
| 自旋与电荷 | nspin=2，nelec_delta=1，nupdown=2；α 16 / β 14 |
| k 点 | 全部 8 个 (±¼,±¼,±¼)，Direct，每点权重 1/8 |
| 对称性开关 | symmetry=-1，保留完整 k 网格；后处理独立找群 |
| SCF | ecutwfc=100 Ry，scf_thr=1e-10，Gaussian σ=1e-4 Ry |
| 输出 | out_wfc_lcao=1，out_hsr=1 12，out_stru=1，out_mul=1，out_app_flag=0 |
| 运行 | 单进程，OMP_NUM_THREADS=1；24 次迭代，计时表 113.19 s |
| 收敛 | 最后 DRHO=5.0903e-11；总能量 −1181.662747465605 eV |

INPUT / KPT / STRU、程序和生成脚本 SHA256 记录在计算目录的 provenance.json；
完整分析和所有输入输出 SHA256 在 stateid-report.json。
仓库内 [回归数据](../../tests/fixtures/nv_non_gamma/README.md) 包含未改写的原始文本：
S(R)、β 自旋前三个 k 点完整波函数使用 gzip 压缩，其他元数据保留原文。
测试先解压到临时目录，再调用正式解析器，不依赖本机 ABACUS 安装。

初次输入检查发现该版本要求 ksdft 而非 ks，已修正后重新启动；
原始输入检查日志保留为 input_validation.log。成功计算的 INPUT 与校验和一致。

## 2. 解析器与单位核对

新增：

- read_structure：单帧显式晶格 STRU / STRU_FINAL；Direct、Cartesian、
  Cartesian_au、Cartesian_angstrom；保留原子顺序与未解释的 m/v/f/mag 字段。
- read_kpoints_info：读取 KPT.info 中实际采用的 Direct 坐标、编号和原始权重。
- read_input_parameters：读取 INPUT / INPUT.info 的原始参数字符串。
- read_scf_log / AbacusReader.read_output：单次 SCF 的版本、commit、
  收敛标志、迭代次数、最终能量。
- read_multik_wavefunctions：在显式 nspin 和离子步下读取新版
  wfk{k}s{spin}g{step}_nao.txt，复用已有系数 parser；校验头部全局 k 编号与坐标。

本次 β 的第一个空间 k 文件名为 wfk1s2g1_nao.txt，
文件头的 index of k points 却是 9：先列 8 个 α，再列 8 个 β。
这两个编号不可混用。occupations 为 ABACUS 原始 wg，
本次占据带每项为 0.125；对 k 和 bands 求和后分别为 16、14。
不能把 0.125 当作每个轨道的真实分数占据，也不能再重复乘 k 权重。

波函数头给的是以 2π/lat0 为单位的笛卡尔坐标，转换为倒格分数坐标：

$$
\mathbf k_{\rm frac}=\mathbf k_{\rm cart}A_{\rm input}^{T},
$$

其中 A_input 是**原始输入 STRU 的无量纲行晶格矩阵**。
STRU_FINAL 可改用新的 lat0 表示同一个物理晶胞，因此不能直接用其无量纲矩阵
转换旧 kvec_c。本次输入/输出物理晶胞差小于 5×10⁻¹¹ Å。

还修正了旧 parser 的字段单位：新版 CSR 的 UcellIO::write_ucell 输出
lat0 × BOHR_TO_A，因此字段现在叫 lattice_constant_angstrom，
不是 lattice_constant_bohr。这是早期 API 的更正；
STRU 的 lattice_constant_bohr 仍确实是 Bohr。
CSR 晶胞头只有约六位有效数字，本次物理晶胞误差为 2.89×10⁻⁶ Å，
只适合粗核对，不用它反算高精度 k。
另外，标有 gamma-only folded matrix 的 CSR 只含 X(Γ)，现在会拒绝非 Γ Fourier；
本次明确 gamma_only=0，保留了完整 93 个 R 块。

## 3. 联合验证的公式

由同次 S(R)、H(R) 构造：

$$
S(k)=\sum_R e^{+2\pi i k\cdot R}S(R),\qquad
H(k)=\sum_R e^{+2\pi i k\cdot R}H(R).
$$

先检查

$$
C_k^\dagger S(k)C_k\simeq I,\qquad
H(k)C_k\simeq S(k)C_k\varepsilon_k.
$$

对操作 g=(W,t)，k'=W^{-T}k，原子映射满足

$$
W\tau_a+t=\tau_b+L_{g,a},\qquad
B_{b\mu,a\nu}(g,k)=T_{b\mu,a\nu}(g)e^{-2\pi i k'\cdot L_{g,a}}.
$$

在小群中，k' 与 k 等价；cell gauge 下同一套 AO 基周期相同。
对候选 band 子空间 C，project_representation 计算

$$
Q=C(C^\dagger SC)^{-1/2},\quad
D(g)=Q^\dagger S(k)B(g,k)Q,\quad \chi(g)=\operatorname{Tr}D(g).
$$

它同时检查 B†SB=S 和 BQ≈QD；之后 analyze_little_group
检查群乘法并与 spgrep 的完整不可约特征标分解。
本例仅用能量差 10⁻⁶ Ry 产生候选簇，再逐个检查这些独立条件；
能量接近本身不是识别结论。

跨 k 点则使用两个独立的子空间：

$$
D_{k'\leftarrow k}=Q_{k'}^\dagger S(k')B(g,k)Q_k,\qquad
B^\dagger S(k')B=S(k).
$$

这里检查跨 k 闭合和幺正性，不把其 trace 当作特征标。
还独立检查 B†H(k')B=H(k)，以验证几何操作确实保持本次自洽哈密顿量。

## 4. 本次结果

两个 [111] 轴上 k 点的小群阶为 6（C₃ᵥ），其余六点的小群阶为 2（Cₛ）。
两自旋、每 k 全部 91 带，共有 1,336 个候选子空间通过闭合及分解；
检查了 48 个跨 k 映射，跨 k 子空间取各自前 16 带。

| 全数据最大残差 | 数值 |
|---|---:|
| max\|C†SC−I\|（逐元素） | 2.29×10⁻⁸ |
| max\|HC−SCε\|（逐元素，Ry） | 2.08×10⁻⁷ |
| AO metric 相对残差（白化基） | 6.51×10⁻¹⁰ |
| 子空间闭合误差 ‖BQ−QD‖S/√d | 3.68×10⁻⁸ |
| 群乘法残差 | 6.52×10⁻¹⁴ |
| 特征标分解重构误差 | 6.52×10⁻¹⁴ |
| 跨 k 子空间闭合误差 | 7.76×10⁻¹⁰ |
| max\|B†H(k')B−H(k)\|（逐元素，Ry） | 4.07×10⁻¹⁰ |

波函数文本只有约 8 位精度，所以不要求原始 C†SC 达到机器精度。
特征标可以因对称结构和 Q 正交化更加精确，但这不能替代闭合残差。

在 k=(−¼,−¼,−¼)，选几个 β 自旋轨道，以文件中的一基 band 编号表示：

| bands | 维数 | χ(E), χ(C₃), χ(σᵥ) | 空间 irrep |
|---|---:|---|---|
| 10 | 1 | [1, 1, −1] | A₂ |
| 14 | 1 | [1, 1, 1] | A₁ |
| 15、16 | 2 | [2, −1, 0] | E |

在 15/16 内作随机复酉混合，特征标不变；只取其中一个轨道会被闭合检查拒绝。
这里用的是**零平移代表元下的点群 Mulliken 名称**，不是擅自赋予某个 BCS
k-line 编号。c3v_labels 现在也接受这一有限 k 的情形；
带平移的有限 k 代表元仍须用包含相位的标签表。
它们是 KS 轨道的空间标签，尚未据此认定某个 LR 根为 ³E；
nupdown=2 固定的是自旋投影，不是独立的 S² 验证。

## 5. Bloch 相位正反对照

此 NV 构型在当前原点下，六个操作的原子返回矢量都为零。
仅在这种放置下观察 B=T，不能验证相位实现正确。

因此从同一批真实 C、S 出发，将第一个 C 原子的代表改为 τ₀+(1,0,0)。
这没有改变物理结构，只改变 AO 的晶胞代表。
若每个 AO 所属原子的移动为 Lμ，定义

$$
V_{\mu\mu}(k)=e^{+2\pi i k\cdot L_\mu},\qquad
C'_k=V(k)C_k,\quad S'(k)=V(k)S(k)V(k)^\dagger.
$$

新的 Bloch 操作必须满足

$$
B'(g,k)=V(k')B(g,k)V(k)^\dagger.
$$

用新原子代表重新调用 builder，而不是手工塞入这个预期 B'，
再投影 β 的 15/16 E 对：

- 特征标变化：8.90×10⁻¹⁶；闭合误差：2.20×10⁻⁹。
- 故意使用不含相位的 T：闭合误差 0.620，白化后的 AO metric 误差 81.1。

这项回归验证相位的符号、施加到源 AO 列的位置，以及 C/S/B 的规范一致性。
另外还有错误 Fourier 符号负例；不能只通过某个实数 character 反推所有约定都正确。

## 6. 如何复现后处理

在 stateid 仓库中：

~~~bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python \
  examples/abacus_non_gamma.py \
  "$HOME/work/LR-TDDFT-cal/grad/18_NV_diamond/111/pbe/stateid_non_gamma_20261007" \
  --output /tmp/nv-non-gamma-report.json

OPENBLAS_NUM_THREADS=1 .venv/bin/python -m unittest discover -s tests -v
~~~

示例需要 stateid[symmetry]（spglib/spgrep），不需要 IrRep 标签依赖。
完整测试：106 项通过，包含本次真实数据和原有合成不变量测试。
没有重新运行 ABACUS 的环境也可直接运行仓库回归。
如需重新计算，可复制 fixture 中 INPUT/STRU/KPT 到新目录，
修改赝势/轨道路径，并用记录的 ABACUS 版本运行；
本次使用 Intel oneAPI 运行库，source /opt/intel/oneapi/setvars.sh 后设置 OMP_NUM_THREADS=1。

下一步优先接入原生 LR/TDA 的全局 ph 索引、spin/k/MPI 排列及振幅，
并把轨道表示、参考态表示与完整 LR 根子空间连起来。
几何磁性信息虽保留在原始文本，当前流程仅验证标量共线计算的幺正空间操作；
不处理 SOC、反幺正磁群，也不自动换胞。
