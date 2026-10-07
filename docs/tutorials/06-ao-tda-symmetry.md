# 06：显式 Γ AO 操作与矩阵无关 TDA 分析

这一步补上交接文档中的 AO 操作构造和 TDA 根投影数学接口。
真实 NV 数据仍缺同几何 C/S/Orbital，原生 MPI amplitudes 仍需布局 manifest。
本教程不宣称已经验证真实 NV 的 a1/e 或 ³E。

## ABACUS 实操：采集同一几何的元数据

以下组合已在本地 Γ、共线双自旋 LCAO/TDA 算例验证。加入 INPUT：

```text
out_wfc_lcao  1
out_mat_hs2   1 12
out_stru      1
out_freq_ion  1
out_app_flag  0
out_mul       1
out_wfc_lr    1
```

`out_mat_hs2` 是兼容名称；支持新名称的版本可改为 `out_hsr 1 12`，
二者选一个。这里 12 是 H/S 文本精度；`out_wfc_lcao 1` 的文本系数在
本地版本仍固定为小数点后 8 位，`out_ndigits` 不会提高其精度。

| 开关 | 所需数据及用途 |
|---|---|
| `out_wfc_lcao 1` | C↑/C↓，轨道在 AO 基中的列系数 |
| `out_mat_hs2 1 12` | H(R)/S(R)；Γ-only 的 S 用作 AO metric |
| `out_stru 1` | 与电子数据对应的几何及晶胞 |
| `out_freq_ion 1`、`out_app_flag 0` | 每个离子步分文件，避免把多帧当作单帧读取 |
| `out_mul 1` | 触发 `Orbital` AO 标签输出；不是为 irrep 使用 Mulliken 数值 |
| `out_wfc_lr 1` | 完整 TDA 振幅；MPI 文件仍需要可信的布局 manifest |

本地版本的第一个已计算离子步对应：
`WFC/wfs1g1_nao.txt`、`WFC/wfs2g1_nao.txt`、`srg1_nao.csr`、
`STRU1`，LR 振幅文件则含 `step0`。因此 `gN`/`STRUN` 对应 LR
`step(N-1)`，不要按字符串中的数字直接配对。`Orbital` 首列为重复的
零基原子编号；全局 AO 顺序是文件行顺序。

这些参数**不是只输出最后一帧**：`out_freq_ion 0` 对 LCAO C/S 仍会
在各电子计算收敛或达到 SCF 步数上限时输出；`out_app_flag 0` 使用
分步文件，`1` 追加。当前负数会被重置为 0，不能用 `-1` 请求末帧。
应从一个干净的输出目录采集，并选择最后一个具有完整 C/S/X 的已计算帧。
尚未收敛或达到离子步数上限时，退出结构可能已经移动而未重新计算；
不能默认 `STRU_FINAL` 与最后的 C/S/X 同几何。归档时同时保留对应编号
STRU、步号、SCF 收敛状态及优化停止原因。

静态 `scf` 可先采集基态 C/S/Orbital，无需为补元数据重跑整条 JT 轨迹。
分析前检查 C†SC、AO metric、轨道和根子空间闭合与群关系；完整 MPI
振幅不能按 rank 文件后缀当 root，也不能未经布局验证直接拼接。
Γ-only 的 S 不能用于任意 k 插值。上述文件采集并不意味着结构解析、
自动群发现或 MPI adapter 已由 stateid 自动完成。

## AO 行标签与主动操作

```python
from stateid.io import read_ao_labels
from stateid.symmetry import build_gamma_ao_operation, analyze_c3v

labels = read_ao_labels("OUT/Orbital")
built = build_gamma_ao_operation(
    fractional_positions, lattice_vectors, species, labels,
    cartesian_rotation, translation=cartesian_translation, symprec=1e-6)
T = built.matrix
```

晶胞是 `(3,3)` 行向量，长度单位由调用者统一提供，分数坐标为 `(nat,3)`。
`rotation` 是实正交矩阵，对笛卡尔列向量作主动旋转；`translation` 是同单位
笛卡尔平移，旋转原点是笛卡尔零点。绕其他原点 o 时传入 `translation=o-R@o`。
镜面无需改写成 proper rotation，s/p/d 多项式会包含相应宇称。
当前显式拒绝 ell>2，不支持 k≠Γ、SOC 或 spinor。

`AOLabel(atom_index, species, ell, zeta, m)` 的原子编号零基、zeta 一基、m 有符号。
`read_ao_labels` 把 Orbital 文件的 component index 0..2l 转为 `(0,+1,-1,...)`。
绝不能把文件首列的重复原子编号当作全局 AO index。
每个原子/ell/zeta 必须有完整 m 壳层；AO 行顺序可任意，但必须与 C/S 完全一致。
映射只能在相同 species 和相同径向壳层之间，禁止猜测每原子 AO 数。

返回值记录每个源原子到目标原子的映射、整数返回矢量和笛卡尔残差：
`R @ source + translation = target + lattice_shift @ lattice`（位置改写为列/行时须一致）。
实际实现的行向量等式是 `source_cart @ R.T + translation = target_cart + shift @ lattice`。
匹配须唯一且双射，旋转须保持周期晶格。`symprec` 是匹配阈值，不是强行恢复
对称性的开关；报告 `atom_errors`，再检查电子参考态是否保持该群。

将六个 `built.matrix` 按 `E/C3/C3^2/sv1/sv2/sv3` 传给 `analyze_c3v`，
由该函数检查 T†ST、所选轨道空间闭合与投影群关系。构造器本身没有 S，
因此不会替你证明电子解或完整 AO 度量的对称性。结构解析和自动群发现仍待实现。

## TDA 总空间表示

```python
from stateid.symmetry import analyze_tda_c3v

result = analyze_tda_c3v(
    X_pair,  # (nph,2)，选择完整根对
    ph_operations,  # 每操作 [(D_occ_up,D_virt_up),(D_occ_down,D_virt_down)]
    reference_occupied_operations=reference_operations,
    # 每操作 [D_all_occ_up,D_all_occ_down]
)
print(result.symmetry.characters, result.symmetry.match.irrep)
print(result.root_orthonormality_error, result.reference_characters)
```

每个 spin 块占据指标在外、虚轨道指标最快，X 行顺序是各完整块依次连接。
函数按每根 `D_occ.conj() @ X.reshape(nocc,nvirt) @ D_virt.T` 收缩，再乘
所有块共用的 `det(D_all_occ_up)*det(D_all_occ_down)` 参考相位。
内存不随 nph² 增长；需要存储轨道表示矩阵、X 和同形的变换后根。
输入根必须欧氏正交归一，函数不会静默重归一阈值裁剪表。
单独一个 E 分量会因泄漏被拒绝，完整根对在随机 U(2) 混合后 character 不变。

参考 occupied 矩阵必须来自**完整整数占据**的单行列式空间，不能以 LR occupied
窗口的行列式代替；激发窗口与完整参考占据空间是独立参数。
空自旋通道可用 `(0,0)` 矩阵，行列式为 1。所有非空轨道表示都会检查酉性及
完整 C3v 群关系。但本函数只收到 D，不收到原始 C/T/S；调用者必须先通过
`analyze_c3v` 的闭合与度量检验，不能把泄漏的轨道空间传入并宣称已验证。

这返回的是包含参考相位的**总空间 irrep**，不是自旋分类。
不支持分数占据、多行列式参考、spin-flip、多 k、稀疏 ph 表或完整 LR X/Y。
输入 `X_pair` 应来自可信的完整 TDA 数据；原生 ABACUS MPI rank 文件
尚未自动恢复，不能把 rank 后缀当 root 或直接拼接。

测试 `tests/test_ao.py` 检查非交换旋转、直接 p/d 函数点值、镜面与 inversion、
双径向壳层、周期返回矢量及非正交 S。`tests/test_tda_symmetry.py` 用小型 dense
kron 对照双 spin 矩阵无关结果，检查复 gauge、A2 参考相位以及错误路径。
