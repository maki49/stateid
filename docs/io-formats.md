# IO 合约与支持范围（0.2.0）

## 统一约定

- C 的形状 `(nao, nband)`；所有波函数按列存储。不同 k 文件分别返回一个 `OrbitalData`。
- 文件 `band_numbers` / `k_index` 保留一基标签，NumPy 下标从 0 开始。
- 波函数能量保留 Ry，交换 LR 能量为 eV。CSR 数值保留输入单位，不自动转换 H 的 Ry。
- `occupations` 保留原始字段，可能包含 spin/k 权重，不自动解释为单行列式整数占据。
- C、S、T 必须对应相同结构、离子步、k 点、AO 基与顺序。文件名不承担这些校验。
- alpha/beta 显式传入；不根据文件名推断 spin，也不自动解释 SOC spinor。

## Γ 与复数 k 点波函数：共用一个解析器

`read_wavefunctions(path, format="auto", spin="beta")` 按头自动判断；
也可指定 `gamma_text` / `k_text`。`read_gamma_wavefunctions` 保留为严格 Γ 包装器。
读取 `wfc_nao_write2file` / `wfc_nao_write2file_complex` 的单帧文本。

共同部分是 nbands、nao 与逐 band 的 `(band)`、`(Ry)`、`(Occupations)`。
Γ 系数为实数；复数 k 文件在共同部分之前有：

```text
1 (index of k points)
0.25 0.1 -0.2
```

复数文件每个 AO 系数是相邻的 `real imag` 两个数，可跨行。
系数与维数、连续 band 编号、有限数值均检查；支持 e/E/d/D 指数。
二进制、一个文件内追加多波函数帧、分布式 shard 和 spinor 语义尚不支持。

`k_cartesian` 保存头中的 ABACUS `kvec_c`，单位 **2π/lat0**。
用 `stateid.realspace.k_cartesian_to_fractional(k_cartesian, lattice_vectors)`
转成 Fourier 所需的倒格分数坐标。lattice_vectors 是以 lat0 为单位的三个直接晶格
行向量；不是带 Å 单位的完整晶胞。新 CSR 中可以取得它，旧 CSR 需要外部结构数据。

## 原生新旧文本 CSR：统一解析

`read_csr(path, frame=None)` 自动识别：

- 旧/LTS：可选 `STEP: ...`，`Matrix Dimension of ...` / `Matrix number of ...`。
- 新版：`--- Ionic Step ... ---`，spin/维数/R 数的注释字段、结构块和 `CSR Format` 标记。

两种文件头共用同一个 R 数据块解析器，支持实数和 C++ `(real,imag)`（允许括号内空格）、
跨行数组、注释、空 R 块、新版空块的全零行指针。
检查非负 nnz、R 唯一性、列索引范围、行指针起点/终点/单调性、有限数值和截断。
CSR 重复的同一矩阵位置由 SciPy 按标准语义求和；重复 R 块会报错。

返回 `RealSpaceMatrix`，用 SciPy CSR 储存 `(nR, nao*nao)` 的完整数值，
AO 对索引按 C 顺序展平。保留原始 `ionic_step`、格式类型、spin_index（若有）、
新版 `lattice_vectors` 与 `lattice_constant_bohr`。

多段/多离子步文件必须显式传 `frame=0,1,...`；这是文件段序号，**不是原始 ionic_step**。
`iter_csr(path)` 可以逐段读取，避免同时保存所有步。
二进制、缺少所述头/新版 CSR 标记的其他方言尚不支持。

## 任意 k 的 Fourier 与重叠矩阵检查

\[
X(k)=\sum_R e^{+2\pi i k\cdot R}X(R).
\]

`data.to_k(kpoints)`：输入 `(3,)` 返回 `(nao,nao)`；输入 `(nk,3)` 返回 `(nk,nao,nao)`。
使用倒格分数坐标，不自动加入原子中心 τ 相位、不除以 nR、不乘 k 权重。
完整 X(R) 保持稀疏，结果为稠密；大 k 网格请分批调用。
每个 X(R) 不必 Hermitian，接口不会补齐下三角或缺失的 −R 块。

`AbacusReader.read_overlap(path, format="abacus_csr", kpoints=..., frame=...)`
复用同一读取/Fourier 路径，在每个 k 检查 Hermitian 与正定。
省略 kpoints 表示 Γ。`atol` 默认 1e-8，传给已有度量校验；
需结合输出精度调整，同时留意它也用于判定近线性依赖，不能任意放大。
一般 H 不应使用 read_overlap，而应走 read_csr→to_k。

`stateid_npy` 仍支持已构造的单个 S(k)，不允许再传 kpoints/frame。

## LR：pickle-free NPZ schema_version=1

```python
import numpy as np
np.savez("lr.npz", schema_version=1, method="TDA", X=np.eye(2, dtype=complex),
         energies_ev=np.array([2.1, 2.1]),
         transitions=np.array([[0, 1, 125, 1, 126], [0, 1, 125, 1, 127]]))
```

X 的形状 `(nph,nroot)`；transitions 为 `(k, spin_hole, occupied_band, spin_particle, virtual_band)`。
索引从 0 开始，spin 0=alpha、1=beta，band 为全局 KS 列号；每行唯一。
TDA 不含 Y，`full_lr` 必须含同形状 Y；不自动归一化，也未实现完整 LR 自旋计算。
导出者须负责全局 ph 映射，不能简单重命名 ABACUS 的 rank 局部振幅文件。

## 仍为占位的入口

`read_output`、`read_lr_eigenvectors(..., format="abacus_lr")` 明确报错。
多 k 小群/空间群操作、跨 k 变换、spinor、自旋翻转分析和完整 LR 自旋仍需后续工作。
`OutputReader` 为其他软件后端保留相同接口。
