# 04：从 S(R) 到多个 k 点的轨道内积

## 为什么 Γ 点必须求和，不能只取 R=0

在 cell gauge 中，Bloch AO 的重叠为
\(S(k)=\sum_R e^{+2\pi i k\cdot R}S(R)\)。因此 Γ 点是所有晶格平移块之和。
R 是整数直接晶格坐标，k 是倒格分数坐标；二者点积无量纲。
这个正号与本地 pyATB `base_data::get_exp_ikR` 一致。

Hermitian 条件是 \(S(-R)=S(R)^\dagger\)，不要求每个 R 块都 Hermitian。
合成测试刻意采用复数且单个 R 块非对称的模型，保证能发现相位符号、转置和
提前补齐三角矩阵的错误；仅用实数 Γ 点测试发现不了这些问题。

## 一次读取，多次计算

```python
import numpy as np
from stateid.io import read_csr, read_wavefunctions
from stateid.realspace import k_cartesian_to_fractional

sr = read_csr("path/to/sr_nao.csr", frame=0)
wf = [read_wavefunctions(path, spin="alpha") for path in
      ["path/to/wfs1k1_nao.txt", "path/to/wfs1k2_nao.txt"]]

# 新版 CSR 提供 dimensionless lattice_vectors；旧版需从同次计算的结构中提供。
kfrac = k_cartesian_to_fractional([w.k_cartesian for w in wf], sr.lattice_vectors)
sk = sr.to_k(kfrac)
for w, s in zip(wf, sk):
    c = w.coefficients
    error = np.linalg.norm(c.conj().T @ s @ c - np.eye(c.shape[1]))
    print(w.k_index, error)
```

以上路径为示意。一定要使用匹配的同一次计算、同一离子步和 AO 顺序。
新文件头中的 k 是笛卡尔坐标，并且可能带 k 点权重的 occupation；
它们分别不能直接当作分数 k 坐标和每轨道电子数。
应按物理 k 坐标关联数据，不能只依据 filename/index（nspin=2 的索引可能覆盖自旋块）。

使用与软件无关的 `RealSpaceMatrix` 可在同一 R 列表上计算 H(k)、S(k) 等一般矩阵。
不同物理矩阵不必强行具有相同 R 列表，各自 Fourier 即可。
大网格可以切片 kpoints 后分批调用，避免最终 `(nk,nao,nao)` 稠密数组过大。

## 已执行的验证（2026-10-06）

自动化测试覆盖新旧格式等价、空块、离子步选择、复数值、错误 CSR、k→k+G 周期性、
正相位的解析模型、非正交晶胞坐标变换、复数波函数和原有 symmetry/spin 功能。

对本地 `/home/fortneu/abacus/abacus-develop` 的文件另做了只读交叉验证：

| 文件 | 格式 / 规模 | 对 pyATB Fourier 上三角的最大误差 |
|---|---|---:|
| `examples/10_hs_matrix/02_out_hsr_multik/sr_nao.csr` | 旧，26 AO / 183 R | 0（本次浮点结果） |
| `tests/03_NAO_multik/scf_out_hsr/sr_nao.csr.ref` | 新，8 AO / 13 R | 0（本次浮点结果） |
| `source/source_io/test/support/SR.csr` | 新，4 AO / 2 R | 0（本次浮点结果） |

比较 k=(0,0,0)、(0.173,−0.219,0.317)、(0.5,0.25,0.125)。
前两份物理矩阵的 Hermitian 残差分别约 2.22e−15、7.76e−18，三个 k 上的
最小本征值分别约 0.0026361、0.8567144，均为正。
第三份是非 Hermitian 的 parser 单元测试数据，只验证读取/Fourier，不能作为物理 S(k)。

pyATB 顶层导入受本机 MPI 动态库缺失影响；交叉验证直接执行其本地源文件中独立的
CSR reader 定义，不加载 MPI/C++ 包初始化。没有把它们复制到本项目或作为运行时依赖。
比较上三角是因为 pyATB 本身只保留这些项；stateid 的完整下三角另由解析模型和
前两份文件的 Hermitian 残差检查。

另读取了 `source/source_io/test/support/wfs1k1_nao.txt` 至 `wfs1k4_nao.txt`：
均得到 `(63,3)` 的复数 C、正确的 k 标签与坐标。这四份 C 与上表的 S 不是同一体系，
**未宣称已验证它们之间的 C†SC**；同源真实 C/S 联合回归是下一步。

## 多 k 支持与物理态分类的区别

任意 k 的文件读取与 S(k) 已支持。对称性归属则只能使用把 k 保持到模倒格矢意义下
不变的小群。一般空间群操作会把 k 送到不同的点，不能直接对单个 k 的轨道套 C3v 表。
还需正确的 AO 原子映射、平移/Bloch 相位以及必要时的小表示；这些自动构造工作尚未实现。
多 k 加权的周期行列式自旋也不能简单等同于逐 k 调用 `determinant_s2` 后加权。
