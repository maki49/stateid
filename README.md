# stateid

**从第一性原理输出识别电子态的物理性质，并解释识别依据。**

`stateid` 是 electronic **state identification** 的简写。相比 `state-identifier`，
它更短，可以直接 `import stateid`；相比 `electronic-state-tools`，它更明确地强调
“这个态是什么、证据有多强”，而非通用性质计算；相比 `state-analyzer`，它突出最终
的物理归属。名字不绑定 ABACUS、基态或某一类缺陷，方便后续扩展。
这是本地项目名称，尚未发布到 PyPI，也未确认公共名称的可注册性。

## 项目定位与当前边界

愿景是把能量、轨道/激发态对称性、自旋、跃迁组分和局域性等证据组织起来，
辅助判断缺陷态、分子态及固体中局域电子态的物理身份。教程同时讲公式、适用条件
和诊断量，而不只是文件操作。首阶段 **ABACUS-first**，数值核心与软件 IO 分离。

当前为 `0.1.0` 初始框架：基态轨道与 Slater 行列式可分析；LR/TDA 已有数据模型、
标准化交换格式、基础对称性和显式算符收缩接口，**还不能直接从完整 ABACUS LR
计算目录自动识别激发态，也不能自动输出 NV⁻ 的 ³E 结论**。

| 模块 | 已实现 | 待实现 |
|---|---|---|
| ABACUS 波函数 IO | 单帧、Γ 点、实数 LCAO 文本 | 复数多 k、二进制、追加帧、spinor |
| ABACUS 元数据 | `OutputReader` / `AbacusReader` 接口 | 日志、STRU、AO 标签与计算版本联动 |
| 重叠矩阵 IO | 显式 `stateid_npy` 稠密 S | 原生新旧 CSR S(R)、S(k) |
| LR eigenvector IO | `stateid_npz` v1 的 X/Y、能量、ph 映射 | ABACUS 原生 LR 文件、MPI 分片重组 |
| symmetry | S-正交化、D/χ、闭合误差、C3v 群关系及 irrep 匹配 | 从结构与 AO 标签自动组装 T(R) |
| 球谐约定 | ABACUS m 顺序与复/实球谐基变换 | 完整 Wigner D、Euler 角和原子映射适配 |
| spin | unrestricted 行列式 ⟨S²⟩；给定 S²_ph 的 TDA 收缩 | 自动构建 S²_ph；完整 LR 的自旋响应方法 |

未实现的格式/API 会明确抛出 `NotImplementedError`（IO 使用其子类
`UnsupportedFormatError`），不会猜测矩阵顺序或用虚构数据继续计算。
初期限于无 SOC、共线自旋、Γ 点/有限体系的幺正空间对称性；不支持磁群、双群、
反幺正操作或连接不同 k 点的操作。数组核心允许复数，这不等于已经支持这些物理情形。

## 安装和运行

Python ≥ 3.10，运行时仅依赖 NumPy；测试使用标准库 unittest。
从仓库根目录执行：

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -e .
python -m unittest discover -s tests -v
python examples/c3v_minimal.py
```

本次 WSL 初始化已创建 `.venv`，复用本机已有 NumPy；可直接从激活环境这一步开始。
纯净环境按上面的安装流程重建即可。预期最小示例输出：

```text
A1: characters=[1. 1. 1.], identified=A1
E: characters=[ 2. -1.  0.], identified=E
two parallel spins: <S^2>=2.0
```

也可直接匹配按共轭类排列的数值特征标：

```python
from stateid.symmetry import match_characters
assert match_characters([1, 1, 1]).irrep == "A1"
assert match_characters([2, -1, 0]).irrep == "E"
assert not match_characters([2, -0.7, 0]).valid
```

这里的三个数是**各类代表元的特征标或类内平均**，不是类内求和。
数组匹配只能检查特征标相容性；完整流程应使用 `analyze_c3v`，同时检查子空间闭合和群关系。

## 物理核心：表示矩阵与 character

轨道按列存储：\( |\psi_n\rangle=\sum_\mu |\phi_\mu\rangle C_{\mu n} \)，
AO 重叠为 \(S_{\mu\nu}=\langle\phi_\mu|\phi_\nu\rangle\)。定义

\[
\hat R|\phi_\nu\rangle=\sum_\mu|\phi_\mu\rangle T_{\mu\nu}(R),\qquad
M_{\mu\nu}(R)=\langle\phi_\mu|\hat R|\phi_\nu\rangle=(ST)_{\mu\nu}.
\]

若 \(C^\dagger SC=I\)，则

\[
\boxed{D(R)=C^\dagger M(R)C=C^\dagger ST(R)C},\qquad
\boxed{\chi(R)=\operatorname{Tr}D(R)}.
\]

代码必须显式选择 `operator_kind="coefficient"`（输入 T）或 `"matrix_element"`
（输入 M）。**本库的 M 定义不等同于 ABACUS 源码中所有名叫 M 的对象。**

对尚未正交归一的子空间，计算 \(G=C^\dagger SC\)，再构造
\(Q=CG^{-1/2}\)，用 Q 代替 C。若 G 奇异则要求重新选择子空间，不静默丢弃轨道。
投影算符的 AO 系数表示为 \(P=QQ^\dagger S\)，闭合误差由

\[
\epsilon_R=\frac{\|T(R)Q-QD(R)\|_S}{\sqrt{d}},\quad
\|A\|_S^2=\operatorname{Tr}(A^\dagger SA)
\]

衡量。只有小误差的完整不变子空间才适合赋予精确 irrep 标签。
实际缺陷结构若有应变或 Jahn–Teller 畸变，可能不再具有精确 C3v 对称性，不能只放宽容差掩盖它。

### 为什么不能只看单根简并轨道

E 是二维不可约表示。对简并对 \(C_E=(c_1,c_2)\)，数值求解器可以返回任何
\(C'_E=C_EU\)，其中 U 为二维酉矩阵。于是

\[
D'(R)=U^\dagger D(R)U,\qquad \operatorname{Tr}D'(R)=\operatorname{Tr}D(R).
\]

单根轨道的形状、对角矩阵元和“x/y 方向”依赖所选基；完整子空间的 character 不依赖它。
只取一个实 E 分量，旋转后通常落到另一个分量，会产生明显闭合误差。

| C3v | E（恒等） | 2C3 | 3σv |
|---|---:|---:|---:|
| A1 | 1 | 1 | 1 |
| A2 | 1 | 1 | −1 |
| E（二维 irrep） | 2 | −1 | 0 |

重数采用 \(n_\alpha=\frac16\sum_c |c|\chi_\alpha(c)^*\chi(c)\)。
程序检查其非负整数性和重构误差；可以识别 A1+E 等可约表示，不强行归为单个 irrep。

## 自旋：⟨S²⟩ 能回答什么

对于整数占据、共线自旋的 unrestricted Slater 行列式：

\[
\frac{\langle\hat S^2\rangle}{\hbar^2}
=M_S^2+\frac{N_\alpha+N_\beta}{2}
-\sum_{ij}|(C_\alpha^\dagger S C_\beta)_{ij}|^2,
\quad M_S=\frac{N_\alpha-N_\beta}{2}.
\]

单重态、双重态、三重态的纯自旋值分别为 0、3/4、2。接口接收两个自旋通道的
**全部已占据空间轨道**；只读一个 KS 轨道或只看 `nspin=2` 不能确认三重态。
UKS 中该式是辅助行列式的诊断，不是精确相互作用多电子态的自旋测量。
期望值本身也不能普遍证明自旋纯度，混合态可能具有相同的平均值。

TDA 的预留路线是，在指定的正交行列式/CSF 基底中构造完整算符后计算
\(X^\dagger S^2_{ph}X/(X^\dagger X)\)。现有 `tda_s2` 只负责这个收缩；
`build_s2_ph` 尚待实现。完整 LR 的 X/Y 具有响应理论的归一化，不能直接套用
普通 CI 向量公式，`lr_s2` 暂时明确拒绝计算。推导见[自旋教程](docs/tutorials/02-spin.md)。

## NV⁻：未来使用示意

目标证据链是 ↓126 的 A1 轨道、↓127/128 的 E 子空间、LR 中 a1→e 的跃迁组分，
再结合参考态、多电子对称性和自旋证据讨论 ³E。轨道 irrep 用小写 a1/e 表述时，
程序仍统一返回表中的 `A1`/`E`，不自动添加自旋上标。

```python
from stateid.io import AbacusReader
from stateid.symmetry import analyze_c3v

# 示意：替换成真实文件；S 和 operations 尚需人工提供/后续适配器生成。
reader = AbacusReader()
down = reader.read_wavefunctions("path/to/down_gamma.txt", spin="beta")
S = reader.read_overlap("path/to/validated_S.npy", format="stateid_npy")
# operations: 完整六个 AO 系数变换 T(R)，需先验证 AO 顺序与旋转约定。
# 此处假设 126/127/128 是文件中从 1 开始的 band 标签；先核对原始输出！
a1 = analyze_c3v(down.coefficients[:, [125]], operations, S,
                 operator_kind="coefficient")
e = analyze_c3v(down.coefficients[:, [126, 127]], operations, S,
                operator_kind="coefficient")
```

上面是未来真实工作流示意，`operations` 未自动生成，不是现在可直接跑的 NV⁻ 示例。
可运行的合成示例在 `examples/c3v_minimal.py`；更完整的证据边界见
[NV⁻ 教程](docs/tutorials/03-nv-minus.md)。

## 路线图与目录

1. **已完成的种子层**：数值核心、C3v、行列式自旋、Gamma 文本、数组交换与合成测试。
2. **下一步优先**：实现 ABACUS 原生 `sr_nao.csr` / 旧式 S(R) 读取，显式选择离子步，
   构造 \(S(\Gamma)=\sum_R S(R)\)，与同一次计算的 C 校验 \(C^\dagger SC\approx I\)。
3. 接入 STRU/数值轨道头信息与 AO 标签，从结构、原子映射、实球谐旋转组装 T(R)。
4. 接入 ABACUS LR 全局 ph 索引、spin block、MPI 分片、X/Y 定义；先 TDA，再完整 LR。
5. 增加真实 NV⁻ 回归案例、更多点群与多电子态分析；稳定数据模型后添加其他软件适配器。

```text
stateid/
├── README.md
├── pyproject.toml
├── src/stateid/
│   ├── io/          # 软件适配器、数据模型、明确的格式边界
│   ├── symmetry/    # 表示、特征标、C3v、球谐与 ph 基础操作
│   └── spin/        # 行列式 S²、TDA 收缩与 LR 预留 API
├── examples/c3v_minimal.py
├── tests/          # 数学不变量、错误输入、IO 格式测试
└── docs/
    ├── abacus-conventions.md
    ├── io-formats.md
    └── tutorials/  # 推导、数值诊断与物理证据链
```

本地 ABACUS 参考路径、commit 和符号约定见[约定说明](docs/abacus-conventions.md)。
开发没有改动 ABACUS 源码，也没有把其 C++ 实现复制进本项目。
贡献新后端时实现 `OutputReader` 协议，不让文件布局渗入物理核心。
项目尚未选择发布许可证；公开分发前由维护者确定。
