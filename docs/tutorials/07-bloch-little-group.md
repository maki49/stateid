# Bloch AO 变换、小群、一般特征标分解和标准标签

本教程对应四个连续步骤：构造 B(g,k)，筛选小群，用 spgrep 生成参考表示，
将实际波函数的表示分解并映射到标准标签。范围为**原胞中的标量、幺正空间对称性**。
基础数学不依赖外部群库；spgrep 和 IrRep 分别作为可选后端。

## 1. 固定坐标和 Bloch 基约定

晶格矩阵 A 的每一行是一个实空间基矢。分数坐标按列变换：

\[
g=\{W|t\},\quad x'=Wx+t,\quad
R_{\rm cart}=A^TWA^{-T},\quad t_{\rm cart}=A^Tt.
\]

倒空间同样使用分数坐标，且不把 2π 吸收到 k 中：

\[
\boxed{k'=W^{-T}k}.
\]

这里 W 是分数坐标矩阵。只有在正交笛卡尔坐标下才能把倒空间操作直接写成
同一个旋转 R；在倾斜晶胞中不能误用 Wk。程序保留未折回的 k'，便于审计。

取 cell gauge：

\[
|\phi_{\alpha m,k}\rangle
=\frac1{\sqrt N}\sum_n e^{2\pi i k\cdot n}|\phi_{\alpha m,n}\rangle.
\]

指数只含晶格矢量 n，不含原子位置 τ。它与本库
\(S(k)=\sum_n e^{+2\pi i k\cdot n}S(n)\) 一致。
ABACUS 波函数头的笛卡尔 k 坐标应先转换，不能直接送入本接口。

## 2. 为什么 Γ 点的 T 需要加相位

对源原子 α，周期原子映射给出

\[
W\tau_\alpha+t=\tau_\beta+L_{g\alpha},\qquad L_{g\alpha}\in\mathbb Z^3.
\]

如果同一壳层的局域实球谐旋转为 \(T^{(\ell)}_{m'm}\)，则

\[
\hat g|\phi_{\alpha m,n}\rangle
=\sum_{m'}T^{(\ell)}_{m'm}
|\phi_{\beta m',Wn+L_{g\alpha}}\rangle.
\]

在 Bloch 求和中令 \(n'=Wn+L_{g\alpha}\)，得到

\[
e^{2\pi i k\cdot n}
=e^{2\pi i k'\cdot n'}e^{-2\pi i k'\cdot L_{g\alpha}}.
\]

因此

\[
\boxed{
B_{\beta m',\alpha m}(g,k)
=e^{-2\pi i k'\cdot L_{g\alpha}}T^{(\ell)}_{m'm}.
}
\]

**相位乘在每个源 AO 列上**，按源原子的 L 取值；不能给整个矩阵随意乘一个相位。
Γ 点相位为 1，恢复原来的 AOOperation.matrix。

- `build_bloch_ao_operation` 接收分数坐标 W、t；复用既有 s/p/d AO builder。
- `build_gamma_ao_operation` 仍接收笛卡尔 R、t，旧接口保持不变。
- `bloch_from_ao_operation` 可缓存原子映射/角动量旋转，仅在各 k 重算相位；
  调用者必须提供与缓存操作一致的 W 和 AO→atom 索引。
- 返回的 `mapped_kpoint`、`ao_operation.atom_mapping`、`lattice_shifts`、
  `atom_errors` 可用于审计。AO 行顺序和各原子的坐标代表不能与 C/S 混用。

ABACUS `cal_Ms` 中同样存在返回晶格矢量相位；这里由上述主动变换定义独立推导，
没有复制其矩阵存储/密度矩阵恢复实现。stateid 的 B 是**系数变换**；
在 `project_representation` 中必须使用 `operator_kind="coefficient"`。
stateid 的矩阵元 M 为 S B，不能仅按 ABACUS 对象名字 M 判断输入类型。

## 3. 小群和非对称型空间群的乘法

\[
G_k=\{g:W^{-T}k-k\in\mathbb Z^3\}.
\]

`find_little_group` 返回所选操作在原列表中的 indices、倒格矢差、操作本身、
乘法表和相位因子。translations 的整数部分**不被删除**。

给定平移子群的代表元 \(g_a,g_b,g_c\)，可能有

\[
g_ag_b=\{I|L_{ab}\}g_c,\quad
L_{ab}=t_a+W_at_b-t_c.
\]

在固定 k 的 Bloch 空间中，正确关系为

\[
\boxed{
D_aD_b=\omega_k(a,b)D_c,\qquad
\omega_k(a,b)=e^{-2\pi i k\cdot L_{ab}}.
}
\]

所以有限代表元的矩阵不一定满足普通点群乘法。例如滑移
\(g=\{m_y|\frac12 a_x\}\) 满足 \(g^2=\{I|a_x\}\)，因此
\(D(g)^2=e^{-2\pi i k_x}I\)。在 \(k_x=1/2\) 处，滑移本征值为 ±i，
不能强制套镜面本征值 ±1。

输入要求是完整的操作群、每个旋转一个平移陪集代表元，定义在**原胞**中。
重复旋转（例如中心化晶胞的多个纯平移）、缺操作、不闭合的操作集会明确报错。
`symmetry_from_structure` 可以用 spglib 从显式结构发现几何操作，但不自动换胞；
若输入为非原胞，须连同 C、S、AO 映射一起正确处理，不能只把结构单独原胞化。
缺陷超胞若没有更小的平移对称性，本身就是该缺陷周期体系的原胞。

几何对称不自动保证自洽电子解对称。还需逐操作验证 metric 和态子空间。

## 4. 从实际系数得到表示

先在非正交 AO metric 下正交化：

\[
G=C_k^\dagger S_kC_k,\qquad Q_k=C_kG^{-1/2}.
\]

对小群操作，cell gauge 下 k+整数倒格矢与 k 使用相同 Bloch AO 基，因而

\[
D_k(g)=Q_k^\dagger S_kB(g,k)Q_k,\quad
\chi_k(g)=\operatorname{Tr}D_k(g).
\]

`analyze_little_group` 会验证：

1. \(B^\dagger S_kB=S_k\)，以正交 AO 坐标下的残差衡量；
2. \(BQ_k-Q_kD_k\) 的 S 范数，防止只选中简并壳层的一部分；
3. D 的幺正性，以及上述带 ω 的群乘法；
4. 数值 character 是否可由参考不可约表示按非负整数重数重构。

只改变简并子空间的酉基底 \(Q\to QU\)，会使 \(D\to U^\dagger DU\)，
但不改变 character。漏掉 B 的相位或把操作顺序弄错，不能靠放宽匹配容差补救。

对于跨 k 映射，使用另一个接口：

\[
D_{k'\leftarrow k}=Q_{k'}^\dagger S_{k'}B(g,k)Q_k,\qquad
B^\dagger S_{k'}B=S_k.
\]

`project_sewing_matrix(C_k, C_kprime, B, S_k, S_kprime)` 返回矩阵、闭合与
metric 误差。它不返回 character，因为两端独立换基时
\(D\to U_{k'}^\dagger DU_k\)，迹一般不再不变。同 k 投影复用相同数值核心。
这个接口提供未来多 k LR 的基础；当前没有因此自动实现多 k LR 根表示。

## 5. spgrep 后端与一般特征标分解

安装：

~~~bash
python -m pip install -e '.[symmetry]'
python examples/bloch_little_group.py
~~~

`spgrep_irreps(rotations, translations, kpoint)` 先筛选小群，再调用 spgrep，
返回 `IrrepSet`。每个矩阵数组为 `(ng, d, d)`，严格对齐 `reference.group`
中的操作顺序。spgrep 要求恒等代表元的平移为零，适配器仅在调用后端时将平移
约化到基本区间，随后乘回整数平移相位，保持调用者的原始 Seitz 代表元。
后端返回值也会经过群关系、character 正交性与维数平方和检查。
使用复数 small irreps；不擅自把共轭表示合并为实表示，也不额外加入时间反演。

参考特征标为 \(\chi_\alpha(g)=\operatorname{Tr}\Gamma_\alpha(g)\)，分解为

\[
n_\alpha=\frac1{N_g}\sum_g\chi_\alpha(g)^*\chi_{\rm num}(g).
\]

对于带同一因子系统的 projective irreps，此内积仍适用。
`decompose_characters` 内部按每个操作处理，不假定普通点群共轭类。
可选的 weights 才用于类大小加权；现有 `match_characters` 复用此通用核心。
它返回 raw_multiplicities、通过验证的整数 multiplicities、重构 residual 和 valid。
character 匹配本身不证明闭合或群关系；实际轨道分析优先用 `analyze_little_group`。

结果中的 `irrep_0` 等仅是**当前结果的局部 ID**，不是标准名称；
不同算法/版本的枚举顺序可能改变。

## 6. 标准标签单独映射

三条入口：

- `c3v_labels(reference)`：Γ 点按实际旋转识别恒等、C3、镜面，再用已有表映射
  A1/A2/E，不按 spgrep 枚举序号猜测；任意坐标轴取向都适用。
- `labels_from_character_table(...)`：外部命名表必须提供同一晶胞/原点下的
  rotations、translations、kpoint、逐操作 characters 和 source。
  自动重排操作，并修正整数平移代表元差：
  \(\chi_{\rm ours}=e^{-2\pi i k\cdot(t_{\rm ours}-t_{\rm table})}\chi_{\rm table}\)。
- `irrep_labels(reference, lattice, positions, numbers, kpoint_name)`：
  可选 IrRep 适配器，处理参考晶胞/原点转换并取得 BCS 标准标签。

~~~bash
python -m pip install -e '.[labels]'
python examples/bloch_little_group.py --irrep-labels
~~~

IrRep 适配器隔离使用其内部 `SpaceGroupIrreps` API，明确固定
`irrep==3.3.0`、`irreptables==3.1.0`，升级时必须重新验证。
仅支持表中收录的 maximal k 点，调用者显式给出名称，例如 GM；
缺表、坐标不符、群不匹配或标签匹配不唯一会报错。一般 k 点仍可用 spgrep 分解，
但不承诺所有 k 点都有标准名称。

`IrrepLabels` 保留 ID→名称、来源及匹配残差，reference 保留用于分析的 k、
操作、因子系统和后端版本。命名与物理分解分离，避免外部标签约定渗入数值核心。
IrRep 等第三方包由用户作为可选依赖安装；本仓库不复制其源码或数据表。

## 7. 最小调用顺序

下面的 C_k/S_k/几何/AO 元数据必须来自同次计算。这是接口示意；
完整可运行的合成例子见 `examples/bloch_little_group.py`。

~~~python
rotations, translations = symmetry_from_structure(lattice, positions, numbers)
reference = spgrep_irreps(rotations, translations, kpoint)
operations = [
    build_bloch_ao_operation(
        positions, lattice, species, ao_labels, w, kpoint,
        translation_fractional=t,
    )
    for w, t in zip(reference.group.rotations, reference.group.translations)
]
result = analyze_little_group(C_k, operations, reference, S_k)

# 先检查 result.match.valid；标准标签是单独的一层。
names = irrep_labels(reference, lattice, positions, numbers, "GM")
named_counts = {
    names.labels[identifier]: count
    for identifier, count in result.match.multiplicities.items()
}
~~~

测试覆盖 Γ 极限、纯平移、源原子列相位、非正交晶胞的倒空间变换、
非交换跨 k 乘法、独立周期高斯 AO 重叠的 metric 协变性、边界小群、
滑移面的 ±i、随机复酉混合、可约表示和错误路径。
IrRep 的真实集成测试包含平移原点的 P-1 和换轴的非对称型 Pc 边界点。

## 8. 后续工作

此阶段完成单电子标量 Bloch 对称性。下一步应把同次 ABACUS C(k)、S(R)、
STRU 和 Orbital 元数据接成可复现真实数据测试。
当前仍不支持自动换胞/能带折叠、spinor/double group、反幺正/磁群、ell>2 AO 旋转、
多 k LR 的 ph 映射及完整 X/Y 根表示。多电子态标签还需要参考态表示，
自旋上标仍需要独立的自旋证据。

外部接口依据：
[spgrep API](https://spglib.github.io/spgrep/api/api_core.html)、
[IrRep 源码](https://github.com/irreducible-representations/irrep/blob/master/irrep/spacegroup_irreps.py)。
实际验证版本为 spgrep 0.7.0、spglib 2.7.0、IrRep 3.3.0、irreptables 3.1.0。
