# 01：从系数到物理对称性

## 1. 正确的内积

非正交 AO 中不能用 c†c 判断归一化；必须用 c†Sc。
对 d 维子空间 C，Gram 矩阵 G=C†SC 的谱分解 G=VgV† 给出
Q=CVg⁻¹ᐟ²V†。Q†SQ=I，而 Q 与 C 张成同一空间。
`orthonormalize` 完成这个过程，秩不足会报错，防止把二维问题静默改为一维。

## 2. 投影而不是盲目给矩阵贴标签

定义 T 为 AO 系数变换，M=ST 为算符矩阵元，则 D=Q†MQ。
Q D 是 TQ 在选定子空间的投影。若 TQ≠QD，选定轨道不构成不变子空间；
此时 D 只是压缩矩阵，不是该群在此空间的精确表示。

`project_representation` 返回 D、Tr D、闭合误差和 AO 度量保持误差。
度量保持用 A†A=S 的 A 构造 \(\widetilde T=ATA^{-1}\)，检查
\(\|\widetilde T^\dagger\widetilde T-I\|_F/\sqrt{N_{AO}}\)，避免 AO 缩放影响诊断。
`analyze_c3v` 检查全部六个操作、类内一致性及群关系，再做匹配。

输入为数值畸变结构时，应先辨别几何对称性破缺与纯数值误差。
默认 character/闭合容差为 1e-5，可显式调整，但调整必须有精度或结构依据。
默认 AO/Gram 秩阈值 1e-10；稠密线性代数目前适合原型验证，尚未优化大基组性能。

## 3. 用二维 E 作为最小例子

令 \(\theta=2\pi/3\)，在实基中

\[
D(C_3)=\begin{pmatrix}\cos\theta&-\sin\theta\\\sin\theta&\cos\theta\end{pmatrix},
\quad D(\sigma_v)=\begin{pmatrix}1&0\\0&-1\end{pmatrix}.
\]

所以 χ(E)=2、χ(C3)=2cosθ=−1、χ(σv)=0。
换成复线性组合，C3 可以对角化为 exp(±2πi/3)，但镜面仍会交换两分量；
一个 C3 本征轨道仍不足以构成 C3v 的一维 irrep。

C3v 操作名在代码中固定为 `E, C3, C3^2, sv1, sv2, sv3`，
其中 sv2=C3·sv1、sv3=C3²·sv1（右端先作用）。检查
r³=I、f²=I、frf=r²，避免把同类标签或操作方向弄错。

运行 `python examples/c3v_minimal.py`。`tests/test_symmetry.py` 进一步测试复数酉旋转、
非正交 AO 基变换、不闭合的单根 E、错误的群关系和 noisy/unmatched characters。

## 4. 重数与证据强度

\[
n_\alpha=\frac1{|G|}\sum_c |c|\chi_\alpha(c)^*\chi(c).
\]

先计算非整数的数值重数，再检查是否接近非负整数及重构误差。
`match_characters([3,0,1])` 返回 A1+E，不把 3 维空间硬说成某一根态。
闭合小、群关系成立、类内一致、重数近整数是互补检验；
“能量很接近”是选择候选子空间的线索，不是 irrep 的定义。

## 5. 从轨道拓展到 TDA 根

对一个对称性纯的一维参考态 \(R|\Phi_0\rangle=\chi_0(R)|\Phi_0\rangle\)，
以 \(|\Phi_i^a\rangle=a_a^\dagger a_i|\Phi_0\rangle\) 为基底，

\[
T_{ph}(R)=\chi_0(R)\,[D_{occ}(R)^*\otimes D_{virt}(R)].
\]

`particle_hole_operation` 实现这一完整张量积、单 spin block 情形，虚轨道索引变化最快。
把一组 TDA 根的 X 当作 ph 空间中的列系数，即可用同样的投影方法研究根子空间。
若只关心跃迁算符可省去 χ0，但识别总多电子态时必须计入参考态。
ABACUS 自旋适配通道、占据数差异、截断激发集合和行列式相位都需要适配器明确处理。
完整 LR 的 X/Y 尚不能走这个普通 Hilbert 空间流程。
