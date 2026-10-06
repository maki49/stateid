# 02：自旋诊断与 LR 的边界

取 ℏ=1。由
\(\hat S^2=\hat S_z^2+\tfrac12(\hat S_+\hat S_-+\hat S_-\hat S_+)\)，
整数占据共线行列式具有 \(M_S=(N_\alpha-N_\beta)/2\)，交换项由两个自旋通道
空间轨道的重叠决定：

\[
\langle S^2\rangle=M_S^2+\frac{N_\alpha+N_\beta}{2}
-\sum_{i\in occ\alpha,j\in occ\beta}|\langle\psi_i^\alpha|\psi_j^\beta\rangle|^2.
\]

在同一个非正交 AO 基底中，交叉重叠是 Cα†SCβ。
当 Nα≥Nβ 时，也可写为 \(S_z(S_z+1)+N_\beta-\sum|O_{ij}|^2\)。
本库采用对 alpha/beta 对称的前一种形式，避免 beta 多数时误用公式。
公式与 [PySCF 的 UHF spin_square 推导](https://pyscf.org/_modules/pyscf/scf/uhf.html#spin_square)一致。

```python
import numpy as np
from stateid.spin import determinant_s2

c = np.eye(2)
assert abs(determinant_s2(c, c)) < 1e-12             # 闭壳层：0
assert abs(determinant_s2(c, np.empty((2, 0))) - 2) < 1e-12  # 两个同向自旋：2
assert abs(determinant_s2(c[:, :1], c[:, 1:]) - 1) < 1e-12   # 破缺对称行列式：1
```

同一 spin block 的占据轨道必须正交归一。Cα†SCβ 则通常不为单位阵；这正是
自旋污染信息，不能把两个自旋块一起正交化而抹掉它。
`determinant_s2` 不接收部分占据权重、不做整数占据猜测、不适用于 SOC spinor。
对 KS 轨道的这个诊断作用于辅助行列式，不能替代真实多电子自旋信息。

## TDA 还缺什么

对明确的正交行列式/CSF 展开 \(|\Psi\rangle=\sum_I X_I|\Phi_I\rangle\)，

\[
\langle S^2\rangle=\frac{X^\dagger S^2_{ph}X}{X^\dagger X},\qquad
(S^2_{ph})_{IJ}=\langle\Phi_I|\hat S^2|\Phi_J\rangle.
\]

已有 `tda_s2` 只收缩用户提供的 Hermitian **绝对** S² 算符，不能传入仅激发增量
并误称为总 S²。未来 `build_s2_ph` 需要全部已占据 alpha/beta 轨道、相关虚轨道、
全局 ph 索引、参考行列式与费米符号约定，必要时包含自旋适配变换。
这里的 CI 式解释须明确是采用的 TDA 状态近似，不能未经证明推广到所有 TDDFT 输出。

完整 LR 的常见归一化是 X†X−Y†Y=1；这个不定度量不等于普通概率权重。
目前 `lr_s2` 拒绝计算，后续应选定并文档化激发态响应密度/自旋期望值的定义与推导，
而不是把 X 和 Y 拼起来交给 `tda_s2`。

## 物理归属的用法

理想 S=0、1/2、1 对应 0、3/4、2。偏离预期可能提示自旋污染、输入占据错误或近似限制。
仅仅“接近 2”不是普遍的三重态纯度证明；需要参考态、自旋通道和计算方法共同约束，
有条件时进一步检查自旋投影或 S² 方差。上标 ³ 不能由单个 a1→e 轨道跃迁直接推断。
