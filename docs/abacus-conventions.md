# 本地 ABACUS 约定核查

核查日期：2026-10-06。参考仓库 `/home/fortneu/abacus/abacus-develop`，HEAD：
`119f1b4bbc064a57a288fd52b3053f11436f465f`。
该工作区已有用户未提交改动；本项目仅作只读参考，没有修改、暂存或提交 ABACUS 文件。
下列重点参考文件在核查时没有工作区改动。

## 旋转与实球谐

`source/source_lcao/module_ri/module_exx_symmetry/symm_rotation.cpp` 中：

- `m2im`（约 325 行）把 m 映射到 `(0,+1,-1,+2,-2,...)`。
- `ovlp_Ylm_Slm`（约 257 行）定义 \(B_{mm'}=\langle Y_l^m|S_l^{m'}\rangle\)。
- `cal_rotmat_Slm`（约 331 行）构造 \(T_l=B^\dagger D_lB\)。
- `get_euler_angle`（约 280 行）的注释明确说明原子位置采用行向量，涉及转置约定。
- 反常旋转使用 \(D_l(R)=(-1)^lD_l(-R)\)，即先用正旋转 \(-R\) 再乘宇称。
- `contruct_2d_rot_mat_ao`（约 449 行）结合原子映射、壳层块和返回晶格矢量的
  \(e^{-2\pi i\mathbf k\cdot\mathbf O}\) 相位。

本库 `complex_to_real_basis` 由以下数学定义独立实现（m>0）：

\[
S_l^m=\frac{Y_l^m+(-1)^mY_l^{-m}}{\sqrt2},\qquad
S_l^{-m}=\frac{-iY_l^m+i(-1)^mY_l^{-m}}{\sqrt2},\quad S_l^0=Y_l^0.
\]

`rotation_in_real_basis` 只作 B†DB，不解释输入旋转的主动/被动意义。
测试中的 l=1 四分之一转动矩阵与本地
`.../test/symm_rotation_test.cpp` 的 `c_dagger_D_c_C41_ref` 对照。
这验证了基变换约定，并不代表已验证整个晶体 AO 旋转构造。

## 最容易混淆的两个 M

ABACUS `Symmetry_rotation` 的 AO 旋转对象也称 M，用于密度矩阵恢复；
源码涉及列主序并行存储、行向量旋转，以及 `rot_matrix_ao` 的具体变换方向。
不能仅因变量名相同，就将它当作本库的
\(M_{\mu\nu}=\langle\phi_\mu|\hat R|\phi_\nu\rangle\)。

本库固定列系数、主动变换 \((\hat R\psi)(r)=\psi(R^{-1}r)\)，输入 T 满足
\(C'=TC\)，矩阵元 M=ST。未来连接 ABACUS 时，必须通过已知 s/p/d 轨道及
非交换旋转检验转置、逆、共轭和原子映射，而不是盲目复制一个矩阵。
还需验证 \(T^\dagger ST=S\) 和 \(T(R_1R_2)=T(R_1)T(R_2)\)。
单靠 C3v 的实 character 不能区分 R 与 R⁻¹ 的全部方向约定。

## IO 参考文件

- `source/source_io/module_wf/write_wfc_nao.cpp`：`wfc_nao_write2file` 的单帧
  Γ 点实数文本，nbands/nlocal 后逐 band 写出 Ry、原始 occupations、AO 系数。
  每轨道输出为一组数，本库读取后转为 `(nao, nband)` 列存储。
- `source/source_esolver/esolver_lr_lcao_tddft.cpp`（约 508 行起）：
  `Excitation_Energy_<label>.dat` 和带 rank 后缀的
  `Excitation_Amplitude_<label>_<rank>.dat`，包括 `openshell` 分支。
  后缀不能被忽略：必须确定 `nloc_per_state`、全局 ph 次序和并行分布。
- `source/source_base/module_out/csr_reader.cpp`、
  `source/source_io/module_hs/write_hs_r.cpp` 和
  `docs/advanced/interface/migration-guide-csr-format.md`：新旧 CSR 头、离子步、
  晶胞和 S(R) 数据块；这是下一项 IO 工作的直接参考。

ABACUS 版本间文件名和布局有变化。本库依据内容和明确的 `format` 选择，
不会仅依据 `WFC_NAO_GAMMA...` / `wfs..._nao.txt` 文件名宣称兼容。
公开参考：[ABACUS 波函数输出文档](https://abacus.deepmodeling.com/en/v3.7.3/advanced/elec_properties/wfc.html)、
[较新版本初始化文档](https://abacus.deepmodeling.com/en/latest/advanced/scf/initialization.html)。
实际接口开发以记录的本地 writer 和真实小型 fixture 为准。
