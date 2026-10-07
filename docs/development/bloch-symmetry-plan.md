# Bloch symmetry implementation plan — 2026-10-07

用户提出四项扩展：B(g,k)、小群、spgrep/通用分解、标准标签。
本轮将它们实现为单电子标量流程，保留原有 C3v/TDA API。

## 已落地

| 阶段 | 实现 | 验收证据 |
|---|---|---|
| 坐标和操作契约 | 分数坐标 W/t/k，k'=W^-T k；显式原胞与 cell gauge | 倾斜晶胞、边界倒格矢、错误/缺失操作 |
| Bloch B | 复用 Γ AO builder；每源 AO 列的返回晶格相位；可缓存 T | Γ 极限、纯平移、反演相位、滑移平方、非交换跨 k 复合 |
| 小群 | 操作筛选、乘法表、晶格余项与 factor system；可选 spglib 找群 | 一般 k/高对称轴/边界点；不闭合与非原胞拒绝 |
| spgrep 后端 | scalar complex small irreps；验证顺序、群关系、正交性和完备性 | C3v 与 nonsymmorphic glide，随机复酉规范 |
| 一般分解 | per-operation characters；复用 class-weighted 旧 API | 整数重数、重构误差、可约表示与负例 |
| 标准标签 | C3v builtin、外部带操作表、固定版本 IrRep adapter | 实际 IrRep P-1 原点移动、Pc 换轴边界标签 |
| 跨 k 投影 | 独立 source/target metric，复用旧投影数值核心 | 独立 gauge；周期高斯 AO overlap 协变性 |
| 文档示例 | 教程07、双语 README、bloch_little_group.py | Γ A1/E/A1+E，边界 glide 与可选标准名 |

原有 70 tests；本轮扩展后 94 tests 全通过，包括真实可选后端调用。
基础依赖与可选依赖分离；不安装后端时基础数学可用，集成测试明确 skip。
禁用四个可选包导入的隔离检查：83 tests 通过，11 个后端集成 tests 明确 skip。
可编辑安装、新旧两个示例、git diff --check 均通过。
本轮未修改 ABACUS 或 pyATB；未复制第三方群表/源码。

## 设计选择

- spgrep 负责生成数学表示，stateid 负责实际 AO 波函数投影和物理诊断。
- 标签和数值分解分离。irrep_0 等是本次结果 ID，绝非固定 A1/Γ1 编号。
- IrRep 的标准命名复用内部 API，因此约束 irrep 3.3.0 / irreptables 3.1.0。
  查不到的 k 点不猜名称；仍可保留 spgrep 数值分解。
- 非对称型空间群的有限代表元满足带相位的乘法，而非普通点群乘法。
- 同 k、跨 k 使用共享投影数值核心；后者没有规范不变的单独 character。
- 不自动把结构转换成原胞，以免与已有 C/S/AO 排列失配。
- AO 壳层仍限 ell<=2；磁群、spinor 和多 k LR 是独立后续工作。

## 下一步优先级

1. 同一真实 ABACUS 计算的 STRU / Orbital / C(k) / S(R) 联合适配与回归。
   必须固定离子步、k 坐标单位、原子坐标代表、AO 行顺序及基函数约定。
2. 用真实非 Γ 波函数检查 B†S(k')B=S(k) 和波函数闭合；
   当前新增非 Γ 数值测试为解析/合成模型，不能替代 ABACUS 实测。
3. 多 k TDA：保留全局 ph 映射，旋转并排列 k，组合 D_occ* 与 D_virt，
   最后投影完整根子空间并结合参考态表示。避免构造 ph² 矩阵。
4. 扩展可用标准标签覆盖与跨版本适配；需要时再处理换胞/折叠、double groups、
   磁群与反幺正；这些不是给当前 scalar API 加一个布尔开关就能安全完成的。

## 验证环境

.venv/bin/python（Python 3.13），NumPy 2.4.6、SciPy 1.17.1、
spgrep 0.7.0、spglib 2.7.0、IrRep 3.3.0、irreptables 3.1.0。
可选包安装到仓库已有 .venv，系统环境没有被卸载或替换。
IrRep 集成调用会产生上游 irreptables 文件句柄 ResourceWarning 和 spglib
DeprecationWarning；不影响当前断言结果，未修改外部包掩盖提示。

运行：`.venv/bin/python -m unittest discover -s tests -v`。
示例：`.venv/bin/python examples/bloch_little_group.py --irrep-labels`。
