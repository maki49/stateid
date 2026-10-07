# ABACUS TDA 跃迁组成

读取 `trans_analysis_*_tda.dat` 的实数幅度和原始 `Excitation rate`：

```python
from stateid.io.abacus_transitions import read_transition_analysis
states = read_transition_analysis('OUT.ABACUS/trans_analysis_updown_tda.dat')
print(states[0].summary())
```

命令行支持多个文件及根选择：

```bash
python -m stateid.io.abacus_transitions file1.dat file2.dat --states 0 1 --axis 1 1 1
```

根编号保持 ABACUS 的零基编号；轨道、k 点保持一基编号，a/b 分别对应 alpha/beta。
续行继承上方根编号。输出 `reported_weight` 只累加文件中打印的权重，不能重新
归一化；ABACUS 按 `abs(X) > ana_thr` 筛选，未打印的跃迁不等于零。调用前须从
INPUT 确认 TDA，不能把完整 LR 仅输出 X 的分析当作 CI 权重。

这项功能无需 MPI 振幅重组。它不读取原始振幅分片，不自动赋予轨道 irrep 或
激发态自旋，也不证明跨结构态跟踪。简并根的逐跃迁分配依赖根和轨道的基选择，
应同时报告完整根对到候选轨道子空间的总权重；打印截断意味着该权重仍是近似值。

真实 NV⁻ 63 原子 PBE 数据（lr-grad/run_data/18_NV_diamond/222/pbe）：

| 结构 | 根 | E/eV | ↓126→↓127 | ↓126→↓128 | 两项合计 |
|---|---:|---:|---:|---:|---:|
| 02_3E_avg 终点 | 0 | 1.79077 | 0.828112 | 0.135772 | 0.963884 |
| 02_3E_avg 终点 | 1 | 1.79086 | 0.135977 | 0.828194 | 0.964171 |
| 03_3E_jt 终点 | 0 | 1.72661 | 0.966687 | 未打印 | ≥0.966687 |
| 03_3E_jt 终点 | 1 | 1.86788 | 未打印 | 0.959353 | ≥0.959353 |

数据为输出舍入值；候选 a1↓→e↓ 激发组分得到支持，但 a1/e 标签仍须独立
轨道对称性分析。TDA 参考的 Ms=1 不能单独证明激发态自旋纯度。

`state.longitudinal_fraction([1,1,1])` 计算偶极沿 NV 轴的平方占比；
两个终点前两根都小于 1e-6（打印偶极精度下），与横向跃迁相容。
零偶极态的方向占比无定义，函数会明确拒绝。
