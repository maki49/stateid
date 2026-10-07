# stateid 开发上下文

先读 README.md 和 docs/development/excited-state-irrep-handoff-2026-10-07.md。
项目用 Python ≥3.10、NumPy/SciPy，src 布局，测试使用 unittest。

目标：ABACUS-first 的可审计电子态识别；物理核心与软件 IO 分离。
列系数 C，AO metric S；主动系数操作 T，矩阵元 M=ST。
轨道/根子空间必须验证 metric、闭合与群关系后再赋 irrep。
单 E 分量、近简并或 band 编号不能替代完整子空间诊断。
总多电子 irrep 需要参考态相位，自旋上标需要独立证据。

当前新增 Γ s/p/d AO builder、Orbital 标签 reader、矩阵无关 TDA C3v 分析。
真实 NV C/S/Orbital、STRU 解析、原生 LR MPI metadata/adapter 仍缺失。
详见 docs/tutorials/06-ao-tda-symmetry.md 和交接文档最新补充。
缺数据或不支持的物理/格式须明确失败，不推断布局或虚构结果。
不修改外部 ABACUS 源码；已有用户未提交改动应保留。

验证：`PYTHONPATH=src python3 -m unittest discover -s tests -v`。
新增数值功能覆盖解析解、不变量、错误输入及独立布局对照。
