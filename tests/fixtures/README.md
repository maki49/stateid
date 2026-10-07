`gamma_real.txt` 是人工构造的三轨道测试数据，按本地 ABACUS
`wfc_nao_write2file` 实数文本布局编写，不是 ABACUS 计算结果，也不是 NV⁻ 数据。
参考版本和源文件见 `docs/abacus-conventions.md`。

`nv_avg_tda.dat` / `nv_jt_tda.dat`：真实 63 原子 PBE NV⁻ 终点输出，
来源 lr-grad/run_data/18_NV_diamond/222/pbe/{02_3E_avg,03_3E_jt}/OUT.ABACUS/
trans_analysis_updown_tda.dat。仅含能量、偶极和阈值截断跃迁表，不含波函数。
