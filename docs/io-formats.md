# IO 合约与支持范围

## 统一约定

- 空间轨道 C 的形状 `(nao, nband)`；复数可以传入物理核心。
- 两个 spin block 必须共享 AO 基底与 AO 顺序；alpha/beta 不由文件名推断。
- 文件 band 标签保存在 `band_numbers`（从 1 开始），NumPy 列号从 0 开始。
- Γ 点波函数能量保留 Ry（`energies_ry`）；交换 LR 能量明确为 eV。
- AO 重叠、旋转和 C 必须来自相同结构、离子步、AO 基和 k 点。
- 不根据占据字段自动生成单行列式：ABACUS 字段可能涉及自旋简并或 k 权重，
  应结合输出设置解释；分数占据不能直接用于当前 `determinant_s2`。

## 原生 ABACUS：已支持一个窄格式

`AbacusReader.read_wavefunctions(path, format="gamma_text", spin="beta")`
读取 `wfc_nao_write2file` 的实数单帧文本：

```text
2 (number of bands)
3 (number of orbitals)
1 (band)
-1.0 (Ry)
1.0 (Occupations)
1.0 0.0 0.0
2 (band)
0.2 (Ry)
0.0 (Occupations)
0.0 1.0 0.0
```

系数允许换行，支持 e/E/d/D 指数；逐 band 验证编号、系数数量和有限数值。
拒绝复数多 k 头、二进制、截断文件、重复/追加帧。
当前回归测试是依 writer 布局人工构造的 fixture，尚未进行真实 ABACUS 端到端验证。

## stateid 数组交换：不是 ABACUS 原生格式

### S：NPY

通过 `np.save("S.npy", S)` 导出完整稠密矩阵，以
`reader.read_overlap(path, format="stateid_npy")` 读取。
验证方阵、Hermitian 与正定性；不读取对象数组，不从 R 空间自动作 Fourier 变换。
AO 来源与顺序由调用者负责，初版文件本身没有附带 AO 元数据。

### LR：NPZ schema_version=1

```python
import numpy as np
np.savez(
    "lr.npz", schema_version=1, method="TDA",
    X=np.eye(2, dtype=complex),        # (nph,nroot)，列是根
    energies_ev=np.array([2.1, 2.1]),
    transitions=np.array([[0, 1, 125, 1, 126], [0, 1, 125, 1, 127]]),
)
```

`transitions` 每行为 `(k, spin_hole, occupied_band, spin_particle, virtual_band)`。
所有索引均从 0 开始；spin 0=alpha、1=beta；band 指全局 KS 列号，不是活性空间局部号。
允许显式保存 spin-flip 行，但本版没有 spin-flip 物理分析器。
同一行不得重复；同 spin 的 occupied/virtual 不能指向同一轨道。
TDA 必须没有 Y；`method="full_lr"` 必须有同形状 Y。
读取时不会自动归一化，保存 full_lr 也不表示已经支持其自旋计算。

这个格式要求导出者完成全局 ph 映射，明确 X/Y 的物理含义。它不是
`Excitation_Amplitude_*_<rank>.dat` 的简单改名。无 pickle；未知键和错误版本明确报错。

## 明确的占位入口

`read_output`、原生 `read_overlap(..., format="abacus_csr")`、
`read_lr_eigenvectors(..., format="abacus_lr")` 均抛出 `UnsupportedFormatError`。
`OutputReader` 是未来其他软件后端的协议。

## 下一项任务：原生 S(R) → S(Γ)

优先实现只读新旧 CSR 检测、选择离子步、解析每个 R 的 values/columns/row pointers；
验证行指针、列范围、矩阵维度、重复 R、复数格式和截断数据。
Γ 点使用 \(S_\Gamma=\sum_R S(R)\)，不能把单个 R=0 块当作总重叠。
一般 k 点未来需要与 ABACUS Bloch 约定一致的
\(S(k)=\sum_R e^{2\pi i k\cdot R}S(R)\) 及原子内相位约定核验。
不要假设每个 S(R) 块自身 Hermitian；应检查 S(-R)=S(R)† 和求和后的 S(Γ)。

验收应使用同一小体系/同一离子步的真实 C、S(R)、STRU、轨道描述文件：
检查 AO 维度、S(Γ) Hermitian/正定、C†S(Γ)C≈I，并对照 ABACUS 或独立解析器。
这会先打通真实输入的度量，随后才能可信地组装与检验 AO 对称操作。
