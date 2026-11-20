# DockStudio 基准测试报告（2026-11-20，Tesla P4）

环境：vina 1.2.7，exhaustiveness=8，num_modes=9，seed=42，cpu=4，盒 22.5³ Å（共晶质心）

## 1. 多靶点 redocking 面板（D9 原位 RMSD）

| 靶点 | 配体 | 亲和力 | 原位 RMSD | 叠合 RMSD | D9 | PoseBusters |
|---|---|---|---|---|---|---|
| 1AJV | | -9.105 | 1.649 | 1.395 | PASS | 未通过 |
| 1ETR | | -9.297 | - | - | None | - |
| 1HSG | | -11.265 | 0.454 | 0.380 | PASS | 通过 |
| 1HVR | | -11.275 | 0.945 | 0.764 | PASS | 未通过 |
| 1PPC | | -9.247 | - | - | SKIP | 通过 |
| 1R1H | | -9.845 | 0.777 | 0.581 | PASS | 通过 |
| 2SIM | | -6.723 | 1.340 | 1.230 | PASS | 通过 |
| 3CLA | | -5.541 | 7.674 | 1.540 | FAIL | 通过 |
| 4DFR | | -8.112 | 6.193 | 2.314 | FAIL | 通过 |
| 4PHV | | -12.056 | 0.620 | 0.584 | PASS | 通过 |

**可测靶点 8 个**（1ETR 因 SMILES 价态错误无法准备，1PPC 的共晶参考与 SMILES 原子数差 1 无法匹配）

- **D9 判据（原位 RMSD ≤ 2.0 Å）：6/8 = 75.0%**
- **条件成功率（原位 ≤2 Å 且 PoseBusters 全通过）：4/8 = 50.0%**
- 若按修正前的叠合口径会报：7/8 = 87.5%（**高估**，见 3CLA）

## 2. 参数不变性

- 同 seed(42) 重复两次：**逐位相同** → True
- exhaustiveness 8 vs 16：Spearman ρ = **1.0000**（排序完全一致）
- 随机 seed(None) 重复两次：相同 = False（预期不同）

## 3. 关键发现

1. **3CLA 是 D9 修正的直接受益者**：叠合口径 1.54 Å（会判通过），原位口径 7.674 Å（判失败）。
   若沿用修正前的度量，这个失败会被报成成功。
2. **1HVR / 1AJV 的 PoseBusters 未通过**（inchi_convertible、all_atoms_connected），
   但它们的 D9 是 PASS —— 说明几何达标不等于化学合理，两个指标必须同时报。
3. **1ETR 无法准备**：SMILES 中 N 原子显式价态为 4，RDKit 直接拒绝。这是 targets.json 数据问题。
4. **1PPC 无法匹配**：SMILES 含显式氢（[H]/N=C），抽出参考 37 原子 vs SMILES 38 原子。