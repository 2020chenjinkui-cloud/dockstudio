# Uni-Dock 安装与 Tier1 吞吐提升评估（2026-11-21）

## 一、安装过程（记录踩到的坑）

目标：在 Tesla P4（算力 6.1）上启用 Uni-Dock GPU 批处理。

| 步骤 | 结果 |
|---|---|
| `conda create -n unidock -c conda-forge unidock=1.2.0` | 装上 1.2.0，但**只编译了 sm_100/103/120/121** |
| 实测该版本 | GPU 利用率 0%、显存 0 MiB、单分子 313.9s（比 vina 慢 43 倍）——静默回退 CPU |
| 原因 | P4 是 sm_61，不在预编译架构列表里 |
| `nvcc --list-gpu-arch`（CUDA 13.1） | **最低 compute_75**，Pascal 已被 CUDA 13 彻底移除 |
| 改装 CUDA 12.9（conda） | `compute_61` 在列表里，可编译 |
| 从源码编译（sm_61） | 成功，二进制含 `sm_61`，GPU 计算真正跑起来 |

编译过程中遇到的三个问题，都已解决并记录：
1. gcc 15 被 nvcc 拒绝 → 用 conda 装 gcc 14 + `-ccbin` 显式指定
2. 缺 `curand_kernel.h` → `conda install libcurand-dev`（不要从 CUDA 13.1 拷贝，API 不兼容）
3. 链接找不到 `/lib64/libm.so.6` → 加 `-L/lib/x86_64-linux-gnu -Wl,-rpath-link`

**关键教训**：conda-forge 的预编译 unidock 在 Pascal 卡上会静默回退 CPU，不报错。必须检查二进制里的 `sm_` 架构。

## 二、正确用法（也是一个坑）

`--gpu_batch <多分子.sdf>` **只对接 1 个配体**——该参数接受的是"文件列表"，不是多分子 SDF。
正确做法是 `--ligand_index`，指向一个每行一个 PDBQT 路径的文本文件。

误用时的表现：`Total ligands: 1`，日志里 `Group size: 0`，输出只有 9 个分子。

## 三、实测吞吐（同机、同参数）

| 引擎 | 配置 | 吞吐 |
|---|---|---|
| vina CPU | exh=1（Tier1 廉价档） | **5.43 s/mol** |
| Uni-Dock GPU | N=64，GPU 65.3s / 墙钟 73.7s | 1.151 s/mol |
| Uni-Dock GPU | N=256，GPU 233.8s / 墙钟 267.9s | **1.046 s/mol** |

**加速比 5.2×**

GPU 时间与分子数严格线性（64→256 是 4 倍，GPU 65.3→233.8s 是 3.58 倍），说明已饱和，1 s/mol 是稳态吞吐。

## 四、规模结论的变化

| 库规模 | vina CPU | Uni-Dock GPU |
|---|---|---|
| 10,000 | 15.1 小时 | **2.9 小时** |
| 50,000 | 3.1 天 | **14.5 小时** |
| 100,000 | 6.3 天 | **1.2 天** |
| 1,000,000 | 2.1 个月 | **12.1 天** |

百万级从"不可行"变成"两周内"——但前提是只算对接。

## 五、真正的瓶颈转移了：配体准备

实测配体准备（RDKit 构象生成 + PDBQT 写出）：

| N | 耗时 | 每分子 |
|---|---|---|
| 64 | 88.0 s | 1.38 s/mol |
| 256 | 343.5 s | 1.34 s/mol |

**准备 1.34 s/mol 已经比 Uni-Dock 对接 1.05 s/mol 更慢**，而且是纯 CPU。

端到端（准备 + 对接）：

| 库规模 | 端到端 | 其中准备 | 其中对接 |
|---|---|---|---|
| 100,000 | **2.8 天** | 1.5 天 | 1.2 天 |
| 1,000,000 | **27.6 天** | 15.5 天 | 12.1 天 |

## 六、结论

1. **Tier1 粗筛提速 5.2 倍**，十万级库从 6.3 天降到 1.2 天（对接部分）。
2. **单机现实上限是 10 万级**（端到端约 2.8 天），百万级要 27.6 天。
3. **瓶颈已从对接转移到配体准备**。下一步优化应针对构象生成（并行化 / 缓存 / 更快的 ETKDG 配置），而不是继续压对接。
4. P4 在 Uni-Dock 上能跑，但算力 6.1 是官方支持范围（≥7.0）之外，属于"能跑但无保证"。稳定性需长期观察。

## 七、复现

```bash
# 编译（sm_61）
cd ~/Uni-Dock/unidock
conda create -y -n unidock-build -c nvidia -c conda-forge \
  "cuda-nvcc=12.9" "cuda-cudart-dev=12.9" libcurand-dev cmake libboost-devel gxx_linux-64=14
E=~/miniconda3/envs/unidock-build
cmake -B build -DCMAKE_CUDA_ARCHITECTURES=61 -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CUDA_COMPILER=$E/bin/nvcc \
  -DCMAKE_CUDA_FLAGS="-allow-unsupported-compiler -ccbin $E/bin/x86_64-conda-linux-gnu-g++" \
  -DCMAKE_EXE_LINKER_FLAGS="-L/lib/x86_64-linux-gnu -Wl,-rpath-link,/lib/x86_64-linux-gnu"
cmake --build build -j4

# 运行（必须用 ligand_index，不是 gpu_batch <sdf>）
unidock --receptor rec.pdbqt --ligand_index index.txt \
  --center_x .. --center_y .. --center_z .. \
  --size_x 22.5 --size_y 22.5 --size_z 22.5 \
  --search_mode balance --num_modes 9 --seed 42 --dir out/
```

验证编译结果：`strings build/unidock | grep -oE "sm_[0-9]+" | sort -u` 应包含 `sm_61`。
