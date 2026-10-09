# Boltz-2 AI 亲和性复核

DockStudio **不内置** Boltz-2（torch + 模型权重合计数十 GB），而是由用户自备安装，
后端通过发现机制调用。AI 复核的语义：把「蛋白序列 + 配体 SMILES」交给 Boltz-2 独立预测
亲和性，作为与对接打分交叉验证的信号。**复核不读取姿态三维坐标**——不同姿态的复核输入相同，
因此通常只需对 top 姿态复核一次（提交时 `pose_indices` 留空=全部姿态，建议手动只选 1 个）。

## 安装方式一：独立 venv（推荐）

```bash
python -m venv boltz-venv
boltz-venv/Scripts/activate        # Windows; Linux/macOS: source boltz-venv/bin/activate
pip install boltz
```

然后告诉 DockStudio 这个 Python 在哪：

```bash
# Windows (PowerShell)
setx DOCKSTUDIO_BOLTZ "C:\path\to\boltz-venv\Scripts\python.exe"
# Linux/macOS
export DOCKSTUDIO_BOLTZ=/path/to/boltz-venv/bin/python
```

DockStudio 每次提交复核任务前会用 `python -c "import boltz"` 验证该环境。

## 安装方式二：conda

```bash
conda create -n boltz python=3.11 -y
conda activate boltz
pip install boltz
# 记录 python 位置：which python / where python，配到 DOCKSTUDIO_BOLTZ
```

## 安装方式三：命令前缀（WSL / Docker / 自定义）

`DOCKSTUDIO_BOLTZ_CMD` 优先级最高，DockStudio 把它当作命令前缀直接拼
`predict input.yaml --out_dir ...`：

```bash
# WSL 里已装好 boltz
export DOCKSTUDIO_BOLTZ_CMD="wsl boltz"
# Docker 镜像
export DOCKSTUDIO_BOLTZ_CMD="docker run --rm --gpus all -v /data:/data ghcr.io/jwohlwend/boltz boltz"
```

（docker 方式需要工作目录已挂载进容器，out_dir 会默认挂在 job 的 ai/<task_id>/ 下。）

## GPU 要求与 CPU 耗时

- Boltz-2 强烈建议 NVIDIA GPU（8GB+ 显存）；CPU 可跑但**单条复合物通常 30 分钟到数小时**，
  且首次需下载数 GB 权重。
- `DOCKSTUDIO_BOLTZ_CMD` 模式可自行追加 `--accelerator gpu` 等参数（拼到前缀末尾）。

## Hugging Face 下载失败

权重/依赖从 huggingface.co 拉取。国内网络失败时 Boltz 会因连接错误报错，
DockStudio 会在错误提示里追加建议：

```bash
export HF_ENDPOINT=https://hf-mirror.com
```

## 隐私说明（MSA 服务器）

复核默认走 **ColabFold 公共 MSA 服务器**（`use_msa_server=true`）：
**蛋白残基序列会上传到公共服务器**。不能用公共服务的场景请：
1. 在提交复核时把 `use_msa_server` 设为 `false`（本地 MSA 慢得多，需自备 mmseqs 数据库）；
2. 或完全在隔离环境用 `DOCKSTUDIO_BOLTZ_CMD` 自定义命令并在前缀中传 boltz 的离线参数。
当 `msa_server_used=true` 时，复核任务结果会带该标记，前端应提示用户。

## 输出字段

Boltz-2 affinity 预测输出 `affinity_*.json`（DockStudio 解析 `affinity_pred_value` 与
`affinity_probability_binary`）：
- `log_ic50`（log10 IC50，μM 基准；越低预测结合越强）
- `binder_probability`（0-1，是否 binder 的概率）

参考：Boltz 官方仓库 README（boltz predict 输入 YAML 与 affinity properties 小节）。
