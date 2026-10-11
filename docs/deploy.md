# DockStudio 部署指南

## 服务器模式（Docker）

最简单的方式：

```bash
cd dockstudio
docker build -f docker/Dockerfile -t dockstudio .
docker run -d -p 8000:8000 -v dockstudio-data:/data dockstudio
```

或 `docker compose up -d`（在 `docker/` 目录下）。容器启动后访问 `http://localhost:8000/api/health` 应返回 `{"status":"ok"}`。

数据保存在 named volume `dockstudio-data`（容器内 `/data`），删容器不丢任务数据。

## 局域网访问

默认容器把 8000 映射到宿主机所有网卡，局域网内其他机器用 `http://<宿主机IP>:8000` 即可访问。
**安全提醒**：M0-M3 的 API 无认证，仅适合可信内网。如需加固，建议在反向代理（nginx）上加
Basic Auth 或 IP 白名单；`/api/file` 有路径沙箱但任务/项目数据任何能访问端口的人都能读删。

## GPU 引擎（gnina / Uni-Dock）

gnina 与 Uni-Dock 仅支持 Linux + NVIDIA GPU。`docker-compose --profile gpu up` 启动 gpu-worker
（需要宿主机装好 NVIDIA 驱动 + NVIDIA Container Toolkit，并把 gnina/unidock Linux 二进制放到
`resources/bin/linux-x86_64/`）。

### 实测兼容性（2026-10，Tesla P4 / 驱动 580.178 / CUDA 12.9 容器环境）

| 引擎 | 官方产物 | 结论 |
|---|---|---|
| gnina v1.3.3 | `gnina.cuda12.8.static` | ✅ 二进制非完全 static，需 `libcudnn.so.9`：从 PyPI 拿 `nvidia_cudnn_cu12` 轮子解压取 `nvidia/cudnn/lib` 加入 `LD_LIBRARY_PATH` 即可运行（其余 CUDA 库 NGC CUDA 12.9 镜像自带）。Tesla P4（CC 6.1）上 1HSG 实测 5 pose、CNN 重打分正常，75 秒完成 |
| Uni-Dock v1.2.0 | `dptechnology/unidock` 镜像 | ❌ P4（CC 6.1）进入 GPU 对接阶段即挂死（sm_61 无内核映像）。官方产物需要 **CC ≥ 7.0（Volta+）**；P4 及以下请勿启用，或自行以 sm_61 从源码编译 |

建议 gnina 的容器启动方式：以 `dptechnology/unidock`（NGC CUDA 12.9）或 `nvidia/cuda:12.8.1-cudnn-runtime-*` 为基础镜像 + 上述 cuDNN 轮子库路径。

### gnina 依赖库完整清单（DockStudio 适配器路径实测，2026-10-11）

`gnina.cuda12.8.static` *并非*完全静态——在 python-slim 类精简环境按依赖逐层补齐的顺序实测如下，
**全部缺一不可**（后 3 个是 ldd 不可见的运行时 dlopen/JIT 依赖，漏了报 `c10::DynamicLibraryError` 或 nvrtc 编译失败）：

```
libcudnn.so.9  libcudart.so.12  libcusparse.so.12  libcufft.so.11
libcublas.so.12  libcublasLt.so.12  libcusolver.so.11  libnvJitLink.so.12
libnvrtc.so.12  libnvrtc-builtins.so.12[.alt]
```

获取方式：前 8 个从任一 NGC CUDA 12.x 容器（如 `dptechnology/unidock`）的
`/usr/local/cuda*/lib64` 拷贝；cudnn 从 PyPI 的 `nvidia_cudnn_cu12` 轮子解压；
然后：

```bash
export LD_LIBRARY_PATH=/path/cudalib:/path/cudnn/lib  # 两个目录拼一起
export DOCKSTUDIO_GNINA=/path/gnina
docker run --gpus all ...  # 或宿主机直跑
```

DockStudio 容器实测命令：

```bash
docker run -d --gpus all -p 8000:8000 \
  -v /path/engines:/eng \
  -e DOCKSTUDIO_GNINA=/eng/gnina \
  -e LD_LIBRARY_PATH=/eng/cudalib:/eng/cudnn/lib \
  dockstudio
# /api/engines 中 gnina available=true；POST /api/dock {"engine":"gnina", ...} 实测出 pose + CNN 重打分
```

注意：DockParams 的 `energy_range` 对 gnina 无效（gnina CLI 无此参数，适配器已自动剔除）。

<a id="wsl2"></a>
## Windows 用户：WSL2 跑 GPU 引擎

1. 安装 WSL2：`wsl --install`（管理员 PowerShell），安装 Ubuntu 24.04。
2. 安装 Windows 版最新 NVIDIA 驱动（WSL2 内核自动透传 GPU，WSL 内**不要**装 Linux 驱动）。
3. WSL 里验证：`nvidia-smi` 能看到显卡。
4. 在 WSL 里按上面"服务器模式"跑 Docker（Docker Desktop 开启 WSL2 backend），或直接：
   - 下载 gnina Linux 二进制 → `export DOCKSTUDIO_GNINA=/path/to/gnina`
   - 下载 Uni-Dock Linux 二进制 → `export DOCKSTUDIO_UNIDOCK=/path/to/unidock`
5. WSL 内启动 `py -3.12 -m dockstudio.cli serve`，Windows 浏览器访问 `http://localhost:8000`。
