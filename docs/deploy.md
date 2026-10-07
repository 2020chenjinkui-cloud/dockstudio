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

<a id="wsl2"></a>
## Windows 用户：WSL2 跑 GPU 引擎

1. 安装 WSL2：`wsl --install`（管理员 PowerShell），安装 Ubuntu 24.04。
2. 安装 Windows 版最新 NVIDIA 驱动（WSL2 内核自动透传 GPU，WSL 内**不要**装 Linux 驱动）。
3. WSL 里验证：`nvidia-smi` 能看到显卡。
4. 在 WSL 里按上面"服务器模式"跑 Docker（Docker Desktop 开启 WSL2 backend），或直接：
   - 下载 gnina Linux 二进制 → `export DOCKSTUDIO_GNINA=/path/to/gnina`
   - 下载 Uni-Dock Linux 二进制 → `export DOCKSTUDIO_UNIDOCK=/path/to/unidock`
5. WSL 内启动 `py -3.12 -m dockstudio.cli serve`，Windows 浏览器访问 `http://localhost:8000`。
