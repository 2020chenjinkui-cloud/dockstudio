# DockStudio backend (server mode). Build from repo root:
#   docker build -f docker/Dockerfile -t dockstudio .
FROM python:3.12-slim

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1

# vina needs libgomp at runtime
RUN apt-get update \
    && apt-get install -y --no-install-recommends curl libgomp1 \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/dockstudio
COPY server ./server
COPY resources/bin/linux-x86_64 ./resources/bin/linux-x86_64

RUN pip install "./server" \
    && chmod +x ./resources/bin/linux-x86_64/vina || true

ENV DOCKSTUDIO_DATA_DIR=/data
RUN useradd -m dockstudio \
    && mkdir -p /data \
    && chown -R dockstudio:dockstudio /data /opt/dockstudio
USER dockstudio
VOLUME ["/data"]

EXPOSE 8000
HEALTHCHECK --interval=30s --timeout=5s --start-period=10s --retries=3 \
    CMD curl -fsS http://127.0.0.1:8000/api/health || exit 1

CMD ["python", "-m", "uvicorn", "dockstudio.api.app:create_app", "--factory", "--host", "0.0.0.0", "--port", "8000"]
