"""Frozen 打包时的后端启动入口：PyInstaller spec 以此为主脚本。"""

from dockstudio.api.app import main

if __name__ == "__main__":
    main()
