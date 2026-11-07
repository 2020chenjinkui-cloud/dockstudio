"""Entrypoint when running `python -m dockstudio.cli`."""

import sys

from . import main

if __name__ == "__main__":
    sys.exit(main())
