#!/usr/bin/env python3
"""Launch the shared radio-continuum pipeline with the uGMRT profile."""

from telescope_entrypoint import run_telescope_cli


def main(argv=None):
    return run_telescope_cli("ugmrt", argv)


if __name__ == "__main__":
    raise SystemExit(main())

