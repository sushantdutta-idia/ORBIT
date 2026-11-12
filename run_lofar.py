#!/usr/bin/env python3
"""Launch the calibrated-target LOFAR continuum imaging pipeline."""

from telescope_entrypoint import run_telescope_cli


def main(argv=None):
    return run_telescope_cli("lofar", argv)


if __name__ == "__main__":
    raise SystemExit(main())

