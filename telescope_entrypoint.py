#!/usr/bin/env python3
"""Shared launcher support for telescope-specific pipeline entry points.

The launchers in this directory deliberately contain no calibration,
flagging, or imaging implementation.  They select a telescope profile and
then delegate to :mod:`pipeline`, whose parser obtains its defaults from the
main ``config.py``/``config.yaml`` configuration.
"""

import sys

import pipeline as PIPELINE


TELESCOPE_ENTRYPOINT_DEFAULTS = {
    # uGMRT is imaged only with CASA tclean. Pinned here, as for LOFAR below,
    # so that an `imager: wsclean` left in a config.yaml shared with a
    # MeerKAT run cannot reach a uGMRT run; asking for WSClean on the command
    # line still fails at preflight rather than being quietly honoured.
    "ugmrt": {"--imager": "casa"},
    # MeerKAT defaults to WSClean through the profile (--imager auto), with
    # its measured katbeam primary beam applied to the result.
    "meerkat": {},
    # The bundled LOFAR path consumes a target MS calibrated upstream (for
    # example by LINC/DP3) and images it with WSClean/EveryBeam.  A beam mode
    # is intentionally not guessed; it remains an explicit CLI/config choice.
    "lofar": {
        "--data-state": "calibrated",
        "--imager": "wsclean",
    },
}


def _option_values(argv, option):
    """Return values explicitly supplied for an argparse-style option."""

    values = []
    index = 0

    while index < len(argv):
        token = str(argv[index])

        if token == "--":
            break

        if token == option:
            if index + 1 < len(argv):
                values.append(str(argv[index + 1]))
                index += 2
                continue

        prefix = option + "="

        if token.startswith(prefix):
            values.append(token[len(prefix):])

        index += 1

    return values


def _has_option(argv, option):
    return bool(_option_values(argv, option))


def parse_for_telescope(telescope, argv=None):
    """Parse the central pipeline CLI with safe instrument defaults."""

    telescope = str(telescope).strip().lower()

    if telescope not in TELESCOPE_ENTRYPOINT_DEFAULTS:
        choices = ", ".join(sorted(TELESCOPE_ENTRYPOINT_DEFAULTS))
        raise ValueError(
            f"Unsupported telescope entry point '{telescope}'. "
            f"Choose one of: {choices}."
        )

    argv = list(sys.argv[1:] if argv is None else argv)
    parser = PIPELINE.build_parser()
    supplied_telescopes = _option_values(argv, "--telescope")

    conflicts = [
        value
        for value in supplied_telescopes
        if value.strip().lower() != telescope
    ]

    if conflicts:
        parser.error(
            f"this launcher is fixed to --telescope {telescope}; "
            f"use pipeline.py for auto-selection or another telescope"
        )

    injected = []

    if not supplied_telescopes:
        # This explicit value intentionally takes precedence over a different
        # telescope left in the shared config from an earlier run.
        injected.extend(["--telescope", telescope])

    for option, value in TELESCOPE_ENTRYPOINT_DEFAULTS[telescope].items():
        if not _has_option(argv, option):
            injected.extend([option, value])

    arguments = parser.parse_args(injected + argv)

    # Retain a final invariant check in case argparse behavior or the central
    # parser changes in the future.
    if arguments.telescope != telescope:
        parser.error(
            f"this launcher requires --telescope {telescope}"
        )

    return arguments


def run_telescope_cli(telescope, argv=None):
    """Run the shared orchestrator for one telescope; return a process code."""

    try:
        arguments = parse_for_telescope(telescope, argv)

        if arguments.write_sbatch is not None:
            PIPELINE.write_sbatch_from_config(arguments)
            return 0

        PIPELINE.PipelineRunner(arguments).run()
        return 0

    except SystemExit:
        raise
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

