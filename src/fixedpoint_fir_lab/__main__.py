"""Run with python -m fixedpoint_fir_lab experiment --output results/demo.json."""

import argparse
import json

from .experiment import DEFAULT_SEED, save_experiment


def main(argv=None):
    parser = argparse.ArgumentParser(description="Standard-library synthetic fixed-point FIR experiment")
    subcommands = parser.add_subparsers(dest="command", required=True)
    experiment = subcommands.add_parser("experiment", help="write numerical evidence and synthetic fixtures")
    experiment.add_argument("--output", required=True, help="JSON report path")
    experiment.add_argument("--seed", type=int, default=DEFAULT_SEED)
    experiment.add_argument("--samples", type=int, default=512)
    args = parser.parse_args(argv)
    try:
        report = save_experiment(args.output, seed=args.seed, samples=args.samples)
    except (ValueError, OverflowError, OSError) as error:
        parser.error(str(error))
    print(json.dumps({"complete": report["experiment_complete"], "records": len(report["cases"]),
                      "outputs_per_record": report["cases"][0]["output_samples"],
                      "int32_bound": report["accumulator"]["conservative_abs_bound"],
                      "shifts": report["parameters"]["output_shifts"]}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
