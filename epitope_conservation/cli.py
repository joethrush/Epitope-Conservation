import argparse
import sys

from .pipeline import run_pipeline


def build_parser():
    """Define the command-line arguments accepted by the analysis tool."""
    parser = argparse.ArgumentParser(
        prog="epitope-conservation",
        description="Rank conserved, accessible epitope candidates between two protein structures.",
    )
    parser.add_argument("state1", help="First PDB or mmCIF structure")
    parser.add_argument("state2", help="Second PDB or mmCIF structure")
    parser.add_argument("--chains1", required=True, help="Comma-separated target chain IDs in state 1")
    parser.add_argument("--chains2", required=True, help="Comma-separated target chain IDs in state 2")
    parser.add_argument("--output", "-o", default="epitope_output", help="Output directory")
    parser.add_argument("--neighborhood-size", type=int, default=8)
    parser.add_argument("--distance", type=float, default=8.0)
    parser.add_argument("--rmsd-scale", type=float, default=4.0)
    parser.add_argument("--weight-backbone", type=float, default=0.60)
    parser.add_argument("--weight-neighborhood", type=float, default=0.20)
    parser.add_argument("--weight-exposure", type=float, default=0.10)
    parser.add_argument("--weight-glycan", type=float, default=0.10)
    parser.add_argument("--force-as-is", action="store_true",
                        help="Override assembly metadata and analyze provided chains as-is")
    return parser


def main(argv=None):
    """Validate command-line input, run the analysis, and report user-facing errors."""
    args = build_parser().parse_args(argv)
    chains1 = [value.strip() for value in args.chains1.split(",") if value.strip()]
    chains2 = [value.strip() for value in args.chains2.split(",") if value.strip()]
    if not chains1 or not chains2:
        raise SystemExit("Both --chains1 and --chains2 must contain at least one chain ID")
    weights = (args.weight_backbone, args.weight_neighborhood,
               args.weight_exposure, args.weight_glycan)
    if sum(weights) <= 0:
        raise SystemExit("Score weights must have a positive sum")
    try:
        run_pipeline(
            state1=args.state1,
            state2=args.state2,
            chains1=chains1,
            chains2=chains2,
            output_dir=args.output,
            neighborhood_size=args.neighborhood_size,
            distance=args.distance,
            rmsd_scale=args.rmsd_scale,
            weights=weights,
            force_as_is=args.force_as_is,
        )
    except (ValueError, KeyError, IOError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
