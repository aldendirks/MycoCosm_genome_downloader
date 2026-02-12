import argparse
from . import prepare, download


def main():
    parser = argparse.ArgumentParser(
        prog="mycocosm_downloader",
        description="MycoCosm genome downloader"
    )

    subparsers = parser.add_subparsers(dest="command", required=True)

    # ---- prepare ----
    prep_parser = subparsers.add_parser("prepare", help="Prepare metadata")
    prep_parser.add_argument("--getcsv", action="store_true")
    prep_parser.add_argument("--getxml", action="store_true")

    # ---- download ----
    dl_parser = subparsers.add_parser("download", help="Download genomes")
    dl_parser.add_argument("--outdir", default="data")

    args = parser.parse_args()

    if args.command == "prepare":
        prepare.run(args)

    elif args.command == "download":
        download.run(args)