# u7_target_screener/batch_run.py

from __future__ import annotations

import argparse
import os
from pathlib import Path

import pandas as pd
from Bio.Seq import Seq

from .aligner import SubsequenceAligner


def run_batch(
    csv_path: str,
    out_dir: str = "screen2",
    seq_low: int = 18,
    seq_high: int = 36,
    threshold: float = 0.8,
) -> None:
    """
    Run pairwise U7-target screening based on a metadata CSV.

    Expected CSV columns:
        - "Comparison #"
        - "Region"
        - "Sequence"
    """
    csv_path = Path(csv_path)
    out_root = Path(out_dir)
    out_root.mkdir(parents=True, exist_ok=True)

    meta_data = pd.read_csv(csv_path)

    # 假设每个 comparison 有两行：一条作为 ref，一条作为 input
    n_pairs = meta_data["Comparison #"].nunique()
    print(f"Found {n_pairs} comparison pair(s)")

    for i in range(n_pairs):
        pair_id = i + 1
        print(f"\nWorking on Comparison #{pair_id}")

        temp_df = meta_data[meta_data["Comparison #"] == pair_id]

        if len(temp_df) != 2:
            print(
                f"  [Warning] Comparison #{pair_id} has {len(temp_df)} rows, "
                f"expected 2. Skipping."
            )
            continue

        # 第一行作为 ref_seq
        ref_seq = Seq(str(temp_df["Sequence"].iloc[0]).upper())
        ref_name = str(temp_df["Region"].iloc[0])

        # 第二行作为 input_seq
        inp_seq = Seq(str(temp_df["Sequence"].iloc[1]).upper())
        inp_name = str(temp_df["Region"].iloc[1])

        pair_name = f"pair{pair_id}_{inp_name}_{ref_name}"
        pair_dir = out_root / pair_name
        pair_dir.mkdir(parents=True, exist_ok=True)

        # Direction 1: input -> ref
        print(f"  Direction 1: {inp_name} (input) -> {ref_name} (ref)")
        align1 = SubsequenceAligner(inp_seq, ref_seq)
        align1.make_query(seq_low=seq_low, seq_high=seq_high)
        align1.init_aligner()
        align1.get_alignment(threshold=threshold)

        out_csv_1 = pair_dir / f"input_{inp_name}_refseq_{ref_name}.csv"
        align1.results.to_csv(out_csv_1, index=False)

        # Direction 2: ref -> input
        print(f"  Direction 2: {ref_name} (input) -> {inp_name} (ref)")
        align2 = SubsequenceAligner(ref_seq, inp_seq)
        align2.make_query(seq_low=seq_low, seq_high=seq_high)
        align2.init_aligner()
        align2.get_alignment(threshold=threshold)

        out_csv_2 = pair_dir / f"input_{ref_name}_refseq_{inp_name}.csv"
        align2.results.to_csv(out_csv_2, index=False)

    print("\nDone.")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Batch U7 antisense target screening."
    )
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Input CSV file (metadata with sequences).",
    )
    parser.add_argument(
        "-o",
        "--out-dir",
        default="screen2",
        help="Output directory root (default: screen2).",
    )
    parser.add_argument(
        "--seq-low",
        type=int,
        default=18,
        help="Minimum subsequence length (inclusive).",
    )
    parser.add_argument(
        "--seq-high",
        type=int,
        default=36,
        help="Maximum subsequence length (exclusive).",
    )
    parser.add_argument(
        "-t",
        "--threshold",
        type=float,
        default=0.8,
        help="Conservation threshold (default: 0.8).",
    )

    args = parser.parse_args()
    run_batch(
        csv_path=args.input,
        out_dir=args.out_dir,
        seq_low=args.seq_low,
        seq_high=args.seq_high,
        threshold=args.threshold,
    )


if __name__ == "__main__":
    main()

