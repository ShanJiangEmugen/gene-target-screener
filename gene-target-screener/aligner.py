# u7_target_screener/aligner.py

from __future__ import annotations

from typing import Dict, List, Tuple, Union

import numpy as np
import pandas as pd
import re
from Bio import Align
from Bio.Align import substitution_matrices
from Bio.Seq import Seq


SeqLike = Union[str, Seq]


class SubsequenceAligner:
    """
    Extract subsequences from an input sequence and align them to a reference
    using Biopython's PairwiseAligner.

    Typical use:
        aligner = SubsequenceAligner(input_seq, ref_seq)
        aligner.make_query(seq_low=18, seq_high=36)
        aligner.init_aligner()
        aligner.get_alignment(threshold=0.8)
        df = aligner.results
    """

    def __init__(self, input_seq: SeqLike, ref_seq: SeqLike) -> None:
        # Normalize to Seq objects and upper-case
        self.input_seq: Seq = Seq(str(input_seq).upper())
        self.ref_seq: Seq = Seq(str(ref_seq).upper())
        self.query: Dict[str, List[Tuple[Seq, str]]] = {}
        self.aligner: Align.PairwiseAligner | None = None
        self.results: pd.DataFrame | None = None

    def make_query(self, seq_low: int = 18, seq_high: int = 36) -> None:
        """
        Slide a window across input_seq to generate candidate subsequences.

        Parameters
        ----------
        seq_low : int
            Minimum subsequence length (inclusive).
        seq_high : int
            Maximum subsequence length (exclusive) — Python-style range.
        """
        query: Dict[str, List[Tuple[Seq, str]]] = {}
        total = 0

        for length in range(seq_low, seq_high):
            temp: List[Tuple[Seq, str]] = []
            for i in range(len(self.input_seq) - length + 1):
                subseq = self.input_seq[i : i + length]

                # Filter out CTGCTGCTG repeat-rich segments
                ctg_hits = [
                    m.start()
                    for m in re.finditer("CTGCTGCTG", str(subseq))
                ]
                if len(ctg_hits) > 0:
                    continue

                # Save subseq and "start/total_length" annotation
                temp.append((subseq, f"{i+1}/{len(self.input_seq)}"))
                total += 1

            query[str(length)] = temp

        print(f"Total {total} sub-sequences were extracted!")
        self.query = query

    def init_aligner(
        self,
        mode: str = "global",
        open_gap: float = -10.0,
        extend_gap: float = -1.0,
        left_open: float = 0.0,
        left_extend: float = 0.0,
        right_open: float = 0.0,
        right_extend: float = 0.0,
        matrix: str = "BLASTN",
    ) -> None:
        """
        Initialize Biopython PairwiseAligner with scoring parameters.

        Notes
        -----
        - Default substitution matrix 'BLASTN':
          match = +2, mismatch = -3
        """
        aligner = Align.PairwiseAligner(mode=mode)
        aligner.open_gap_score = open_gap
        aligner.extend_gap_score = extend_gap

        # Boundary gap penalties (fixed version)
        aligner.left_open_gap_score = left_open
        aligner.left_extend_gap_score = left_extend
        aligner.right_open_gap_score = right_open
        aligner.right_extend_gap_score = right_extend

        # Nucleotide scoring matrix
        aligner.substitution_matrix = substitution_matrices.load(matrix)

        self.aligner = aligner
        print("Aligner initialized")

    def get_alignment(self, threshold: float = 0.8) -> None:
        """
        Align all generated subsequences to the reference and collect hits.

        Parameters
        ----------
        threshold : float
            Minimum fraction of matches (diagonal substitutions / length)
            to keep a subsequence.

        Result
        ------
        self.results : pandas.DataFrame
            Columns:
                - Position
                - Sequence
                - Length
                - Matched
                - Conservation
                - Score
                - Score Max
                - Score-wise Conservation
        """
        if self.aligner is None:
            raise RuntimeError(
                "Aligner not initialized. Call init_aligner() first."
            )
        if not self.query:
            raise RuntimeError(
                "No subsequences generated. Call make_query() first."
            )

        result = pd.DataFrame(
            columns=[
                "Position",
                "Sequence",
                "Length",
                "Matched",
                "Conservation",
                "Score",
                "Score Max",
                "Score-wise Conservation",
            ]
        )

        count = 0
        total = 0

        for length_key, subseqs in self.query.items():
            for subseq, pos_str in subseqs:
                total += 1
                alignments = self.aligner.align(self.ref_seq, subseq)

                for alignment in alignments:
                    # Number of matches from substitution matrix diagonal
                    matches = np.diagonal(alignment.substitutions).sum()
                    length = len(subseq)

                    # Fraction of conserved positions
                    conservation = matches / length

                    # Simple filters:
                    # 1) conservation >= threshold
                    # 2) score at least as large as length
                    if (conservation >= threshold) and (
                        alignment.score >= length
                    ):
                        score_max = length * 2  # BLASTN: match = +2
                        scorewise_conservation = alignment.score / score_max

                        result.loc[count] = [
                            pos_str.split("/")[0],
                            str(subseq),
                            length,
                            matches,
                            conservation,
                            alignment.score,
                            score_max,
                            scorewise_conservation,
                        ]
                        count += 1

        print(f"Got {count} subsequence(s) from {total} segments")
        self.results = result

