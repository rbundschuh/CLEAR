#!/usr/bin/env python3

import argparse
import bisect
import sys
from collections import defaultdict
from typing import Dict, List, Tuple


def main(annotation_path: str, genomecov_path: str) -> None:
    """
    Generate `.dat` files for CLEAR `fitter.py` from `bedtools genomecov` outputs.

    Parameters
    ----------
    annotation_path : str
        Path to the annotation file (genePred format).
    genomecov_path : str
        Path to the profiled BED file (output from `bedtools genomecov`).
    """

    # initialize coverage data structures
    coverage = defaultdict(lambda: [[], [], []])
    coverage_lengths = defaultdict(lambda: 0)

    with open(genomecov_path, 'r') as f:
        for line in f:
            line = line.strip().split()
            if len(line) < 4:
                continue

            coverage[line[0]][0].append(int(line[1]))
            coverage[line[0]][1].append(int(line[2]))
            coverage[line[0]][2].append(int(float(line[3])))

    for key in coverage.keys():
        coverage_lengths[key] = len(coverage[key][0])

    models = []

    with open(annotation_path, 'r') as f:
        for line in f:
            line = line.strip().split()

            if len(line) == 0 or len(line) == 1:
                continue

            # older `genePred` formats had an unused first column (16 total), while
            # newer formats do not (15 total); this hack ensures backwards compatibility
            if len(line) == 15:
                line = ['0'] + line

            seq_name = line[1]
            chr_name = line[2]
            strand = ('+' in line[3])
            name = line[12]
            exons = []

            starts = [int(x) for x in line[9].split(',')[:-1]]
            stops = [int(x) for x in line[10].split(',')[:-1]]

            if not (len(starts) == len(stops)):
                print("START/END MISMATCH!!!", file=sys.stderr)
                sys.exit(0)

            length = 0
            for i in range(len(starts)):
                exon = (starts[i], stops[i])
                length += stops[i] - starts[i]
                exons.append(exon)

            to_add = (seq_name, name, chr_name, strand, exons, length)
            models.append(to_add)

    results = sorted(test_gene(models, coverage))[::-1]

    printed = set()

    # output results
    for a in sorted(results)[::-1]:

        # avoid printing duplicates
        if a[3] in printed:
            continue

        printed.add(a[3])
        a_str = [str(x) for x in a]
        print("\t".join(a_str))


def get_coverage(coverage: Dict[str, List[List[int]]], chromosome: str, locus: int) -> int:
    """
    Get coverage value for a given chromosome and locus.

    Parameters
    ----------
    coverage : Dict[str, List[List[int]]]
        List containing coverage data for each chromosome and locus.
    chromosome : str
        Chromosome name.
    locus : int
        Genomic locus.

    Returns
    -------
    int
        Coverage value at the genomic locus.
    """

    # try lower
    lower = bisect.bisect_left(coverage[chromosome][0], locus) - 1

    try:
        if coverage[chromosome][0][lower] <= locus and coverage[chromosome][1][lower] > locus:
            return coverage[chromosome][2][lower]
        lower += 1
        if coverage[chromosome][0][lower] <= locus and coverage[chromosome][1][lower] > locus:
            return coverage[chromosome][2][lower]
    except IndexError:
        return 0
    return 0


def test_gene(models: List[Tuple[str, str, str, bool, List[Tuple[int, int]], int]], coverage: Dict[str, List[List[int]]]) -> List[Tuple[float, float, int, str, bool]]:
    """
    Test genes and calculate coverage statistics.

    Parameters
    ----------
    coverage : Dict[str, List[List[int]]]
        List containing coverage data for each chromosome and locus.
    models : List[Tuple[str, str, str, bool, List[Tuple[int, int]], int]]
        List of Tuples containing gene information (seq_name, name, chr_name, strand, exons, length).

    Returns
    -------
    List[Tuple[float, float, int, str, bool]]
        List of Tuples containing coverage statistics for each gene.
    """
    results = []
    processed = 0
    for model in models:
        i, total, mu = 0, 0, 0

        for exon in model[4]:
            for j in range(exon[0], exon[1]):
                c = get_coverage(coverage, model[2], j)

                total += c
                mu += c * i
                i += 1

        if total != 0:
            mu = 2 * float(mu) / float(total)
            mu /= float(model[5])
            mu -= 1
        else:
            mu = None

        if mu is not None and not model[3]:  # stranded?
            mu *= -1

        results.append(
            (float(total) / float(model[5]), mu, model[5], model[1], model[3]))

        processed += 1
        if processed % 1000 == 0:
            print(
                f"{processed} processed, {len(models)-processed} to go! ( {100*float(processed)/len(models)}% )",
                file=sys.stderr
            )

    return results


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Generate `.dat` files for CLEAR `fitter.py` using `bedtools genomecov` outputs."
    )

    parser.add_argument(
        "-a", "--annotation",
        type=str,
        help="Path to the annotation (`genePred` format)."
    )

    parser.add_argument(
        "-b", "--bed",
        type=str,
        help="Path to the profiled BED file (output from `bedtools genomecov`)."
    )

    args = parser.parse_args()

    main(args.annotation, args.bed)
