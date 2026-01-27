#!/usr/bin/env python3

import argparse
import math
import os
import sys
from glob import glob
from typing import List

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages

from fitter import find_state_change, load_file, sub_file

matplotlib.use('Agg')





def main(inputs: List[str], output: str, page_width: int, page_height: int) -> None:
    """
    Generate scatter plots from ``.dat`` files and write them to a PDF.

    Parameters
    ----------
    inputs : List[str]
        One or more glob patterns for input ``.dat`` files.
    output : str
        Output PDF filename.
    page_width : int
        Number of plots per row.
    page_height : int
        Number of rows per page.
    """

    # collect files from provided glob patterns
    files = []
    for pattern in inputs:
        files.extend(glob(pattern))
    files = sorted(set(files))

    if not files:
        print("No input .dat files matched the given patterns.", file=sys.stderr)
        sys.exit(1)

    # load inputs
    file_data = [load_file(f) for f in files]
    file_data_sub = [sub_file(d) for d in file_data]
    limits = [find_state_change(d) for d in file_data]

    pdf = PdfPages(output)

    per_page = page_width * page_height
    page_count = math.ceil(len(files) / float(per_page))

    for page_num in range(int(page_count)):
        fig, axes = plt.subplots(
            nrows=page_height, ncols=page_width, figsize=(14, 7), dpi=300
        )

        for i in range(page_width * page_height):
            dat_index = page_num * page_width * page_height + i

            x_cord = i % page_width
            y_cord = i // page_width

            ax = axes[y_cord, x_cord]

            if dat_index >= len(files):
                ax.axis('off')
                ax.get_xaxis().set_visible(False)
                ax.get_yaxis().set_visible(False)
                continue

            ax.set_title(os.path.basename(
                files[dat_index]).replace('.dat', ''))

            ys = np.asarray(file_data_sub[dat_index], dtype=float)
            xs = np.arange(len(ys), dtype=int)
            ax.scatter(
                xs,
                ys,
                s=2,
                alpha=0.25,
                linewidths=0,
                color='k',
                rasterized=True,
            )

            if 0 <= limits[dat_index] < len(ys):
                ax.axvline(x=limits[dat_index], color='r', linestyle='-', linewidth=2)

            RIG_FACTOR = 0.1
            ax.set_ylim(-1 - RIG_FACTOR, 1 + RIG_FACTOR)

            ax.set_xlim(-1, max(1, len(ys)))

        fig.tight_layout()
        pdf.savefig(fig)

    pdf.close()


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Generate CLEAR scatter plots from .dat files."
    )

    parser.add_argument(
        "-i",
        "--inputs",
        nargs="+",
        default=["*.dat"],
        help="One or more glob patterns for input .dat files (default: *.dat)",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="CLEAR_scatter.pdf",
        help="Output PDF filename (default: CLEAR_scatter.pdf)",
    )
    parser.add_argument(
        "--page-width",
        type=int,
        default=3,
        help="Number of plots per row (default: 3)",
    )
    parser.add_argument(
        "--page-height",
        type=int,
        default=2,
        help="Number of rows per page (default: 2)",
    )

    args = parser.parse_args()

    main(args.inputs, args.output, args.page_width, args.page_height)
