#!/usr/bin/env python3

import argparse
import math
import os
import sys
from glob import glob
from typing import List

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from fitter import find_state_change, load_file, sub_file

matplotlib.use('Agg')





def main(inputs: List[str], output: str, page_width: int, page_height: int) -> None:
    """
    Generate violin plots from ``.dat`` files and write them to a PDF.

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

            pos, data = [], []
            x_max = 15
            for j in range(0, x_max):
                inc = 1000
                p = inc * j
                pos.append(p + (inc / 2))
                a = file_data_sub[dat_index][p:p+inc]
                if not a:
                    a = [0]
                data.append(a)

            # draw the violin plot, note that the bw_method parameter defines the 'definition' of the KDE
            ax.violinplot(
                data,
                pos,
                widths=850,
                showmeans=True,
                showextrema=True,
                showmedians=False,
                bw_method=0.2,
            )
            del pos, data

            if limits[dat_index] < x_max * 1000:
                ax.axvline(x=limits[dat_index], color='r',
                           linestyle='-', linewidth=2)

            RIG_FACTOR = 0.1
            ax.set_ylim(-1 - RIG_FACTOR, 1 + RIG_FACTOR)

        fig.tight_layout()
        pdf.savefig(fig)

    pdf.close()


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Generate CLEAR violin plots from .dat files."
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
        default="CLEAR_violins.pdf",
        help="Output PDF filename (default: CLEAR_violins.pdf)",
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
