import glob
import pickle
import os
import sys

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

import persistent_homology



def load_data(pik_file):
    '''
    '''
    with open(pik_file, "rb") as IN:
        data = pickle.load(IN)

    return data



def make_latex_table(data):
    '''
    '''
    table = """ \\begin{center}
    \\begin{tabular}{|| c c c ||}
    \\hline
    Birth & Death & Basis \\\\
    \\hline\\hline
    """
    for row in data:
        table += " & ".join(row)
        table += " \\\\"
        table += "\n"
        table += "\\hline"
    table += "\n"
    table += "\\end{tabular}\n"
    table += "\\end{center}"

    return table


def birth_death_plot_many(many_intervals):
    ''' Birth Death Plots
    '''

    fig, axs = plt.subplots(3, 4, figsize = (15, 15))

    labels = ["H0", "H1", "H2"]
    row = 0
    col = 0
    lims = [
        -1,  # min of both axes
        12,  # max of both axes
    ]
    for title, intervals in many_intervals:
        max_weight = 0
        for key in intervals:
            birth = [e[0] for e in intervals[key]]
            death = [e[1] for e in intervals[key]]
            max_weight = max(max_weight, max(death))
            axs[row, col].scatter(birth, death, alpha = 0.5, s = 20, label = labels[key])
        axs[row, col].set_xlabel("Birth")
        axs[row, col].set_ylabel("Death")
        # now plot both limits against each other
        axs[row, col].plot(lims, lims, 'k-', alpha=0.75, zorder=0)
        axs[row, col].set_aspect('equal')
        axs[row, col].set_xlim(lims)
        axs[row, col].set_ylim(lims)
        axs[row, col].set_title("{}".format(title) , fontweight = "bold")
        axs[row, col].legend(loc = "lower right")
        # max weight
        axs[row, col].axhline(y = max_weight)
        axs[row, col].axvline(x = max_weight)

        col += 1
        if col == 4:
            row += 1
            col = 0
        if row == 3:
            row = 0

    plt.suptitle("Birth Death Diagram", fontweight = "bold")
    plt.tight_layout()
    plt.savefig("birth_death.png")



def main(out_dir):
    '''
    '''

    all_intervals = []
    ## load picked data from disk
    for f in sorted(glob.glob(os.path.join(out_dir, "*.data"))):
        print(f)
        data = load_data(f)
        ## retrieve important datastructures
        fname, vr_complex, pivots, pivots_idx = data

        intervals, homol_ranks, homol = persistent_homology.persistence(
            pivots, pivots_idx,
            vr_complex.simplices_for_lookup,
            vr_complex.weights)

        sample = os.path.basename(fname).replace(".csv", "")
        all_intervals.append((sample, intervals))

    birth_death_plot_many(all_intervals)

    for sample, intervals in all_intervals:
        print(sample)
        persistent_homologies = persistent_homology.check_for_persistent_homologies(intervals, 10)
        homologies = list(((str(round(p[1], 2)), str(round(p[2], 2)), str(p[3])) for p in persistent_homologies))
        print(len(homologies))
        print(make_latex_table(homologies))
        print("-------------------------------")

if __name__ == '__main__':
    main(sys.argv[1])
