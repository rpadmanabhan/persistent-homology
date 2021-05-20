import os
import sys

import numpy as np
import pickle

import persistent_homology
import numerics


dump_to_disk = True


def save_data(pik_file, data):
    with open(pik_file, "wb") as f:
        pickle.dump(data, f)


def main(infile):
    '''
    '''
    data = []
    with open(infile, "r") as IN:
        for line in IN:
            contents = line.strip("\n").split(",")
            data.append(list(map(float, contents[1:])))

    data = np.array(data).reshape((len(data), len(contents[1:])))


    vr_complex = persistent_homology.VRComplex(vertices = data, n = 4)
    nrows, ncols = vr_complex.boundary_matrix.shape
    indptr       = vr_complex.boundary_matrix.indptr
    indices      = vr_complex.boundary_matrix.indices
    pivots, pivots_idx = numerics.reduc(nrows, ncols, indptr, indices)

    intervals, homol_ranks, homol = persistent_homology.persistence(
        pivots, pivots_idx,
        vr_complex.simplices_for_lookup,
        vr_complex.weights)

    ## dump to disk - expensive to compute pivots on whole data with no cutoff
    if dump_to_disk:
        pik_file = infile.replace(".csv", ".saved.data")
        save_data(
            pik_file,
            [os.path.basename(infile), vr_complex, dict(pivots), dict(pivots_idx)]
        )

if __name__ == '__main__':
    main(sys.argv[1])
