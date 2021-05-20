

import numpy as np
import numba

@numba.njit
def symm_diff(arr1, arr2):
    ''' Compute the symmetric difference between two sorted numpy arrays
    i.e. return non-intersecting elements in sorted order
    '''

    # Traverse both arrays
    # simultaneously.
    i = 0
    j = 0
    out = []
    while i < arr1.size and j < arr2.size:

        ## Check which array has the smaller element
        ## add it to our out list
        if (arr1[i] < arr2[j]):
            out.append(arr1[i])
            i += 1
        elif (arr2[j] < arr1[i]):
            out.append(arr2[j])
            j += 1
        ## don't add if duplicate
        else:
            i += 1
            j += 1

    # Add remaining element of either array
    while i < arr1.size:
        out.append(arr1[i])
        i += 1

    while j < arr2.size:
        out.append(arr2[j])
        j += 1

    return np.array(out)


int_array = numba.types.int32[:]
@numba.njit
def reduc(nrows, ncols, indptr, indices):
    '''
    '''
    ## Track any column with a pivot_row
    ## pivot_row_idx -> col_value
    pivots = numba.typed.Dict.empty(
        key_type   = numba.types.int64,
        value_type = int_array,
    )
    ## pivot_row_idx -> col_idx
    pivots_idx = numba.typed.Dict.empty(
        key_type   = numba.types.int64,
        value_type = numba.types.int64,
    )
    ## iterate over cols
    for j in range(0, ncols):
        row_idx_j = indices[indptr[j]: indptr[j + 1]]

        ## Print progress
        if j % 10000 == 0:
            print(j)

        ## Reduce this column till there is a pivot and some k < j
        ## has a pivot on the same row
        while len(row_idx_j) != 0:
            if row_idx_j[-1] not in pivots_idx: break
            row_idx_k = pivots[row_idx_j[-1]]
            row_idx_j = symm_diff(row_idx_j, row_idx_k)

        ## Keep track of updated non-zero columns for future
        if len(row_idx_j) != 0:
            pivots[row_idx_j[-1]] = row_idx_j
            pivots_idx[row_idx_j[-1]] = j

    return pivots, pivots_idx

