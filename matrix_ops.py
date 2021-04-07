import numpy as np

# Team : Bill Lee, Raghavendra Padmanabhan and Francisco Vargas
#--------------------------------------------------------------------------------------

def row_echelon_Z2(A, bounds = None, reduce_cols = False):
    ''' Inplace reduce a Z mod 2 matrix to row echelon form
    :param np.array(dtype = np.uint8): A
    Can also operate on a specific Block of a matrix.
    :param np.array(dtype = np.uint8) array: A
    :param tuple: bounds: (start_row, end_row, start_col, end_col) (default: (1, nrows, 1, ncols))
    :param bool: reduce_cols: Use column operations to 0 out cols to the right of a pivot
    --------------------------------------
    :returns rank of the matrix A
    :rtype int
    '''
    nrows, ncols = A.shape
    if bounds is not None:
        start_row, end_row, start_col, end_col = bounds
    else:
        start_row, end_row, start_col, end_col = 1, nrows, 1, ncols

    cur_row = start_row - 1
    for j in range(start_col - 1, end_col):
        ## Find pivot row, can go with the first row having a 1 since everything is either 1 or 0
        pivot_row = None
        for i in range(cur_row, end_row):
            if A[i][j]:
                pivot_row = i
                break

        ## all 0s => No pivots
        if pivot_row is None: continue

        ## swap pivot_row with cur_row if necessary
        if cur_row != pivot_row:
            A[[cur_row, pivot_row], :] = A[[pivot_row, cur_row], :]
            pivot_row = cur_row

        ## Next pivot must start here
        cur_row += 1

        ## 0 out rows below the pivot
        for i in range(pivot_row + 1, end_row):
            if A[i][j]:
                ## Use XOR for Z mod 2 arithmetic
                A[i] = A[i] ^ A[pivot_row]

        if reduce_cols:
            ## 0 out cols to the right of the pivot
            for k in range(j + 1, end_col):
                if A[cur_row - 1][k]:
                    A[:, k] = A[:, k] ^ A[:, j]

    return cur_row


def reduced_row_echelon_Z2(A, bounds = None):
    ''' Inplace reduce a Z mod 2 matrix to row reduced echelon form.
    Can also operate on a specific Block of a matrix.
    :param np.array(dtype = np.uint8) array: A
    :param tuple: bounds: (start_row, end_row, start_col, end_col) (default: (1, nrows, 1, ncols))
    --------------------------------------
    :returns rank of the matrix A
    :rtype int
    '''
    nrows, ncols = A.shape
    if bounds is not None:
        start_row, end_row, start_col, end_col = bounds
    else:
        start_row, end_row, start_col, end_col = 1, nrows, 1, ncols

    ## First reduce the matrix to its row echelon form
    rank = row_echelon_Z2(A, (start_row, end_row, start_col, end_col))

    ## Now to reduced row echelon form
    cur_col = end_col - 1
    ## Start from the last non-zero row
    for i in range(rank - 1, -1, -1):
        ## Find pivot column
        pivot_col = None
        for j in range(start_col - 1, end_col):
            if A[i, j]:
                pivot_col = j
                break
        assert pivot_col is not None, "Bug ! Could not find pivot in row: {}\nrow_idx: {}".format(A[i], i)

        ## Zero out all elements above this column
        for k in range(0, i):
            if A[k, j]:
                A[k] = A[k] ^ A[i]

    return rank


def inverse_Z2(A):
    ''' Return inverse of a matrix in Z mod 2
    :param np.array(dtype = np.uint8) array: A
    --------------------------------------
    :returns inverse of A
    :rtype np.array(dtype = np.uint8) array: A
    '''
    nrows, ncols = A.shape
    assert nrows == ncols, "A is not square !"
    ## Augment the matrix with identity
    A_aug = np.hstack([A, np.eye(nrows, dtype = np.uint8)])
    ## reduced row echelon form
    reduced_row_echelon_Z2(A_aug, (1, nrows, 1, ncols))
    ## return 2nd block
    return A_aug[:, nrows:]


def smith_normal_form(A):
    ''' Compute the smith normal form of a matrix A
    :param np.array(dtype = np.uint8) array: A
    ---------------------------------------
    :returns: Decomposes A as : P^-1 D A = A Q  and returns (rank_A, P^-1, D, Q)
    :rtype: tuple of (rank_A, np.array(dtype=np.uint8), np.array(dtype=np.uint8), np.array(dtype=np.uint8))
    '''

    ## Followed algorithm from here :
    ## http://people.maths.ox.ac.uk/nanda/cat/Lecture%2003%20Homology.pdf
    #--------------------------------------------------------------------------------------
    ## B = I_nxn  A
    ##     0_mxn  I_mxm
    ## Note: 0_mxn will never get modified.
    nrows, ncols = A.shape
    B = np.block(
        [[np.eye(nrows, dtype = np.uint8), A],
         [np.zeros(shape = (ncols, nrows), dtype = np.uint8),
          np.eye(ncols, dtype = np.uint8)]])
    ## Use Row/Col ops to transform A to a diagonal matrix
    # 1. Reduce it to echelon form, use column ops to reduce right of pivot cols to 0
    rank_A = row_echelon_Z2(B, (1, nrows, nrows + 1, nrows + ncols), reduce_cols = True)
    # 2. Swap pivot cols with non pivot ones to get a diagonal matrix
    prev_pivot_col = None
    for j in range(nrows, nrows + ncols):
        for i in range(0, nrows):
            if not B[i][j]: continue
            if i != j and prev_pivot_col is not None:
                # swap with col next to previous i=j pivot
                B[:, [prev_pivot_col + 1, j]] = B[:, [j, prev_pivot_col + 1]]
                prev_pivot_col += 1
            else:
                prev_pivot_col = j
            break

    P = B[0: nrows, 0: nrows]
    Q = B[nrows:, nrows:]
    D = B[0: nrows:, nrows:]
    ## Note:
    ## basis for im(A) corresponds to the first `rank_A` cols of P_inv
    ## basis for ker(A) corresponds to the last `(n - rank_A)` cols of Q
    return (rank_A, inverse_Z2(P), D, Q)


def basis_img_ker_Z2(d):
    ''' Compute basis for the image/kernel of a boundary matrix
    '''
    nrows, ncols = d.shape
    ## Smith Normal Form for the boundary matrix
    rank, P_inv, D, Q = smith_normal_form(d)
    ## return basis for image - first rank cols of P_inv and kernel - last nrows - rank cols of Q
    return P_inv[:, :rank], Q[:, rank:]


def basis_homology_Z2(d_k, d_k1):
    ''' Compute basis for the kth homology group
    :param d_k: Boundary matrix for C_k
    :param d_k1: Boundary matrix for C_k+1
    '''
    nrows_k, ncols_k = d_k.shape

    ## Form a block matrix with basis vectors of the img of C_k+1 and ker of C_k
    ## G_k = [B_k | Z_k]
    img_d_k, ker_d_k = basis_img_ker_Z2(d_k)
    img_d_k1, ker_d_k1 = basis_img_ker_Z2(d_k1)
    rank_ker_d_k = ker_d_k.shape[1]
    rank_img_d_k1 = img_d_k1.shape[1]
    if rank_img_d_k1 == 0:
        img_d_k1 = np.array(
            img_d_k1, dtype=np.uint8).reshape(ker_d_k.shape[0], 0)
    G_k = np.hstack([img_d_k1, ker_d_k])
    ## Compute reduced row echelon form of G_k
    reduced_row_echelon_Z2(G_k)
    ## E_k = reduced_row_echelon(G_k) = [B'_k | Z'_k]

    ## Find pivot columns of Z'_k
    ## Start from the bottom row of E_k, the last rank_k - rank_k1 pivot cols is our desired answer
    nrows_G_k, ncols_G_k = G_k.shape
    pivot_cols = []
    for i in range(nrows_G_k - 1, -1, -1):
        ## Find the pivot column
        for j in range(0, ncols_G_k):
            if G_k[i, j]:
                ## We need the col number relative to ker_d_k sub-block
                pivot_col = j - rank_img_d_k1
                assert pivot_col >= 0
                pivot_cols.append(pivot_col)
                break
        if len(pivot_cols) == rank_ker_d_k - rank_img_d_k1: break

    ## Corresponding columns in Z_k form a basis for H_k(C_k, d_k)
    return ker_d_k[:, pivot_cols]
