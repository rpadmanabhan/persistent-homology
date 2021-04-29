# stdlib
import collections
import itertools
import functools
import math
import operator

# 3rd party
import numpy as np
import scipy
import scipy.special
import scipy.spatial
import scipy.sparse
import sortednp
import numba
import matplotlib
import matplotlib.pyplot as plt

# ours
import numerics


def birth_death_plot(intervals, plotname):
    ''' Birth Death Plots
    '''

    fig, axs = plt.subplots(figsize = (7, 7))
    labels = ["H0", "H1", "H2"]
    for key in intervals:
        birth = [e[0] for e in intervals[key]]
        death = [e[1] for e in intervals[key]]
        axs.scatter(birth, death, alpha = 0.2, s = 10, label = labels[key])
    axs.set_xlabel("Birth")
    axs.set_ylabel("Death")

    lims = [
        -1,  # min of both axes
        np.max([axs.get_xlim(), axs.get_ylim()]),  # max of both axes
    ]

    # now plot both limits against each other
    axs.plot(lims, lims, 'k-', alpha=0.75, zorder=0)
    axs.set_aspect('equal')
    axs.set_xlim(lims)
    axs.set_ylim(lims)
    axs.set_title("Birth Death Diagram" , fontweight = "bold")
    axs.legend(loc = "best")

    plt.tight_layout()
    plt.savefig(plotname)



def persistence(pivots, pivots_idx, simplices, weights, homol_dim_cutoff = 2):
    ''' Compute intervals for homologies
    '''

    ## bookeeping
    non_zero_cols  = set()
    homol_ranks    = collections.defaultdict(int)
    homol          = collections.defaultdict(list)
    intervals      = collections.defaultdict(list)

    ## process all pivot rows
    for pivot_row_idx, pivot_col_idx in pivots_idx.items():

        non_zero_cols.add(pivot_col_idx)

        homol_basis    = simplices[pivot_row_idx]
        homol_dim      = len(homol_basis) - 1
        homol[homol_dim].append((homol_basis, pivot_row_idx, pivot_col_idx))
        homol_ranks[homol_dim] += 1
        homol_interval = (
            weights[pivot_row_idx], weights[pivot_col_idx],
            pivot_row_idx, pivot_col_idx, homol_basis)

        intervals[homol_dim].append(homol_interval)


    ## check all zero cols
    for col in range(len(simplices)):
        if col in non_zero_cols: continue
        if col not in pivots_idx:
            homol_basis  = simplices[col]
            ## Skip these,
            ## we don't have higher dimensional simplices in the filtration boundary matrix
            ## to know if they ever die
            if len(homol_basis) - 1 > homol_dim_cutoff: continue

            homol_dim      = len(homol_basis) - 1
            homol[homol_dim].append((homol_basis, col, col))
            homol_ranks[homol_dim] += 1
            homol_interval = (
                weights[col], weights[-1],
                col, -1, homol_basis)

            intervals[homol_dim].append(homol_interval)


    return intervals, homol_ranks, homol


class VRComplex:
    '''
    '''
    def __init__(self, *args, **kwargs):
        '''
        '''
        ## Max radius ball
        self.epsilon = kwargs.get("epsilon", None)
        ## Max dim simplex
        self.max_dim = kwargs["n"]
        ## The input data, i.e. vertices
        self.X = kwargs["vertices"]

        self._pairwise_dist = None
        self.simplices      = []
        self.weights        = []

        ## Init info to construct all possible simplices
        self._get_info()
        ## construct VR complex
        self._construct_VR_complex()


    def _ncr(self, n, r):
        ''' See: https://stackoverflow.com/questions/4941753/is-there-a-math-ncr-function-in-python/4941846
        '''
        r     = min(r, n - r)
        numer = functools.reduce(operator.mul, range(n, n - r, -1), 1)
        denom = functools.reduce(operator.mul, range(1, r + 1), 1)
        return numer // denom


    def _get_info(self):
        '''
        '''
        ## can use math.comb(n, r) if Python3.8 but not everyone has it so use the custom func
        N = self.X.shape[0]
        self._num_possible_simplices = 0
        self._offset_idx = [] ## To get easier indexing later

        for r in range(1, self.max_dim + 1):
            cnt = self._ncr(N, r)
            self._num_possible_simplices += cnt
            self._offset_idx.append(self._num_possible_simplices)

        print("Total Number of possible simplices: {}".format(self._num_possible_simplices))


    def _compute_pairwise_dist(self):
        '''
        '''
        ## Update metric here
        self._pairwise_dist = scipy.spatial.distance.pdist(
            self.X, metric = "euclidean"
        )
        print("Max distance b/w any two points is: {}".format(self._pairwise_dist.max()))


    def _ret_pairwise_dist(self, idx):
        '''
        '''
        m = self.X.shape[0]
        ## see : https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.pdist.html
        return self._pairwise_dist[m * idx[:,0] + idx[:,1] - ((idx[:,0] + 2) * (idx[:,0] + 1)) // 2]


    def _add_vertices(self):
        ''' Add vertices labeled 0 to num_points - 1
        '''
        self._vertices = np.arange(self.X.shape[0], dtype = np.uint16)



    def _generate_combinations(self):
        '''
        '''
        for r in range(1, self.max_dim + 1):
            for simplex in itertools.combinations(self._vertices, r):
                yield simplex


    def _comb_index(self, n, r):
        ''' Return the indices of all possble r-combinations from a array of size n
        '''
        count = scipy.special.comb(n, r, exact = True)
        index = np.fromiter(itertools.chain.from_iterable(
            itertools.combinations(range(n), r)), np.uint16, count = count * r)
        return index.reshape(-1, r)


    def _maximal_vr(self):
        ''' Create a maximal Vietoris-Rips Complex
        '''
        ## Construct max simplices using the idea from Zomorodian 2010
        ## The maximal clique here is just all the vertices
        ## Create an array from 0 to npoints -  1
        self._add_vertices()
        ## Take all possible combinations of all possible dimensions
        for r in range(1, self.max_dim + 1):
            print("Adding {}-simplices".format(r))
            indices = self._comb_index(self._vertices.size, r)
            self.simplices.append(indices)
            print("Num {}-simplices added: {}".format(r, len(indices)))


    def _compute_weights(self):
        ''' Assign weights to all simplices in the Complex
        '''
        for n, nsimplices in enumerate(self.simplices):
            if n == 0:
                self.weights.append(np.zeros(len(nsimplices), dtype = float))
            elif n == 1:
                self.weights.append(self._ret_pairwise_dist(nsimplices))
            ## weight for any > 1-simplex is the max of all its edges
            else:
                indices = self._comb_index(n + 1, 2)
                self.weights.append(
                    np.amax(
                        self._ret_pairwise_dist(nsimplices[:, indices].reshape(
                            -1, 2)).reshape(len(nsimplices), self._ncr(n + 1, 2)), axis = 1))


    def _sort_simplices(self):
        ''' Sort simplices by weight
        '''
        for i in range(1, len(self.simplices)):
            sorted_idx = self.weights[i].argsort()
            self.weights[i]   = self.weights[i][sorted_idx]
            self.simplices[i] = self.simplices[i][sorted_idx]
        self.weights = np.hstack(self.weights)


    def _construct_VR_complex(self):
        ''' Construct a VR complex with upto n-simplices
        '''
        ## Compute pairwise distance between vertices/points
        self._compute_pairwise_dist()
        print("Computed Pairwise distances")
        print("--------------------------------------------------")
        self._maximal_vr()
        print("Done enumerating simplices.")
        print("--------------------------------------------------")
        self._compute_weights()
        print("Done computing weights for simplices.")
        print("--------------------------------------------------")
        self._sort_simplices()
        print("Done sorting simplices.")
        print("--------------------------------------------------")
        self._filtration_boundary_matrix()
        print("Done creating filtration boundary matrix")
        print("--------------------------------------------------")


    def _issubset(self, a, b):
        '''Return whether sequence `a` is a subset of sequence `b`'''
        return len(np.setdiff1d(a, b)) == 0



    def _apply_distance_cutff(self, row_indices, col_indices):
        ''' Apply distance cutoff to reduce number of simplices for boundary matrix
        reduction to be tractable
        '''


    def _create_simplices_for_lookup(self):
        '''
        '''
        i = 0
        self.simplices_for_lookup = []
        while i < self.weights.size:
            for simp in self.simplices:
                for s in simp:
                    self.simplices_for_lookup.append(s)
                    i += 1


    def _filtration_boundary_matrix_with_cutoff(self):
        '''
        '''
        idx_to_keep = self.weights <= self.epsilon
        self.weights = self.weights[idx_to_keep]
        self._create_simplices_for_lookup()
        self.simplices_for_lookup = list(itertools.compress(self.simplices_for_lookup, idx_to_keep))

        l = 1
        offsets = []
        for i, s in enumerate(self.simplices_for_lookup):
            if len(s) != l:
                offsets.append(i)
                l = len(s)

        offsets.append(len(self.simplices_for_lookup))
        print(offsets)

        row_indices = []
        col_indices = []
        for idx1, s1 in enumerate(self.simplices_for_lookup):
            for idx2 in range(offsets[len(s1)], offsets[-1]):
                s2 = self.simplices_for_lookup[idx2]
                if len(s1) >= len(s2): continue
                if sortednp.issubset(s1, s2):
                    row_indices.append(idx1)
                    col_indices.append(idx2)

        self.boundary_matrix = scipy.sparse.csc_matrix(
            (np.ones_like(row_indices, dtype = np.uint8), (row_indices, col_indices)),
            shape = (self.weights.size, self.weights.size), dtype = np.uint8)


    def _filtration_boundary_matrix(self):
        ''' Construct a filtration boundary matrix based on the ordered simplices
        '''
        if self.epsilon is not None:
            self._filtration_boundary_matrix_with_cutoff()
            return

        ## Get row and column indices for entries with a 1.
        row_indices       = []
        col_indices       = []
        num_conn_1simplex = []
        cum_sum           = []
        total_num_conn    = 0
        for i in range(1, self.max_dim):
            temp = self._ncr(self.X.shape[0] - 1, i)
            num_conn_1simplex.append(temp)
            cum_sum.append(total_num_conn)
            total_num_conn += temp

        for idx, simp1 in enumerate(self.simplices[0:-1]):
            for i, s1 in enumerate(simp1):
                if len(s1) == 1:
                    indices = np.empty(total_num_conn, dtype = np.uint32)
                    n = self.X.shape[0]
                    prev_idx = 0
                    for j, simp2 in enumerate(self.simplices[1:]):
                        indices[prev_idx: prev_idx + num_conn_1simplex[j]] = np.where(simp2 == s1)[0] + self._offset_idx[j]
                        prev_idx = prev_idx + num_conn_1simplex[j]
                        n -= 1
                    col_indices.append(indices)
                else:
                    ## Hard coding this for upto 4 simplices - need to find a way to make the elements inside this function call a callable
                    if len(s1) == 2:
                        foo = sortednp.kway_intersect(
                            col_indices[s1[0]][cum_sum[len(s1) - 1]:],
                            col_indices[s1[1]][cum_sum[len(s1) - 1]:])
                    elif len(s1) == 3:
                        foo = sortednp.kway_intersect(
                            col_indices[s1[0]][cum_sum[len(s1) - 1]:],
                            col_indices[s1[1]][cum_sum[len(s1) - 1]:],
                            col_indices[s1[2]][cum_sum[len(s1) - 1]:])
                    elif len(s1) == 4:
                        foo = sortednp.kway_intersect(
                            col_indices[s1[0]][cum_sum[len(s1) - 1]:],
                            col_indices[s1[1]][cum_sum[len(s1) - 1]:],
                            col_indices[s1[2]][cum_sum[len(s1) - 1]:],
                            col_indices[s1[3]][cum_sum[len(s1) - 1]:])
                    assert len(foo) != 0
                    col_indices.append(foo)


        row_indices = np.concatenate(
            [np.full(len(col_indices[i]),  i, dtype = np.uint32) \
             for i in range(0, self._offset_idx[-2])])
        col_indices = np.hstack(col_indices)

        self.boundary_matrix = scipy.sparse.csc_matrix(
            (np.ones_like(row_indices, dtype = np.uint8), (row_indices, col_indices)),
            shape = (self._num_possible_simplices, self._num_possible_simplices), dtype = np.uint8)

        # For large data sets think of a better way - maybe use some interval tree datastructure
        self._create_simplices_for_lookup()




def test():
    '''
    '''

    ## Sanity check to Edelsbrunner 2014 A Short Course .. Filtration Boundary Matrix 13.1
    a = np.array(
        [[0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0],
         [0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0],
         [0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
         [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0],
         [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0],
         [0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0],
         [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0],
         [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0],
         [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0],
         [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
         [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0],
         [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
         [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
         [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
         [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]])

    mat = scipy.sparse.csc_matrix(a)
    nrows, ncols = mat.shape
    indptr = mat.indptr
    indices = mat.indices
    print(indptr, indices)
    print(numerics.reduc(nrows, ncols, indptr, indices))


    ## A small test with data simulated from a circle -
    ## I got it by reading this blogpost : http://outlace.com/TDApart1.html
    # number of points to generate
    n = 15
    # generate space of parameter
    theta = np.linspace(0, 2.0 * np.pi, n)
    a, b, r = 0.0, 0.0, 5.0
    x = a + r * np.cos(theta)
    y = b +  r* np.sin(theta)
    x2 = np.random.uniform(-0.75,0.75, n) + x #add some "jitteriness" to the points
    y2 = np.random.uniform(-0.75,0.75, n) + y

    ## Do persistent homology
    vertices = np.array(list(zip(x2, y2)))
    vr_complex = VRComplex(vertices = vertices, n = 4)
    print(len(vr_complex.simplices))
    nrows, ncols = vr_complex.boundary_matrix.shape
    indptr = vr_complex.boundary_matrix.indptr
    indices = vr_complex.boundary_matrix.indices

    pivots, pivots_idx = numerics.reduc(nrows, ncols, indptr, indices)

    plotname = "circle_birth.death.png"
    intervals, homol_ranks, homol = persistence(
        pivots, pivots_idx, vr_complex.simplices_for_lookup,  vr_complex.weights)
    birth_death_plot(intervals, plotname)


if __name__ == '__main__':
    test()




