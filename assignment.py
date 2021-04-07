## stdlib
import itertools
import sys

## our modules
import abst_simplcl_cmplx
import matrix_ops

## 3rd party
import numpy as np
import scipy
import scipy.linalg


# Team : Bill Lee, Raghavendra Padmanabhan and Francisco Vargas
#--------------------------------------------------------------------------------------


class DataCoding1:
    ''' Simplicial Complex for the data in Coding assignment 1
    '''
    def __init__(self, *args, **kwargs):
        '''
        '''
        self.vertex_labels      = kwargs["vertex_labels"]
        self.vertex_connections = kwargs["vertex_connections"]
        ## Build the simplicial complex
        self.vertices = set(
            (abst_simplcl_cmplx.Vertex(label = label) for label in self.vertex_labels))
        self.simplicial_complex = abst_simplcl_cmplx.ASC(vertices = self.vertices)
        for connections in self.vertex_connections:
            self.simplicial_complex.add_connections(connections)

## Part A
partA = DataCoding1(
    vertex_labels = set(
        {"Cow", "Rabbit", "Horse", "Dog", "Fish",
         "Dolphin", "Oyster", "Broccoli", "Fern", "Onion", "Apple"}),
    vertex_connections = [
        ("Cow", "Rabbit"), ("Cow", "Horse"), ("Cow", "Dog"),
        ("Rabbit", "Horse"), ("Rabbit", "Dog"), ("Horse", "Dog"),
        ("Fish", "Dolphin"), ("Fish", "Oyster"), ("Dolphin", "Oyster"),
        ("Broccoli", "Fern"), ("Broccoli", "Onion"), ("Broccoli", "Apple"),
        ("Fern", "Onion"), ("Fern", "Apple"),
        ("Onion", "Apple"),
        ("Cow", "Rabbit", "Horse"), ("Cow", "Rabbit", "Dog"),
        ("Cow", "Horse", "Dog"), ("Rabbit", "Horse", "Dog"),
        ("Fish", "Dolphin", "Oyster"),
        ("Broccoli", "Fern", "Onion"),
        ("Broccoli", "Fern", "Apple"),
        ("Broccoli", "Onion", "Apple"),
        ("Fern", "Onion", "Apple")])

## Part B
partB = DataCoding1(
    vertex_labels = set(
        {"Cow", "Rabbit", "Horse", "Dog", "Fish",
         "Dolphin", "Oyster", "Broccoli", "Fern", "Onion", "Apple"}),
    vertex_connections = [
        ("Cow", "Rabbit"), ("Cow", "Fish"), ("Cow", "Oyster"),
        ("Cow", "Broccoli"), ("Cow", "Onion"), ("Cow", "Apple"),
        ("Rabbit", "Fish"), ("Rabbit", "Oyster"), ("Rabbit", "Broccoli"),
        ("Rabbit", "Onion"), ("Rabbit", "Apple"), ("Fish", "Oyster"),
        ("Fish", "Broccoli"), ("Fish", "Onion"), ("Fish", "Apple"),
        ("Oyster", "Broccoli"), ("Oyster", "Onion"), ("Oyster", "Apple"),
        ("Broccoli", "Onion"), ("Broccoli", "Apple"), ("Onion", "Apple"),
        ("Horse", "Dog"), ("Horse", "Dolphin"), ("Horse", "Fern"),
        ("Dog", "Dolphin"), ("Dog", "Fern"), ("Dolphin", "Fern"),
        ("Cow", "Broccoli", "Apple"), ("Cow", "Onion", "Apple"),
        ("Rabbit", "Broccoli", "Apple"), ("Rabbit", "Onion", "Apple"),
        ("Fish", "Broccoli", "Apple"), ("Fish", "Onion", "Apple"),
        ("Oyster", "Broccoli", "Apple"), ("Oyster", "Onion", "Apple")])



## Took this code as is from the top answer : https://stackoverflow.com/questions/17129290/numpy-2d-and-1d-array-to-latex-bmatrix
def bmatrix(a):
    """Returns a LaTeX bmatrix

    :a: numpy array
    :returns: LaTeX bmatrix as a string
    """
    if len(a.shape) > 2:
        raise ValueError('bmatrix can at most display two dimensions')
    lines = str(a).replace('[', '').replace(']', '').splitlines()
    rv = [r'\begin{bmatrix}']
    rv += ['  ' + ' & '.join(l.split()) + r'\\' for l in lines]
    rv +=  [r'\end{bmatrix}_{' + '{} \\times {}'.format(a.shape[0], a.shape[1]) + '}']
    return '\n'.join(rv)


def matrix_with_labels(a, row_labels, col_labels, matrix_name):
    ''' Returns a LaTeX blockarray matrix with row and column headers
    '''
    lines = str(a).replace('[', '').replace(']', '').splitlines()
    latex_output = [r'$\begin{small}']
    latex_output += [r'\mathbf{' + matrix_name + r'}' + r' = \begin{blockarray}{r*{' + str(len(col_labels)) + r'}{c}}']
    col_labels = tuple(r'\rotatebox{-90}{' + str(e).replace("'", "") + r'}' for e in col_labels)
    latex_output += [' & ' + ' & '.join(col_labels) + r'\\']
    latex_output += [r'\begin{block}{ r!{\,}(' + 'c'*len(col_labels) + r')}']
    for i, l in enumerate(lines):
        temp = l.split()
        row_label = str(row_labels[i]).replace("'", "")
        temp.insert(0, row_label)
        latex_output.append('  ' + ' & '.join(temp) + r'\\')
    latex_output += [r'\end{block}']
    latex_output += [r'\end{blockarray}_{' + str(len(row_labels)) + r' \times ' + str(len(col_labels)) + r'}']
    latex_output += [r'\end{small}$']

    return '\n'.join(latex_output)


def submission2(data, label):
    ''' Submission for coding assignment 2 - Boundaries of the p-chains.
    '''
    ## Boundaries in C_2 are the boundaries of all 2-simplices, i.e. all edges or 1-simplices
    print("Boundaries in $C_2$ are: \n")
    for simplex in data.simplicial_complex.ret_all_simplices(2):
        out = []
        for b in itertools.combinations(list(simplex), 2):
            out.append(str(tuple(b)))
        out = " + ".join(out)
        print("Boundary of the 2-simplex: {} is: {}\n".format(tuple(simplex), out))
    ## Boundaries in C_1 are the boundaries of all 1-simplices, i.e. all vertices or  0-simplices
    print("\nBoundaries in $C_1$ are: \n")
    for simplex in data.simplicial_complex.ret_all_simplices(1):
        v1, v2 = list(simplex)
        print("Boundary of the 1-simplex: {} is: {} + {}\n".format(tuple(simplex), v1, v2))
    ## Boundaries of C_0 are the boundaries of all 0-simplices, i.e. just 0
    print("\nBoundaries in $C_0$ are: \n")
    print("\nBoundary of a vertex is 0, so all boundaries in $C_0$ are just 0.")


def submission3(data, label):
    ''' Submission for coding assignment 3 - Boundary Matrices
    '''
    ## Using numpy as it prints larger matrices better
    print("Boundary Matrices for {}:\n".format(label))
    for d in [2, 1, 0]:
        ## Get the boundary matrix and print it
        matrix, c_n_1_generators, c_n_generators = data.simplicial_complex.ret_boundary_matrix(d)
        ## store the generators as a tuple(dim > 1) or str(dim==1) for cleaner printing
        c_n_1_generators = [tuple(sorted(e)) if len(e) > 1 else str(next(iter(e))) for e in c_n_1_generators]
        c_n_generators   = [tuple(sorted(e)) if len(e) > 1 else str(next(iter(e))) for e in c_n_generators]

        matrix_name = "\\partial_{}".format(d)
        print(matrix_with_labels(matrix, c_n_1_generators, c_n_generators, matrix_name))
        print("\n")


        c_n_vecs = {}
        ## Get all generators in vector form
        for idx, generator in enumerate(c_n_generators):
            c_n_vecs[generator] = np.zeros((len(c_n_generators), 1), dtype = np.bool)
            c_n_vecs[generator][idx] = 1
        ## Boundaries for some linear combinations of our generators
        print("Some select examples of boundaries computed using this matrix :\n")
        for k in range(1, 4):
            for linear_combination in itertools.combinations(c_n_generators, k):
                vector = np.zeros((len(c_n_generators), 1), dtype = np.bool)
                for generator in linear_combination:
                    vector = vector ^ c_n_vecs[generator]

                boundary = (matrix @ vector)%2
                boundary_text = " + ".join(
                    [str(e) for e, mask in zip(c_n_1_generators, boundary) if mask])
                if boundary_text == "":
                    boundary_text = "0"
                linear_combination_text = " + ".join((str(e).replace("'", "") for e in linear_combination))
                print("Boundary of {} is: {}".format(linear_combination_text, boundary_text.replace("'", "")))
                print("\n")
                ## There are too many combinations so just print 1 for each k > 1
                if k > 1:
                    break
        print("\n")


def submission4(data, label):
    ''' Submission for coding assignment 4 - Kernels, Homologies, Ranks
    '''
    rank = {"ker": [0, 0, 0], "img": [0, 0, 0], "homol": [0, 0, 0]}
    print("\n{}\n".format(label))
    print("\nKernel of the boundary matrices:\n")
    for d in [2, 1, 0]:
        print("\nC_{}:\n".format(d))
        matrix, c_n_1_generators, c_n_generators = data.simplicial_complex.ret_boundary_matrix(d)
        c_n_generators   = [tuple(sorted(e)) if len(e) > 1 else str(next(iter(e))) for e in c_n_generators]
        basis_img, basis_ker = matrix_ops.basis_img_ker_Z2(matrix)
        rank["ker"][d] = np.shape(basis_ker)[1]
        rank["img"][d] = np.shape(basis_img)[1]
        for ker_col in basis_ker.T:
            ker_text = " + ".join(
                [str(e) for e, mask in zip(c_n_generators, ker_col) if mask]).replace("'", "")
            if ker_text == "":
                ker_text = "0"
            print(ker_text)

    print("\nHomologies:\n")
    for d in [2, 1, 0]:
        print("\nC_{}:\n".format(d))
        matrix_d, c_n_1_generators, c_n_generators = data.simplicial_complex.ret_boundary_matrix(d)
        c_n_generators   = [tuple(sorted(e)) if len(e) > 1 else str(next(iter(e))) for e in c_n_generators]
        matrix_d1, _, _ = data.simplicial_complex.ret_boundary_matrix(d + 1)
        basis_homol = matrix_ops.basis_homology_Z2(matrix_d, matrix_d1)
        rank["homol"][d] = np.shape(basis_homol)[1]
        for homol_col in basis_homol.T:
            homol_text = " + ".join(
                [str(e) for e, mask in zip(c_n_generators, homol_col) if mask]).replace("'", "")
            if homol_text == "":
                homol_text = "0"
            print(homol_text)

    print("\nRanks of Image, Kernel and Homology:\n")
    for d in [2, 1, 0]:
        print("\nC_{}:\n".format(d))
        print("Rank of the Image is: {}".format(rank["img"][d]))
        print("Rank of the Kernel is: {}".format(rank["ker"][d]))
        print("Rank of the Homology is: {}".format(rank["homol"][d]))
    print("\n-------------------------------------------------------------------------\n")

if __name__ == '__main__':
    if len(sys.argv) == 1:
        print("Specify param coding2 or coding3 for specific output")
        sys.exit(0)
    if sys.argv[1].lower() == "coding2":
        print("PartA\n")
        submission2(partA, "Part A")
        print("\n-------------------------------------------------------------------\n")
        print("PartB\n")
        submission2(partB, "Part B")
        print("\n-------------------------------------------------------------------\n")
    elif sys.argv[1].lower() == "coding3":
        submission3(partA, "Part A")
        submission3(partB, "Part B")
    elif sys.argv[1].lower() == "coding4":
        submission4(partA, "Part A")
        submission4(partB, "Part B")

