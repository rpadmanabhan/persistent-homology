## stdlib
import os
import sys
import unittest

## 3rd party
import numpy as np

## our modules
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import matrix_ops


class TestRowEchelon(unittest.TestCase):
    '''
    '''
    def test_row_echelon_Z2(self):
        ''' Checks for correct computation of row echelon form in Z2
        '''
        A = np.array(
            [[1, 1, 0, 0], [1, 0, 1, 0], [0, 1, 1, 0],
             [1, 0, 0, 1], [0, 1, 0, 1], [0, 0, 1, 1]],
            dtype = np.uint8)
        rank = matrix_ops.row_echelon_Z2(A)
        self.assertEqual(rank, 3)
        self.assertEqual(
            (A ==
             np.array([[1, 1, 0, 0], [0, 1, 1, 0], [0, 0, 1, 1],
                       [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]],
                      dtype = np.uint8)).all(), True)


class TestReducedRowEchelon(unittest.TestCase):
    '''
    '''
    def test_reduced_row_echelon_Z2(self):
        ''' Checks for correct computation of reduced row echelon form in Z2
        '''
        A = np.array(
            [[1, 1, 0, 0], [1, 0, 1, 0], [0, 1, 1, 0],
             [1, 0, 0, 1], [0, 1, 0, 1], [0, 0, 1, 1]],
            dtype = np.uint8)
        rank = matrix_ops.reduced_row_echelon_Z2(A)
        self.assertEqual(rank, 3)
        self.assertEqual(
            (A ==
             np.array([[1, 0, 0, 1], [0, 1, 0, 1], [0, 0, 1, 1],
                       [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]],
                      dtype = np.uint8)).all(), True)


class TestInverse(unittest.TestCase):
    '''
    '''
    def test_inverse_Z2(self):
        ''' Checks for correct computation of matrix inverse in Z2
        '''
        A = np.array(
            [[1, 0, 0, 0], [1, 1, 0, 0], [1, 1, 1, 0], [1, 1, 1, 1]],
            dtype = np.uint8)
        self.assertEqual(
            (matrix_ops.inverse_Z2(A) ==
             np.array([[1, 0, 0, 0], [1, 1, 0, 0], [0, 1, 1, 0], [0, 0, 1, 1]],
                      dtype = np.uint8)).all(), True)


class TestSmithNormalForm(unittest.TestCase):
    '''
    '''
    def test_smith_normal_form_Z2(self):
        ''' Checks for correct computation of Smith Normal From in Z2
        '''
        A = np.array(
            [[1, 0, 0, 1], [1, 1, 0, 0], [0, 1, 1, 0], [0, 0, 1, 1]],
            dtype = np.uint8)
        rank_A, P_inv, D, Q = matrix_ops.smith_normal_form(A)
        self.assertEqual(rank_A, 3)
        self.assertEqual(
            (D ==
             np.array([[1, 0, 0, 0], [0, 1, 0, 0],
                       [0, 0, 1, 0], [0, 0, 0, 0]])).all(), True)
        self.assertEqual(
            (P_inv ==
             np.array([[1, 0, 0, 0], [1, 1, 0, 0],
                       [0, 1, 1, 0], [0, 0, 1, 1]])).all(), True)
        self.assertEqual(
            (Q ==
             np.array([[1, 0, 0, 1], [0, 1, 0, 1],
                       [0, 0, 1, 1], [0, 0, 0, 1]])).all(), True)


class TestBasisImKer(unittest.TestCase):
    '''
    '''
    def test_basis_img_ker_Z2(self):
        ''' Checks for correct computation of the basis for images and kernels
        '''
        A = np.array(
            [[1, 0, 0, 1], [1, 1, 0, 0], [0, 1, 1, 0], [0, 0, 1, 1]],
            dtype = np.uint8)
        im, ker = matrix_ops.basis_img_ker_Z2(A)
        self.assertEqual(
            (im == np.array([[1, 0, 0], [1, 1, 0],
                             [0, 1, 1], [0, 0, 1]], dtype = np.uint8)).all(), True)
        self.assertEqual(
            (ker == np.array([[1], [1], [1], [1]], dtype = np.uint8)).all(), True)


class TestHomology(unittest.TestCase):
    '''
    '''
    def test_homology_Z2(self):
        ''' Checks for correct computation of the basis for homology
        '''
        d_0 = np.array([[0, 0, 0, 0]], dtype = np.uint8)
        d_1 = np.array(
            [[1, 0, 0, 1], [1, 1, 0, 0], [0, 1, 1, 0], [0, 0, 1, 1]],
            dtype = np.uint8)

        basis_H0 = matrix_ops.basis_homology_Z2(d_0, d_1)
        self.assertEqual(
            (basis_H0 == np.array([[1], [0], [0], [0]], dtype = np.uint8)).all(), True)

        d_2 = np.array([[0]], dtype = np.uint8)
        basis_H1 = matrix_ops.basis_homology_Z2(d_1, d_2)
        self.assertEqual(
            (basis_H1 == np.array([[1], [1], [1], [1]], dtype = np.uint8)).all(), True)


if __name__ == '__main__':
    unittest.main()
