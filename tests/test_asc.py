## stdlib
import os
import sys
import unittest

## our modules
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import abst_simplcl_cmplx



class TestASCInit(unittest.TestCase):

    def setUp(self):
        '''
        '''
        self.vertices = set((abst_simplcl_cmplx.Vertex(label = "Cow"),
                            abst_simplcl_cmplx.Vertex(label = "Rabbit"),
                            abst_simplcl_cmplx.Vertex(label = "Horse"),
                            abst_simplcl_cmplx.Vertex(label = "Dog")))

    def test_init_ASC(self):
        ''' Test if vertices are correctly initialized
        '''
        simplicial_complex = abst_simplcl_cmplx.ASC(vertices = self.vertices)
        self.assertEqual(
            simplicial_complex.root.connections,
            self.vertices
        )

    def tearDown(self):
        '''
        '''
        pass


class TestASCInsertConnection(unittest.TestCase):

    def setUp(self):
        '''
        '''
        self.vertices = set((abst_simplcl_cmplx.Vertex(label = "Cow"),
                            abst_simplcl_cmplx.Vertex(label = "Rabbit"),
                            abst_simplcl_cmplx.Vertex(label = "Horse"),
                            abst_simplcl_cmplx.Vertex(label = "Dog")))
        self.simplicial_complex = abst_simplcl_cmplx.ASC(vertices = self.vertices)


    def test_insert_connections(self):
        ''' Insert some connections and test if correct structure is maintained
        '''
        self.simplicial_complex.add_connections(("Cow", "Rabbit"))
        self.simplicial_complex.add_connections(("Cow", "Horse"))
        self.simplicial_complex.add_connections(("Rabbit", "Horse"))
        self.simplicial_complex.add_connections(("Cow", "Rabbit", "Horse"))

        self.assertEqual(self.simplicial_complex.vertex_tracker[("Rabbit", "Horse", 2)].level, 2)
        self.assertEqual(self.simplicial_complex.vertex_tracker[("Cow", "Horse", 1)].level, 1)
        self.assertEqual(self.simplicial_complex.vertex_tracker[("Rabbit", "Horse", 1)].level, 1)


        self.assertEqual(self.simplicial_complex.ret_all_simplices(2), [["Cow", "Rabbit", "Horse"]])
        self.assertEqual(sorted(self.simplicial_complex.ret_all_simplices(1)), sorted([["Cow", "Rabbit"], ["Cow", "Horse"], ["Rabbit", "Horse"]]))
        self.assertEqual(sorted(self.simplicial_complex.ret_all_simplices(0)), sorted([["Cow"], ["Rabbit"], ["Horse"], ["Dog"]]))
        print(self.simplicial_complex)

    def tearDown(self):
        '''
        '''
        pass


if __name__ == '__main__':
    unittest.main()
