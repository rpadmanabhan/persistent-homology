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


    def tearDown(self):
        '''
        '''
        pass

class TestExampleA(unittest.TestCase):
    def setUp(self):
        '''
        '''
        self.vertices = set((abst_simplcl_cmplx.Vertex(label = "Cow"),
                             abst_simplcl_cmplx.Vertex(label = "Rabbit"),
                             abst_simplcl_cmplx.Vertex(label = "Horse"),
                             abst_simplcl_cmplx.Vertex(label = "Dog"),
                             abst_simplcl_cmplx.Vertex(label = "Fish"),
                             abst_simplcl_cmplx.Vertex(label = "Dolphin"),
                             abst_simplcl_cmplx.Vertex(label = "Oyster"),
                             abst_simplcl_cmplx.Vertex(label = "Broccoli"),
                             abst_simplcl_cmplx.Vertex(label = "Fern"),
                             abst_simplcl_cmplx.Vertex(label = "Onion"),
                             abst_simplcl_cmplx.Vertex(label = "Apple")))
        self.simplicial_complex = abst_simplcl_cmplx.ASC(vertices = self.vertices)
        self.simplicial_complex.add_connections(("Cow", "Rabbit"))
        self.simplicial_complex.add_connections(("Cow", "Horse"))
        self.simplicial_complex.add_connections(("Rabbit", "Horse"))
        self.simplicial_complex.add_connections(("Rabbit", "Dog"))
        self.simplicial_complex.add_connections(("Horse", "Dog"))
        self.simplicial_complex.add_connections(("Fish", "Dolphin"))
        self.simplicial_complex.add_connections(("Fish", "Oyster"))
        self.simplicial_complex.add_connections(("Dolphin", "Oyster"))
        self.simplicial_complex.add_connections(("Broccoli", "Fern"))
        self.simplicial_complex.add_connections(("Broccoli", "Onion"))
        self.simplicial_complex.add_connections(("Broccoli", "Apple"))
        self.simplicial_complex.add_connections(("Fern", "Onion"))
        self.simplicial_complex.add_connections(("Fern", "Apple"))
        self.simplicial_complex.add_connections(("Onion", "Apple"))
        self.simplicial_complex.add_connections(("Cow", "Rabbit", "Horse"))
        self.simplicial_complex.add_connections(("Cow", "Rabbit", "Dog"))
        self.simplicial_complex.add_connections(("Cow", "Horse", "Dog"))
        self.simplicial_complex.add_connections(("Rabbit", "Horse", "Dog"))
        self.simplicial_complex.add_connections(("Fish", "Dolphin", "Oyster"))
        self.simplicial_complex.add_connections(("Broccoli", "Fern", "Onion"))
        self.simplicial_complex.add_connections(("Broccoli", "Fern", "Apple"))
        self.simplicial_complex.add_connections(("Broccoli", "Onion", "Apple"))
        self.simplicial_complex.add_connections(("Fern", "Onion", "Apple"))


    def test_exampleA(self):
        '''
        '''
        # TO DO: Add better checks here, pull out the harded coded stuff above to something more generic and add example B
        self.assertEqual(len(self.simplicial_complex.ret_all_simplices(2)), 9)
        self.assertEqual(len(self.simplicial_complex.ret_all_simplices(1)), 14)
        self.assertEqual(len(self.simplicial_complex.ret_all_simplices(0)), len(self.vertices))



if __name__ == '__main__':
    unittest.main()
