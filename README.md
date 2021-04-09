# instalation

run in the root dir
```
pip install -e .
```

# abstract-simplicial-complex

## Representation of an abstract simplicial complex as a tree datastructure in Python3.

#### Usage - Assuming python3.5+, should work for any python3
```python
from tda import abst_simplcl_cmplx

## Initialize an abstract simplicial complex
simplicial_complex = abst_simplcl_cmplx.ASC(
          vertices = [abst_simplcl_cmplx.Vertex(label = "Cow"), abst_simplcl_cmplx.Vertex(label = "Rabbit"),
                      abst_simplcl_cmplx.Vertex(label = "Rabbit")])
## Add some connections
simplicial_complex.add_connections(("Cow", "Rabbit"))
simplicial_complex.add_connections(("Cow", "Horse"))
simplicial_complex.add_connections(("Rabbit", "Horse"))
simplicial_complex.add_connections(("Cow", "Rabbit", "Horse"))

## Return back the 2-simplex just inserted above
simplicial_complex.ret_all_simplices(2)
```

#### Run tests
```
python3 -m unittest discover tests/ -v
```
