import siconos.algebra as alg
import numpy as np

# ------ Vector ----

def test_vector():
    
    # Create a vector
    size = 12
    
    v = alg.SiconosVector(size=size, storage_type=alg.UblasType.dense)
    v.print()
    assert v.size() == size
    
    # Numpy

    w = np.array(alg.SiconosVector(size=size, storage_type=alg.UblasType.dense), copy=False)
    w[...] = 3.
    assert np.sum(w) == 12*3.

    z = np.add(w,v)
    assert np.sum(z) == 12*3.


# --- Matrix ---

def test_matrix():
    
    mat = alg.SiconosMatrix(row=3, col=4, storage_type=alg.UblasType.dense)
    print(mat)

if __name__ == '__main__':
    test_vector()
    test_matrix()
