from neopz import *
vec = PZVecDouble(2, 3)
vec.Resize(2)
vec_item = vec[1]
print(vec_item)
mat = PZMatrix(2, 2, 3.)
mat_item = mat.GetVal(1, 1)
print(mat_item)
print(PZQuad.SideDimension(9))
