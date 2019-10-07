from neopz import *
vec = TPZVecDouble(2, 3)
vec.Resize(2)
vec_item = vec[1]
print(vec_item)
mat = TPZMatrix(2, 2, 3.)
mat_item = mat.GetVal(1, 1)
print(mat_item)
print(TPZQuadrilateral.NSideNodes(7))
