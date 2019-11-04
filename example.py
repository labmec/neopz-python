from neopz import *

print("Testando sides do triangulo")
for sideId in range(TPZQuadrilateral.NumSides()):
    print(TPZQuadrilateral.NSideNodes(sideId))

trans = TPZTriangle.SideToSideTransform(TPZTriangle.NumSides() - 1, 3)
print(trans)

print("Testando print de containers")

stack = TPZStackInt(3, 4)
print(stack)

matrix = TPZFMatrix(3, 1, 1.)
print(matrix)

vec = TPZVecDouble(3, 1.)
print(vec)
