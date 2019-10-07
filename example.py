from neopz import *

print("Testando sides do quadrilatero")
for sideId in range(TPZQuadrilateral.NumSides()):
    print(TPZQuadrilateral.NSideNodes(sideId))

print("Testando sides do cubo")
for sideId in range(TPZCube.NumSides()):
    print(TPZCube.NSideNodes(sideId))
