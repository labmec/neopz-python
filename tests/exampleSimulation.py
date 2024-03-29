from neopz import*
gmesh = TPZGeoMesh()
read = TPZGmshReader()
gmesh = read.GeometricGmshMesh4("tests/geometric-mesh/simple_2D_coarse.msh", gmesh)
gmesh.BuildConnectivity()
mat = TPZMatPoisson3d(1,2)
cmesh = TPZCompMesh(gmesh)
val = cmesh.InsertMaterialObject(mat)

D_Type = 0
N_Type = 1
val1 = TPZFMatrix(1,1,0.0)
val2 = TPZFMatrix(1,1,0.0)
bc_1 = mat.CreateBC(mat, 2,N_Type,val1,val2)
cmesh.InsertMaterialObject(bc_1)

val3 = TPZFMatrix(1,1,30.0)
bc_2 = mat.CreateBC(mat, 4,D_Type,val1,val3)
cmesh.InsertMaterialObject(bc_2)

val4 = TPZFMatrix(1,1,10.0)
bc_3 = mat.CreateBC(mat, 5,D_Type,val1,val4)
cmesh.InsertMaterialObject(bc_3)

cmesh.SetDimModel(2)

cmesh.SetAllCreateFunctionsContinuous()
cmesh.SetDefaultOrder(1)
cmesh.AutoBuild()

an = TPZAnalysis(cmesh, 1)
struc_mat = TPZSymetricSpStructMatrix(cmesh)
an.SetStructuralMatrix(struc_mat)
stepsol = TPZStepSolver()
stepsol.SetDirect(ECholesky)
an.SetSolver(stepsol)
an.Assemble()
an.Solve()

name = str("resultado.vtk")
scalnames = TPZVec_string(1)
vecnames = TPZVec_string(1)
scalnames[0]="state"
vecnames[0]="Flux"

an.DefineGraphMesh(2,scalnames,vecnames, name)
an.PostProcess(0,2)








