from neopz import*
import math as math

# Exact solution for Poisson 2D
def F_source(x = TPZManVecDouble(), sol = TPZManVecDouble(), dsol = TPZFMatrix() ):
	dsol.Resize(3,3)
	sol.Resize(3)
	dsol.Zero()
#    sol.SetItem(0,1.0)



# MAIN
gmesh = TPZGeoMesh()
read = TPZGmshReader()
gmesh = read.GeometricGmshMesh4("tests/geometric-mesh/simple_2D.msh", gmesh)
gmesh.BuildConnectivity()
mat = TPZMatLaplacian(1,2)
mat.SetForcingFunctionExact(F_source)
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

#val5 = TPZFMatrix(2,1,0.0)
#val5.SetItem(0,0, 0)
#val5.SetItem(1,0, 1)
#bc_uy_fixed = mat.CreateBC(mat,6,D_Type, val1, val5)
#cmesh.InsertMaterialObject(bc_uy_fixed)

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
print(an.Solution())
# an.SetExact(Sol_exact)
#an.SetExact2()

# an.SetForcingFunction(F_source)

scalnames = TPZVecString(2)
vecnames = TPZVecString(1)
scalnames[0]="state"
scalnames[1]="ExactSolution"
# vecnames[0]="ExactFlux"
vecnames[0]="Flux"

# an.DefineGraphMesh(2,scalnames,vecnames, "resultado.vtk")
# # print("\n\nQuales tipos de maricas?")
# an.PostProcess(0,2)





# Functions to do post process
post = TPZPostProcAnalysis()
post.SetCompMesh(cmesh,False)
post_mat_id = TPZVecInt(1,mat.Id())

tensname = TPZVecString(1,"state")
names = TPZVecString(3)
names[0] = scalnames[0]
names[1] = vecnames[0]
names[2] = scalnames[1]


post.SetPostProcessVariables(post_mat_id,names)
smatrix = TPZFStructMatrix(post.Mesh())
post.SetStructuralMatrix(smatrix)
post.TransferSolution()
#print("\nIt's on a destructor!!!!\n\n")
post.DefineGraphMesh(2,scalnames,vecnames, "resultado.vtk")
post.PostProcess(0,2)

































# # Exact solution for Poisson 2D




# def F_source(x = TPZManVecDouble() , f = TPZManVecDouble(2)):
#     Pi = math.pi
    
#     f.resize(2)
    
#     xv = x[0]
#     yv = x[1]
#     # zv = x[2]

#     f_x = + 8.0*Pi*Pi*math.cos(2.0*Pi*yv)*math.sin(2.0*Pi*xv)
#     f_y = - 8.0*Pi*Pi*math.cos(2.0*Pi*xv)*math.sin(2.0*Pi*yv)

    
#     f[0] = f_x # x direction
#     f[1] = f_y # y direction
