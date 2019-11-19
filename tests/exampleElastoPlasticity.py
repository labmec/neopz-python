from neopz import *
import math

# Functions to add geometry
gmesh = TPZGeoMesh()
read = TPZGmshReader()
gmesh = read.GeometricGmshMesh4("tests/geometric-mesh/wellbore.msh", gmesh)
gmesh.BuildConnectivity()

# Functions to add computational mesh
mat = TPZMatElastoPlastic2DMC(1,1)
cmesh = TPZCompMesh(gmesh)
cmesh.SetDefaultOrder(1)
dim = gmesh.Dimension()

# Material properties
cohesion = 5.0
phi = 20*math.pi/180
E = 2000
nu = 0.2

# Functions to add elasoplastic material
LER = TPZElasticResponse()
LER.SetEngineeringData(E,nu)
LEMC = TPZPlasticStepPVMC()
LEMC.SetElasticResponse(LER)
MC = TPZYCMohrCoulombPV()
MC = LEMC.YC()
MC.SetUp(phi, phi, cohesion, LER)
default_mem = TPZElastoPlasticMem()
default_mem.SetElasticResponse(LER)
stress = TPZTensor()
stress.Zero()
default_mem.SetStress(stress)
default_mem.SetPlasticState(LEMC.fN())
mat.SetPlasticityModel(LEMC)
mat.SetDefaultMem(default_mem)
cmesh.InsertMaterialObject(mat)

# Functions to add boundary conditions
# D_Type = 0
# N_Type = 1
val1 = TPZFMatrix(dim,dim,0.0)
val2 = TPZFMatrix(dim,dim,0.0)

val2.SetItem(0,0, -10)
bc_inner = TPZBndCondWithMem(mat, 2, 6, val1, val2)
cmesh.InsertMaterialObject(bc_inner)

val2.SetItem(0,0, -0.00625)
bc_outer = TPZBndCondWithMem(mat, 3, 6, val1, val2)
cmesh.InsertMaterialObject(bc_outer)

val2.SetItem(0,0, 1)
val2.SetItem(1,0, 0)
bc_ux_fixed = TPZBndCondWithMem(mat, 4, 3, val1, val2)
cmesh.InsertMaterialObject(bc_ux_fixed)

val2.SetItem(0,0, 0)
val2.SetItem(1,0, 1)
bc_uy_fixed = TPZBndCondWithMem(mat, 5, 3, val1, val2)
cmesh.InsertMaterialObject(bc_uy_fixed)

cmesh.SetDimModel(2)
cmesh.SetAllCreateFunctionsContinuousWithMem()
cmesh.ApproxSpace().CreateWithMemory(1)
cmesh.AutoBuild()

# Functions to do analysis
an = TPZAnalysis(cmesh,1)
# struc_mat = TPZSymetricSpStructMatrix(cmesh)
# struc_mat = TPZParFrontStructMatrix(cmesh)
# struc_mat = TPZSkylineStructMatrix(cmesh)
struc_mat = TPZSkylineNSymStructMatrix(cmesh)

an.SetStructuralMatrix(struc_mat)
stepsol = TPZStepSolver()

# stepsol.SetDirect(ELDLt)
# stepsol.SetDirect(ECholesky)
stepsol.SetDirect(ELU)

an.SetSolver(stepsol)
sol = an.Solution()
sol.Zero()
an.LoadSolution(sol)
du  = an.Solution()
an.Assemble()

nit = 5
tol = 1e-7

for it in range(nit):

	an.Assemble()
	print("iteration # ",it)

# New structure
	an.SetStructMatrixDecomposed(False)
	rn = an.Rhs()
	rn *= 1.0
	an.SetRhs(rn)
# New structure

	# an.PrintRhs()

	an.Solve()

	# an.PrintMatrix()

	dmsh = an.Mesh()
	ddu  = dmsh.Solution()
	du   = ddu+du
	norm_du = du.Norm()
	print("norm_du",norm_du)

	an.LoadSolution(du)
	# print(ddu)
	# print(du)
	an.Assemble()
	# print(ddu)
	# print(du)


	norm_res = an.NormRhs()
	stop_criterion = (norm_res < tol and norm_du < tol)
	print("Nonlinear process: delta_u norm =  ", norm_du)
	print("Nonlinear process: residue norm =  ", norm_res)

	if stop_criterion:
		an.AcceptPseudoTimeStepSolution();
		norm_res = an.NormRhs()
		print("Nonlinear process converged with residue norm =  ", norm_res)
		print("Number of iterations = ", it+1)
		break
	else:
		an.AcceptPseudoTimeStepSolution();
		print("# The process needs more number of iterations? #")

# Functions to do post process
post = TPZPostProcAnalysis()
post.SetCompMesh(cmesh,True)
post_mat_id = TPZVecInt(1,mat.Id())

name = str("elastoplasticWellbore.vtk")
scalnames = TPZVecString(1)
vecnames = TPZVecString(1)
tensnames = TPZVecString(3)
scalnames[0]="FailureType"
vecnames[0] ="Displacement"
tensnames[0]="Stress"
tensnames[1]="Strain"
tensnames[2]="StrainPlastic"

names = TPZVecString(5)
names[0] = scalnames[0]
names[1] = vecnames[0]
names[2] = tensnames[0]
names[3] = tensnames[1]
names[4] = tensnames[2]


post.SetPostProcessVariables(post_mat_id,names)
smatrix = TPZFStructMatrix(post.Mesh())
post.SetStructuralMatrix(smatrix)
post.TransferSolution()
post.DefineGraphMesh(2,scalnames,vecnames, tensnames, name)
post.PostProcess(0,2)

print("# End of Simulation! #")








