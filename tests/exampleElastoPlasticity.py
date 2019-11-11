from neopz import *
import math
gmesh = TPZGeoMesh()
read = TPZGmshReader()
gmesh = read.GeometricGmshMesh4("tests/geometric-mesh/wellbore.msh", gmesh)
gmesh.BuildConnectivity()
mat = TPZMatElastoPlastic2DMC(1,2)
cmesh = TPZCompMesh(gmesh)
cmesh.SetDefaultOrder(1)
dim = gmesh.Dimension()

# Material properties
cohesion = 5
phi = 20*math.pi/180
E = 2000
nu = 0.2

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

an = TPZAnalysis(cmesh,1)
struc_mat = TPZSymetricSpStructMatrix(cmesh)
# struc_mat = TPZParFrontStructMatrix(cmesh)
# struc_mat = TPZSkylineStructMatrix(cmesh)
# struc_mat = TPZSkylineNSymStructMatrix(cmesh)

an.SetStructuralMatrix(struc_mat)
stepsol = TPZStepSolver()

stepsol.SetDirect(ELDLt)
# stepsol.SetDirect(ECholesky)
# stepsol.SetDirect(ELU)


an.SetSolver(stepsol)

sol = an.Solution()
sol.Zero()
an.LoadSolution(sol)
du  = an.Solution()

an.Assemble()

nit = 2
tol = 1e-7


for it in range(nit):
	print("iteration # ",it)
	an.Solve()
	dmsh = an.Mesh()
	ddu  = dmsh.Solution()
	du   = ddu+du
	norm_du = du.Norm()
	print("norm_du",norm_du)

	an.LoadSolution(du)
	an.Assemble()

	norm_res = an.NormRhs()
	stop_criterion = norm_res < tol
	print("Nonlinear process: residue norm =  ", norm_res)

	if stop_criterion:
		an.AcceptPseudoTimeStepSolution();
		break
	else:
		print("# The process needs more number of iterations? #")


# name = str("elastoplasticWellbore.vtk")
# scalnames = TPZVecString(1)
# vecnames = TPZVecString(1)
# tensnames = TPZVecString(1)
# scalnames[0]="FailureType"
# vecnames[0] ="Displacement"
# tensnames[0]="Stress"

# an.DefineGraphMesh(2,scalnames,vecnames, tensnames, name)
# an.PostProcess(0,2)
print("# End of Simulation! #")








