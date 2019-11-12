from neopz import*


gmesh = TPZGeoMesh()
read = TPZGmshReader()
gmesh = read.GeometricGmshMesh4("tests/geometric-mesh/simple_2D_coarse.msh", gmesh)
#gmesh = read.GeometricGmshMesh4("tests/geometric-mesh/TopeMeshPT.msh", gmesh)
gmesh.BuildConnectivity()


#Creacion de la malla computacional de flujo
dim = gmesh.Dimension();
cfluxmesh = TPZCompMesh(gmesh)
cfluxmesh.SetDefaultOrder(1)
val1 = TPZFMatrix(1,1,0.0)
val2 = TPZFMatrix(1,1,0.0)
#id dimension
volume = TPZMixedDarcyFlow(1,dim)
volume.SetPermeability(1.0);
cfluxmesh.InsertMaterialObject(volume)
N_Type = 1
D_Type = 0

nonfluxid = 2
bc_nonflux_mat = volume.CreateBC(volume, nonfluxid,N_Type,val1,val2)
cfluxmesh.InsertMaterialObject(bc_nonflux_mat)

inletid = 4
val3 = TPZFMatrix(1,1,30.0)
bc_inlet_mat = volume.CreateBC(volume, inletid,D_Type,val1,val3)
cfluxmesh.InsertMaterialObject(bc_inlet_mat)

val4 = TPZFMatrix(1,1,10.0)
outletid = 5
bc_outlet_mat = volume.CreateBC(volume, outletid,D_Type,val1,val4)
cfluxmesh.InsertMaterialObject(bc_outlet_mat)
cfluxmesh.SetDimModel(dim)
cfluxmesh.SetAllCreateFunctionsHDiv()
cfluxmesh.AutoBuild()
cfluxmesh.InitializeBlock()

#Creacion de la malla computacional de presion
cpressuremesh = TPZCompMesh(gmesh)
cpressuremesh.SetDefaultOrder(1)
cpressuremesh.InsertMaterialObject(volume)
cpressuremesh.SetDimModel(dim)
cpressuremesh.SetDefaultOrder(1)
cpressuremesh.SetAllCreateFunctionsContinuous()
cpressuremesh.ApproxSpace().CreateDisconnectedElements(1);
cpressuremesh.AutoBuild()
cpressuremesh.InitializeBlock()

#Creacion de la malla multifisica
multcmesh = TPZMultiphysicsCompMesh(gmesh)
multcmesh.InsertMaterialObject(volume)
multcmesh.InsertMaterialObject(bc_nonflux_mat)
multcmesh.InsertMaterialObject(bc_inlet_mat)
multcmesh.InsertMaterialObject(bc_outlet_mat)
multcmesh.SetDimModel(dim)
vecmesh = TPZManVecCompMesh(2)
vecmesh[0] = cfluxmesh
vecmesh[1] = cpressuremesh

active_aprox_space = TPZVecIntn(2)
active_aprox_space[0] =1
active_aprox_space[1] =1
multcmesh.BuildMultiphysicsSpace(active_aprox_space, vecmesh)

#Creacion del analisis y resolucion
an = TPZAnalysis(multcmesh, 1)
struc_mat = TPZSymetricSpStructMatrix(multcmesh)

an.SetStructuralMatrix(struc_mat)
stepsol = TPZStepSolver()
stepsol.SetDirect(ELDLt)
an.SetSolver(stepsol)
an.Assemble()
an.Solve()

name = str("resultadomult.vtk")
scalnames = TPZVecString(1)
vecnames = TPZVecString(1)
scalnames[0]="p"
vecnames[0]="q"

an.DefineGraphMesh(dim,scalnames,vecnames, name)
an.PostProcess(0,dim)

#Creacion de la malla de tranasporte
ctransportmesh = TPZCompMesh(gmesh)
matId=1
sol= TPZVecDouble(1,0.0)
voltransport = TPZL2Projection(matId, dim,1,sol)
ctransportmesh.InsertMaterialObject(voltransport)
ctransportmesh.SetDefaultOrder(0)
typ=0

val1 = TPZFMatrix(1,1,0.0)
val2 = TPZFMatrix(1,1,0.0)

bc_nonflux_mat = voltransport.CreateBC(voltransport, nonfluxid,typ,val1,val2)
ctransportmesh.InsertMaterialObject(bc_nonflux_mat)

inletid = 4
val3 = TPZFMatrix(1,1,1.0)
bc_inlet_mat = voltransport.CreateBC(voltransport, inletid,typ,val1,val3)
ctransportmesh.InsertMaterialObject(bc_inlet_mat)

val4 = TPZFMatrix(1,1,0.0)
outletid = 5
bc_outlet_mat = voltransport.CreateBC(voltransport, outletid,1,val1,val4)
ctransportmesh.InsertMaterialObject(bc_outlet_mat)

ctransportmesh.SetAllCreateFunctionsDiscontinuous()
ctransportmesh.AutoBuild()
ctransportmesh.InitializeBlock()

#Creacion de la malla multifisica de transporte
multtransportcmesh = TPZMultiphysicsCompMesh(gmesh)
mat_id=1
voltrans = TPZTracerFlow(mat_id, 0)
voltrans.SetPorosity(1.0)
multtransportcmesh.InsertMaterialObject(voltrans)
val1 = TPZFMatrix(1,1,0.0)
val2 = TPZFMatrix(1,1,1.0)
type_inlet = 0
type_outlet = 1
bc_1 =voltrans.CreateBC(voltrans, 2, type_inlet, val1,val2)
multtransportcmesh.InsertMaterialObject(bc_1)
bc_2 =voltrans.CreateBC(voltrans, 4, type_inlet, val1,val2)
multtransportcmesh.InsertMaterialObject(bc_2)
val2 = TPZFMatrix(1,1,0.0)
bc_3 =voltrans.CreateBC(voltrans, 5, type_outlet, val1,val2)
multtransportcmesh.InsertMaterialObject(bc_3)
multtransportcmesh.SetDimModel(dim)


vecmesht = TPZManVecCompMesh(3)
vecmesht[0] = cfluxmesh
vecmesht[1] = cpressuremesh
vecmesht[2] = ctransportmesh

active_aprox_spacet = TPZVecIntn(3)
active_aprox_spacet[0] =0
active_aprox_spacet[1] =0
active_aprox_spacet[2] =1
multtransportcmesh.BuildMultiphysicsSpace(active_aprox_spacet, vecmesht)
multtransportcmesh.Reference().ResetReference()
multtransportcmesh.LoadReferences()
nels = multtransportcmesh.NElements()
iel =0
while iel < nels :
    compel = multtransportcmesh.Element(iel)
    compel.CreateInterfaces()
    iel = iel+1


an_transp = TPZAnalysis(multtransportcmesh, 1)
struc_mat = TPZSpStructMatrix(multtransportcmesh)

an_transp.SetStructuralMatrix(struc_mat)
stepsol = TPZStepSolver()
stepsol.SetDirect(ELU)
an_transp.SetSolver(stepsol)


mass_matrix_Q = 1
mattrans = multtransportcmesh.FindMaterial(1)
mattrans.SetMassMatrixAssembly(mass_matrix_Q)




print("Computing mass matrix")
an_transp.Assemble()
print("Mass matrix is computed")


neq = an_transp.Mesh().NEquations()
MatrixAn = an_transp.GetMatrixSolver()


#copia los elementos de la diagonal
nfilas = MatrixAn.Rows()
MatrixDiagonal = TPZFMatrix(nfilas,1,0.0)
count = 0
while count < nfilas :
     MatrixDiagonal.SetItem(count, 0, MatrixAn.GetVal(count,count) )
     count = count+1

mass_matrix_Q = 0
dt = 10.0
n_steps = 25
mattrans.SetTimeStep(dt)
mattrans.SetMassMatrixAssembly(mass_matrix_Q)
an_transp.Assemble()
F_Inlet = an_transp.Rhs()

print("FInletFinalYORCH2")

print("Computing transport operator K=M+T, and F_Inlet")

s_n = TPZFMatrix(neq,1,0.0)
last_state_mass = TPZFMatrix(neq,1,0.0)
s_np1 = TPZFMatrix(neq,1,0.0)


count = 0
i=0
while count < n_steps:
    
    for i in range(neq):
        valortoset = MatrixDiagonal.GetVal(i,0)*s_n.GetVal(i,0)
        last_state_mass.SetItem(i,0, valortoset)
        
        
    
    print("valorafterset2: ")
    print(last_state_mass)
    an_transp.UpdateRhs(F_Inlet, last_state_mass)
    print("posprocesamos")

    count +=1
    an_transp.Solve()
    s_np1 = an_transp.Solution()
    print("SN =")
    print(s_n)
    print("SNP1 =")
    print(s_np1)
    an_transp.LoadSolution(s_np1)


    name = str("Transport.vtk")
    scalnames = TPZVecString(2)
    vecnames = TPZVecString(1)
    scalnames[0]="Sw"
    scalnames[1]="So"


    multtransportcmesh.FindMaterial(1).SetDimension(dim)
    an_transp.DefineGraphMesh(2,scalnames,vecnames, name)
    an_transp.PostProcess(0,2)

    s_n = s_np1





#mat = TPZMatPoisson3d(1,2)
#cmesh = TPZCompMesh(gmesh)
#val = cmesh.InsertMaterialObject(mat)
#
#D_Type = 0
#N_Type = 1
#val1 = TPZFMatrix(1,1,0.0)
#val2 = TPZFMatrix(1,1,0.0)
#bc_1 = mat.CreateBC(mat, 2,N_Type,val1,val2)
#cmesh.InsertMaterialObject(bc_1)
#
#val3 = TPZFMatrix(1,1,30.0)
#bc_2 = mat.CreateBC(mat, 4,D_Type,val1,val3)
#cmesh.InsertMaterialObject(bc_2)
#
#val4 = TPZFMatrix(1,1,10.0)
#bc_3 = mat.CreateBC(mat, 5,D_Type,val1,val4)
#cmesh.InsertMaterialObject(bc_3)
#
#cmesh.SetDimModel(2)
#
#cmesh.SetAllCreateFunctionsContinuous()
#cmesh.SetDefaultOrder(1)
#cmesh.AutoBuild()
#
#an = TPZAnalysis(cmesh, 1)
#struc_mat = TPZSymetricSpStructMatrix(cmesh)
#an.SetStructuralMatrix(struc_mat)
#stepsol = TPZStepSolver()
#stepsol.SetDirect(ECholesky)
#an.SetSolver(stepsol)
#an.Assemble()
#an.Solve()
#
#name = str("resultado.vtk")
#scalnames = TPZVecString(1)
#vecnames = TPZVecString(1)
#scalnames[0]="state"
#vecnames[0]="Flux"
#
#an.DefineGraphMesh(2,scalnames,vecnames, name)
#an.PostProcess(0,2)
