from neopz import*

def Sol_exact(x = TPZManVecDouble(), sol = TPZManVecDouble(), dsol = TPZFMatrix() ):
    
    dsol.Resize(3,3)
    sol.Resize(3)
    dsol.Zero()
    
    xv = x[0]
    yv = x[1]
    zv = x[2]
        
    Pi = math.pi
            
    v_x =  math.cos(2*Pi*yv)*math.sin(2*Pi*xv)
    v_y =  -(math.cos(2*Pi*xv)*math.sin(2*Pi*yv))
            
    sol[0]=v_x
    sol[1]=v_y
            
            # vx direction
    dsol.SetItem(0,0,2*Pi*math.cos(2*Pi*xv)*math.cos(2*Pi*yv))
    dsol.SetItem(0,1,2*Pi*math.sin(2*Pi*xv)*math.sin(2*Pi*yv))
            
            
            # vy direction
    dsol.SetItem(1,0,-2*Pi*math.sin(2*Pi*xv)*math.sin(2*Pi*yv))


gmesh = TPZGeoMesh()
read = TPZGmshReader()
gmesh = read.GeometricGmshMesh4("tests/geometric-mesh/simple_2D.msh", gmesh)
gmesh.BuildConnectivity()
mat = TPZMatLaplacian(1,2)
mat.SetForcingFunctionExact(Sol_exact)
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
an.SetExact2()

name = str("resultado.vtk")
scalnames = TPZVecString(1)
vecnames = TPZVecString(1)
#scalnames[0]="state"
scalnames[0]="ExactSolution"
vecnames[0]="Flux"

an.DefineGraphMesh(2,scalnames,vecnames, name)
an.PostProcess(0,2)




dsol.SetItem(1,1,-2*Pi*math.cos(2*Pi*xv)*math.cos(2*Pi*yv))




