from neopz import *
import numpy as np

# Function to add boundary elements
def AddBoundaryElementsCook(gmesh):
	# Initializing leftset and rightset
	leftset = np.array([])
	rightset = np.array([])
	dim = gmesh.Dimension()
	nnodes = gmesh.NNodes()
	nelem = gmesh.NElements();
	for inodes in range(nnodes):
		nodevec = gmesh.NodeVec()
		xco = TPZVec_double(3)
		nodevec[inodes].GetCoordinates(xco)
		if (np.absolute(xco[0]) < 1e-3):
			leftset = np.append(leftset,[inodes])
		if (np.absolute(xco[0]-48) < 1e-3):
			rightset = np.append(rightset,[inodes])
	# Loop in elements for creating GeoElBC
	for el in range(nelem):
		gel = gmesh.Element(el)
		if (gel.Dimension() != dim-1):
			print("Error!")
		nsides = gel.NSides()
		for isides in range(nsides):
			if (gel.SideDimension(isides) != dim-1):
				continue
			nsidenodes = gel.NSideNodes(isides)
			nfoundleft = 0;
			nfoundright = 0;
			for isn in range(nsidenodes):
				nodeindex = gel.SideNodeIndex(isides,isn)
				if(nodeindex in leftset):
					nfoundleft=nfoundleft+1;
				if(nodeindex in rightset):
					nfoundright=nfoundright+1;
			if (nfoundright == nsidenodes):
				gelbc = TPZGeoElBC(gel,isides,EBC1)
			elif (nfoundleft == nsidenodes):
				gelbc = TPZGeoElBC(gel,isides,EBC2)
			else:
				gelside = TPZGeoElSide(gel,isides)
				neighbour = gelside.Neighbour()
				if (neighbour == gelside):
					gelbc = TPZGeoElBC(gelside,EBC3)

def InsertMaterialObjects3D(cmesh, scalarproblem):
	material = TPZElasticity3D(Emat1)
	material.SetMaterialDataHook(1000.0,0.49999)
	nstate = 3
	cmesh.InsertMaterialObject(material)
	val1 = TPZFMatrix(nstate,nstate,0)
	val2 = TPZFMatrix(nstate,1,0)
	bskeleton = material.CreateBC(material,ESkeleton,1,val1,val2)
	cmesh.InsertMaterialObject(bskeleton)
	val2.SetItem(1,0, 10.)
	print(val1)
	print(val2)
	bcond1 = material.CreateBC(material, EBC1, 1, val1, val2)
	cmesh.InsertMaterialObject(bcond1)
	val2.SetItem(1,0, 0.0)
	print(val1)
	print(val2)
	bcond2 = material.CreateBC(material,EBC2,0,val1,val2)
	cmesh.InsertMaterialObject(bcond2)
	bcond2 = material.CreateBC(material,EBC3,1,val1,val2)
	cmesh.InsertMaterialObject(bcond2)
	cmesh.AutoBuild()

# Data
Emat1 = 1;
EBC1 = 2;
EBC2 = 5;
EBC3 = 6;
ESkeleton = 13;

minporder = 1
maxporder = 2

minrefskeleton = 0
maxrefskeleton = 1

for porder in range(minporder,maxporder,1):
	for irefskeleton in range(minrefskeleton,maxrefskeleton,1):
		# Reading a UNSW file and creating the geometry
		elpartitions = TPZManVector_int64_t()
		scalingcenterindices = TPZVec_int64_t()
		gmesh = TPZGeoMesh()
		gmesh = gmesh.ReadUNSWSBGeoFile("tests/geometric-mesh/CooksMembrane_sbfemesh_16_1_1.txt",ESkeleton,elpartitions,scalingcenterindices)
		print("# Geometric Mesh ready!")
		# Creating the boundary conditions
		print("Creating Boundary Conditions")
		AddBoundaryElementsCook(gmesh)
		nelem = gmesh.NElements()
		elpartitions.Resize(nelem)
		matidtranslation = MapIntInt();
		matidtranslation[ESkeleton] = Emat1
		print("# Boundary Conditions created!")
		# Creating the computational mesh
		print("Building SBFEM Computational Mesh")
		build = TPZBuildSBFem(gmesh, ESkeleton, matidtranslation)
		build.SetPartitions(elpartitions,scalingcenterindices)
		print("# SBFEM Geometric Mesh ready!")
		# Priting the geometric mesh
		vtk = TPZVTKGeoMesh();
		print("Ploting VTK file to Geometric Mesh!")
		vtk.PrintGMeshVTK(gmesh,"GMeshVTK_SBFEM.vtk",1)
		print("# VTK file done!")
		print("Creating Computational Mesh")
		sbfem = TPZCompMesh(gmesh)
		sbfem.SetDefaultOrder(porder)
		scalarproblem = 0
		InsertMaterialObjects3D(sbfem, scalarproblem)
		build.BuildComputationalMeshFromSkeleton(sbfem)
		#print(sbfem)
		print("# Computational Mesh ready!")
		# Starting Analysis
		print("Starting Analysis")
		an = TPZAnalysis(sbfem,1)
		strmat = TPZSymetricSpStructMatrix(sbfem)
		# strmat.SetNumThreads(8)
		an.SetStructuralMatrix(strmat)
		step = TPZStepSolver()
		step.SetDirect(ECholesky)
		an.SetSolver(step)
		print("Entering Assemble...")
		an.Assemble()
		print("Assemble is done! Entering Solve...")
		print("Number of equations: ", sbfem.NEquations())
		an.Solve()
		print("# Solve is done! PostProcessing result...")
		# Post processing solution
		name = str("resultado.vtk")
		scalnames = TPZVec_string(3)
		vecnames = TPZVec_string(1)
		vecnames[0]="State"
		scalnames[0]="StressX"
		scalnames[1]="StressY"
		scalnames[2]="StressZ"
		an.DefineGraphMesh(3,scalnames,vecnames, name)
		an.PostProcess(1,3)
		print("# END OF SIMULATION! #")
