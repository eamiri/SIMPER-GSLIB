#include "SimperInclude.h"

bool Inputs(string propsInputFile, string meshInputFile);
Properties InputProperties(string filePath);
Mesh InputMesh(string filePath);
void SgsimParameterFile();
void GSLIBRunSGSIM();
void UpscaleGSLIBtoSIMPER(string filePath);
void AddcoorParameterFile(int realizationNumber);

int noel;
int nond;
int ndoe;
int nGP;
Mesh MESH;
VectorXd GSLIBCoeffs;
MatrixXd GSLIBGrid;
Properties PROPS;
GaussPoints GP;

bool Inputs(string propsInputFile, string meshInputFile)
{
	PROPS = InputProperties(propsInputFile);
	MESH = InputMesh(meshInputFile);
	if (PROPS.Soil.IsGSLIB)
	{
		SgsimParameterFile();
		GSLIBRunSGSIM();
	}
	
	UpscaleGSLIBtoSIMPER("GSLIB/GSLIB_DATA.out");

	return true;
}

Properties InputProperties(string filePath)
{
	Properties props;
	ifstream propsFile;
	propsFile.open(filePath.c_str());
	if (propsFile.fail())
	{
		cout << "Cannot find file " << filePath << endl;
	}
	else
	{
		string line;

		getline(propsFile, line);
		getline(propsFile, line);
		getline(propsFile, line);
		istringstream isDataLine(line);
		isDataLine
			>> props.Soil.HeatCapacity
			>> props.Soil.ThermalConductivity
			>> props.Soil.Density
			>> props.Soil.Porosity
			>> props.Soil.ResidualWaterSaturation
			>> props.Soil.rSFC
			>> props.Soil.Wpar
			>> props.Soil.Mpar
			>> props.Soil.IsGSLIB;
		getline(propsFile, line);
		getline(propsFile, line);
		getline(propsFile, line);
		isDataLine.clear();
		isDataLine.str(line);
		isDataLine
			>> props.Fluid.SHeatCapacity
			>> props.Fluid.SThermalConductivity
			>> props.Fluid.SDensity
			>> props.Fluid.LHeatCapacity
			>> props.Fluid.LThermalConductivity
			>> props.Fluid.LDensity
			>> props.Fluid.LatentHeat;
		getline(propsFile, line);
		getline(propsFile, line);
		getline(propsFile, line);
		isDataLine.clear();
		isDataLine.str(line);
		isDataLine
			>> props.Nonisothermal.IsIsothermal
			>> props.Nonisothermal.HeatCapacity
			>> props.Nonisothermal.TempTransition
			>> props.Nonisothermal.ThermalConductivity
			>> props.Nonisothermal.TempSolid
			>> props.Nonisothermal.TempLiquid;
		getline(propsFile, line);
		getline(propsFile, line);
		getline(propsFile, line);
		isDataLine.clear();
		isDataLine.str(line);
		isDataLine >> props.BCs.InputFile;
		getline(propsFile, line);
		while (getline(propsFile, line))
		{
			isDataLine.clear();
			isDataLine.str(line);
			int bcNode;
			double bcValue;
			if (!(isDataLine >> bcNode >> bcValue))
			{
				break;
			}

			props.BCs.Node.push_back(bcNode);
			props.BCs.Value.push_back(bcValue);
		}

		getline(propsFile, line);
		isDataLine.clear();
		isDataLine.str(line);
		isDataLine >> props.GQ.NumberOfPoints;
		GP.GP(props.GQ.NumberOfPoints);
		nGP = props.GQ.NumberOfPoints;
		getline(propsFile, line);
		getline(propsFile, line);
		getline(propsFile, line);
		isDataLine.clear();
		isDataLine.str(line);
		isDataLine
			>> props.Solution.MaxTimestep
			>> props.Solution.DeltaTime
			>> props.Solution.TolPsi
			>> props.Solution.MaxIterations
			>> props.Solution.NewmarkGamma
			>> props.Solution.PlotInterval
			>> props.Solution.IsGMSH
			>> props.Solution.IsInputBC;
		getline(propsFile, line);
		getline(propsFile, line);
		getline(propsFile, line);
		isDataLine.clear();
		isDataLine.str(line);
		isDataLine
			>> props.GSLIB.CorrelationLength
			>> props.GSLIB.NumberOfCells
			>> props.GSLIB.GridSize
			>> props.GSLIB.NumberOfRealizations;
		/*getline(propsFile, line);
		while (getline(propsFile, line))
		{
			isDataLine.clear();
			isDataLine.str(line);
			int plotNode;
			if (!(isDataLine >> plotNode))
			{
				break;
			}

			props.PlotNodes.push_back(plotNode);
		}*/
	}

	return props;
}

void SgsimParameterFile()
{
	srand(time(0));
	long seedGSLIB = 305 * abs(rand());

	/*for (int i = 0; i < 5; ++i) {
		seedGSLIB = abs(seedGSLIB << 15) | (rand() & 0x7FFF);
	}*/

	FILE *inputFileGSLIB = fopen("GSLIB/SgsimInput.par", "w");

	fprintf(inputFileGSLIB, "Parameters for SGSIM\n\n");
	fprintf(inputFileGSLIB, "START OF PARAMETERS\n");
	fprintf(inputFileGSLIB, "nodata\n");
	fprintf(inputFileGSLIB, "1 2 0 3 5 0                             -columns for X,Y,Z,vr,wt,sec.var.\n");
	fprintf(inputFileGSLIB, "-1E+21 1E+21\n");
	fprintf(inputFileGSLIB, "0                                       -transform the data (0=no, 1=yes)\n");
	fprintf(inputFileGSLIB, "sgsim.trn\n");
	fprintf(inputFileGSLIB, "0                                       -consider ref. dist (0=no, 1=yes)\n");
	fprintf(inputFileGSLIB, "vmodel.var\n");
	fprintf(inputFileGSLIB, "1 2                                     -columns for vr and wt\n");
	fprintf(inputFileGSLIB, "0 15                                    -zmin,zmax (tail extrapolation)\n");
	fprintf(inputFileGSLIB, "1 0                                     -lower tail option\n");
	fprintf(inputFileGSLIB, "1 15                                    -upper tail option\n");
	fprintf(inputFileGSLIB, "0                                       -debug level (0-3)\n");
	fprintf(inputFileGSLIB, "GSLIB/SGSIM_nodata.dbg\n");
	fprintf(inputFileGSLIB, "GSLIB/SGSIM_output.out\n");
	fprintf(inputFileGSLIB, "%i                                       -number of realizations to generate\n", PROPS.GSLIB.NumberOfRealizations);
	fprintf(inputFileGSLIB, "%i 0 %e                              -nx, xmin, xsize\n", PROPS.GSLIB.NumberOfCells + 1, PROPS.GSLIB.GridSize);
	fprintf(inputFileGSLIB, "%i 0 %e                              -ny, ymin, ysize\n", PROPS.GSLIB.NumberOfCells + 1, PROPS.GSLIB.GridSize);
	fprintf(inputFileGSLIB, "1 0 1                                   -nz, zmin, zsize\n");
	fprintf(inputFileGSLIB, "%ld                               -random number seed\n", seedGSLIB);
	fprintf(inputFileGSLIB, "0 8                                     -Min and max original data for sim\n");
	fprintf(inputFileGSLIB, "12                                      -number of simulated nodes to use\n");
	fprintf(inputFileGSLIB, "1                                       -assign data to nodes (0=no, 1=yes)\n");
	fprintf(inputFileGSLIB, "1 3                                     -multiple grid search (0=no, 1=yes), num\n");
	fprintf(inputFileGSLIB, "0                                       -maximum data per octant (0=not used)\n");
	fprintf(inputFileGSLIB, "10 10 10                                -maximum search radii (hmax, hmin, vert)\n");
	fprintf(inputFileGSLIB, "0 0 0                                   -angles for search ellipsoid\n");
	fprintf(inputFileGSLIB, "51 51 11                                -size of covariance lookup table\n");
	fprintf(inputFileGSLIB, "0 0 1                                   -kType: 0=SK,1=OK,2=LVM,3=EXDR,4=COLC\n");
	fprintf(inputFileGSLIB, "nodata\n");
	fprintf(inputFileGSLIB, "4                                       -column\n");
	fprintf(inputFileGSLIB, "1 0                                     -nst, nugget NOFILE\n");
	fprintf(inputFileGSLIB, "2 1 0 0 0                               -it, cc, ang1, ang2, ang3\n");
	fprintf(inputFileGSLIB, "%e %e 1.6                           -a_hmax, a_hmin, a_vert\n", PROPS.GSLIB.CorrelationLength, PROPS.GSLIB.CorrelationLength);
	fflush(inputFileGSLIB);
}

void AddcoorParameterFile(int realizationNumber)
{
	FILE *inputFileGSLIB = fopen("GSLIB/AddcoorInput.par", "w");

	fprintf(inputFileGSLIB, "Parameters for SGSIM\n\n");
	fprintf(inputFileGSLIB, "START OF PARAMETERS\n");
	fprintf(inputFileGSLIB, "GSLIB/SGSIM_output.out\n");
	fprintf(inputFileGSLIB, "GSLIB/ADDCOOR_output.out\n");
	fprintf(inputFileGSLIB, "%i                                   -realization number to add coordinate\n", realizationNumber);
	fprintf(inputFileGSLIB, "%i 0 %e                  -nx, xmin, xsize\n", PROPS.GSLIB.NumberOfCells + 1, PROPS.GSLIB.GridSize);
	fprintf(inputFileGSLIB, "%i 0 %e                  -ny, ymin, ysize\n", PROPS.GSLIB.NumberOfCells + 1, PROPS.GSLIB.GridSize);
	fprintf(inputFileGSLIB, "1 0 1                               -nz, zmin, zsize\n");
	fflush(inputFileGSLIB);
}

void GSLIBRunSGSIM()
{
	cout << endl << "=== GSLIB UNCONDITIONAL SIMULATION ===" << endl;
	int ExecuteGSLIB = system("GSLIB\\GSLIBSimulation.exe \"GSLIB/SgsimInput.par\""); // simulation
	AddcoorParameterFile(1);
	ExecuteGSLIB = system("GSLIB\\GSLIBAddCoordinates.exe \"GSLIB/AddcoorInput.par\""); // adding coordinates
	cout << "=== END GSLIB UNCONDITIONAL SIMULATION ===" << endl;

	FILE *inputFileGSLIB = fopen("../Results/GSLIB_SIMULATION.plt", "w");
	GSLIBCoeffs.resize((PROPS.GSLIB.NumberOfCells + 1) * (PROPS.GSLIB.NumberOfCells + 1));
	GSLIBCoeffs.setZero();
	GSLIBGrid.resize((PROPS.GSLIB.NumberOfCells + 1) * (PROPS.GSLIB.NumberOfCells + 1), 2);
	GSLIBGrid.setZero();
	string AddcoorOutputFilePath = "GSLIB/ADDCOOR_output.out";
	ifstream gslibFile;
	gslibFile.open(AddcoorOutputFilePath.c_str());
	string line;
	getline(gslibFile, line);
	getline(gslibFile, line);
	getline(gslibFile, line);
	getline(gslibFile, line);
	getline(gslibFile, line);
	int index = 0;
	double gslibCoeff;
	double xCoord;
	double yCoord;
	double zCoord;
	while (getline(gslibFile, line))
	{
		if (!(gslibFile >> xCoord >> yCoord >> zCoord >> gslibCoeff))
		{
			break;
		}

		GSLIBGrid(0, index) = xCoord;
		GSLIBGrid(1, index) = yCoord;
		GSLIBCoeffs(index) = gslibCoeff;
	}

}

void UpscaleGSLIBtoSIMPER(string filePath)
{
	if (PROPS.Soil.IsGSLIB)
	{
		double maxGslibCoeff = GSLIBCoeffs.maxCoeff();
		double minGslibCoeff = GSLIBCoeffs.minCoeff();
		for (int n = 0; n < MESH.NumberOfNodes; n++)
		{
			double xNode = MESH.Nodes[n].Coordinates.x;
			double yNode = MESH.Nodes[n].Coordinates.y;
			double gslibCoeff;
			double distanceP = INFINITY;
			for (int i = 0; i < (PROPS.GSLIB.NumberOfCells + 1) * (PROPS.GSLIB.NumberOfCells + 1); i++)
			{
				double distance = sqrt((xNode - GSLIBGrid(i, 0)) * (xNode - GSLIBGrid(i, 0)) + (yNode - GSLIBGrid(i, 1)) * (yNode - GSLIBGrid(i, 1)));
				if (distance < distanceP)
				{
					distanceP = distance;
					gslibCoeff = GSLIBCoeffs(i);
				}
			}

			if (gslibCoeff > 0)
			{
				MESH.Nodes[n].GSLIBCoeff = gslibCoeff / maxGslibCoeff;
			}
			else
			{
				MESH.Nodes[n].GSLIBCoeff = gslibCoeff / abs(minGslibCoeff);
			}

			break;
		}

		double GSLIBCoeffE;
		VectorXi elementDofs;
		VectorXd GSLIBCoeffsNode;
		FILE *plotHeatCapacity = fopen("../Results/HeatCapacity.plt", "w");
		for (int e = 0; e < MESH.NumberOfElements; e++)
		{
			elementDofs = MESH.GetElementDofs(e, ndoe);
			GSLIBCoeffsNode = MESH.GetNodalValues(GSLIBCoeffs, elementDofs);
			GSLIBCoeffE = 0;
			for (int ig = 1; ig < 4; ig++)
			{
				GSLIBCoeffE += abs(1 + GSLIBCoeffsNode(ig));
			}

			MESH.Elements[e].SoilHeatCapacity = PROPS.Soil.HeatCapacity * GSLIBCoeffE / 4;

			fprintf(plotHeatCapacity, "variables =\"X\" \"Y\" \"<i>c</i><sub>soil</sub>\"\n");
			fprintf(plotHeatCapacity, "ZONE N = %5.0d, E = %5.0d, ZONETYPE = FEQuadrilateral, DATAPACKING = POINT\n", ndoe, 1);
			for (int inod = 0; inod < ndoe; inod++)
			{
				fprintf(plotHeatCapacity, "%e\t%e\t%e\n", MESH.Elements[e].Nodes[inod].Coordinates.x, MESH.Elements[e].Nodes[inod].Coordinates.y, MESH.Elements[e].SoilHeatCapacity);
			}

			fprintf(plotHeatCapacity, "1 2 3 4");
			fflush(plotHeatCapacity);
		}
	}
	else
	{
		cout << "Homogeneous Media" << endl;
		for (int e = 0; e < MESH.NumberOfElements; e++)
		{
			MESH.Elements[e].SoilHeatCapacity = PROPS.Soil.HeatCapacity;
		}
	}	
}

Mesh InputMesh(string filePath)
{
	Mesh mesh;
	Node node;
	Element element;

	double *maxX = &mesh.MaxX;
	double *minX = &mesh.MinX;
	double *maxY = &mesh.MaxY;
	double *minY = &mesh.MinY;
	double *maxZ = &mesh.MaxZ;
	double *minZ = &mesh.MinZ;

	*maxX = -1.0* INFINITY;
	*minX = INFINITY;
	*maxY = -1.0 * INFINITY;
	*minY = INFINITY;
	*maxZ = -1.0 * INFINITY;
	*minZ = INFINITY;

	if (PROPS.Solution.IsGMSH)
	{
		ifstream meshFile;
		meshFile.open(filePath.c_str());
		if (meshFile.fail())
		{
			cout << "Cannot find file " << filePath << endl;
		}
		else
		{
			string line;

			getline(meshFile, line);
			getline(meshFile, line);
			getline(meshFile, line);
			getline(meshFile, line);
			getline(meshFile, line);

			istringstream isSurfaceNumber(line);
			double nSurface;

			isSurfaceNumber >> nSurface;

			for (int i = 0; i < nSurface; i++)
			{
				getline(meshFile, line);
			}

			getline(meshFile, line);
			getline(meshFile, line);
			getline(meshFile, line);

			// Nodes Data
			while (getline(meshFile, line))
			{
				istringstream isNodeData(line);
				int n;
				double x, y, z;

				if (!(isNodeData >> n >> x >> y >> z))
				{
					break;
				}

				node.n = n;
				node.Coordinates.x = x;
				node.Coordinates.y = y;
				node.Coordinates.z = z;

				*minX = (x < *minX) ? x : *minX;
				*maxX = (x > *maxX) ? x : *maxX;
				*minY = (y < *minY) ? y : *minY;
				*maxY = (y > *maxY) ? y : *maxY;
				*minZ = (z < *minZ) ? z : *minZ;
				*maxZ = (z > *maxZ) ? z : *maxZ;

				mesh.Nodes.push_back(node);
			}

			getline(meshFile, line);
			getline(meshFile, line);

			//Elements Data
			while (getline(meshFile, line))
			{
				istringstream isElementData(line);
				int e, X1, X2, X3, X4, n1, n2, n3, n4;

				if (!(isElementData >> e >> X1 >> X2 >> X3 >> X4 >> n1 >> n2 >> n3 >> n4))
				{
					break;
				}

				element.Nodes.clear();
				element.e = e;
				node.n = n1;
				node.Coordinates = mesh.Nodes[n1 - 1].Coordinates;
				element.Nodes.push_back(node);
				node.n = n2;
				node.Coordinates = mesh.Nodes[n2 - 1].Coordinates;
				element.Nodes.push_back(node);
				node.n = n3;
				node.Coordinates = mesh.Nodes[n3 - 1].Coordinates;
				element.Nodes.push_back(node);
				node.n = n4;
				node.Coordinates = mesh.Nodes[n4 - 1].Coordinates;
				element.Nodes.push_back(node);

				mesh.Elements.push_back(element);
			}

			mesh.NumberOfElements = mesh.Elements.size();
			mesh.NumberOfNodes = mesh.Nodes.size();
			mesh.ElementNumberOfNodes = 4;
			noel = mesh.NumberOfElements;
			nond = mesh.NumberOfNodes;
			ndoe = mesh.ElementNumberOfNodes;
		}
	}
	else
	{
		ifstream meshFile;
		meshFile.open(filePath.c_str());
		if (meshFile.fail())
		{
			cout << "Cannot find file " << filePath << endl;
		}
		else
		{
			string line;

			getline(meshFile, line);
			getline(meshFile, line);

			// Nodes Data
			while (getline(meshFile, line))
			{
				istringstream isNodeData(line);
				int n;
				double x, y, z;

				if (!(isNodeData >> n >> x >> y >> z))
				{
					break;
				}

				node.n = n;
				node.Coordinates.x = x;
				node.Coordinates.y = y;
				node.Coordinates.z = z;

				*minX = (x < *minX) ? x : *minX;
				*maxX = (x > *maxX) ? x : *maxX;
				*minY = (y < *minY) ? y : *minY;
				*maxY = (y > *maxY) ? y : *maxY;
				*minZ = (z < *minZ) ? z : *minZ;
				*maxZ = (z > *maxZ) ? z : *maxZ;

				mesh.Nodes.push_back(node);
			}

			getline(meshFile, line);
			getline(meshFile, line);

			//Elements Data
			while (getline(meshFile, line))
			{
				istringstream isElementData(line);
				int e, n1, n2, n3, n4;

				if (!(isElementData >> e >> n1 >> n2 >> n3 >> n4))
				{
					break;
				}

				element.Nodes.clear();
				element.e = e;
				node.n = n1;
				node.Coordinates = mesh.Nodes[n1 - 1].Coordinates;
				element.Nodes.push_back(node);
				node.n = n2;
				node.Coordinates = mesh.Nodes[n2 - 1].Coordinates;
				element.Nodes.push_back(node);
				node.n = n3;
				node.Coordinates = mesh.Nodes[n3 - 1].Coordinates;
				element.Nodes.push_back(node);
				node.n = n4;
				node.Coordinates = mesh.Nodes[n4 - 1].Coordinates;
				element.Nodes.push_back(node);

				mesh.Elements.push_back(element);
			}

			mesh.NumberOfElements = mesh.Elements.size();
			mesh.NumberOfNodes = mesh.Nodes.size();
			mesh.ElementNumberOfNodes = 4;
			noel = mesh.NumberOfElements;
			nond = mesh.NumberOfNodes;
			ndoe = mesh.ElementNumberOfNodes;
		}		
	}

	return mesh;
}