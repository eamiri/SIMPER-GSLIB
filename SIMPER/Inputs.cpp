#include "SimperInclude.h"

bool Inputs(string propsInputFile, string meshInputFile);
Properties InputProperties(string filePath);
void GSLIBInput(string filePath);
Mesh InputMesh(string filePath);

int noel;
int nond;
int ndoe;
int nGP;
Mesh MESH;
VectorXd GSLIBCoeffs;
Properties PROPS;
GaussPoints GP;

bool Inputs(string propsInputFile, string meshInputFile)
{
	PROPS = InputProperties(propsInputFile);
	MESH = InputMesh(meshInputFile);
	GSLIBInput("../Inputs/" + PROPS.GSLIBInputFile);

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
		isDataLine.clear();
		isDataLine.str(line);
		isDataLine
			>> props.GSLIBInputFile;
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

void GSLIBInput(string filePath)
{
	if (PROPS.Soil.IsGSLIB)
	{
		GSLIBCoeffs.resize(MESH.NumberOfNodes);
		GSLIBCoeffs.setZero();
		ifstream gslibFile;
		gslibFile.open(filePath.c_str());
		if (gslibFile.fail())
		{
			cout << "Cannot find file " << filePath << endl;
		}
		else
		{
			cout << "Heterogeneous Media" << endl;
			cout << "GSLIB input file: " << filePath << endl;
			string line;
			getline(gslibFile, line);
			getline(gslibFile, line);
			getline(gslibFile, line);
			for (int i = 0; i < MESH.NumberOfNodes; i++)
			{
				istringstream isGSLIBData(line);
				double gslibCoeff;
				if (!(gslibFile >> gslibCoeff))
				{
					break;
				}

				GSLIBCoeffs(i) = gslibCoeff;
			}
		}

		double maxGslibCoeff = GSLIBCoeffs.maxCoeff();
		double minGslibCoeff = GSLIBCoeffs.maxCoeff();

		for (int n = 0; n < MESH.NumberOfNodes; n++)
		{
			if (GSLIBCoeffs[n] > 0)
			{
				MESH.Nodes[n].GSLIBCoeff = GSLIBCoeffs[n] / maxGslibCoeff;
			}
			else
			{
				MESH.Nodes[n].GSLIBCoeff = GSLIBCoeffs[n] / abs(minGslibCoeff);
			}
		}

		double GSLIBCoeffE;
		VectorXi elementDofs;
		VectorXd GSLIBCoeffsNode;
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