#include "IBConditions.h"

IBConditions::IBConditions (Mesh *m, Properties *p)
{
	MESH = &(*m);
	PROPS = &(*p);
	SetInitialCondition();
	SetBoundaryConditions();
	InputBC("../Inputs/" + PROPS->BCs.InputFile);
}

void IBConditions::SetBoundaryNodes(vector<int> boundaryNodes)
{
	for (int iBC = 0; iBC < boundaryNodes.size(); iBC++)
	{
		MESH->Nodes[boundaryNodes[iBC]].IsBoundary = true;
	}
}

vector<int> IBConditions::FindDirichletBoundaryNodes()
{
	vector<int> dirDOF;
	for (int n = 0; n < MESH->NumberOfNodes; n++)
	{
		double yNode = MESH->Nodes[n].Coordinates.y;
		double xNode = MESH->Nodes[n].Coordinates.x;
		if ((yNode == MESH->MinY) || (yNode == MESH->MaxY))
		{
			dirDOF.push_back(n);
		}
	}

	return dirDOF;
}

vector<int> IBConditions::FindTopBoundaryNodes()
{
	vector<int> dirDOF;
	for (int n = 0; n < MESH->NumberOfNodes; n++)
	{
		double yNode = MESH->Nodes[n].Coordinates.y;
		double xNode = MESH->Nodes[n].Coordinates.x;
		if (yNode == MESH->MaxY)
		{
			dirDOF.push_back(n);
		}
	}

	return dirDOF;
}

vector<int> IBConditions::FindBottomBoundaryNodes()
{
	vector<int> dirDOF;
	for (int n = 0; n < MESH->NumberOfNodes; n++)
	{
		double yNode = MESH->Nodes[n].Coordinates.y;
		double xNode = MESH->Nodes[n].Coordinates.x;
		if (yNode == MESH->MinY)
		{
			dirDOF.push_back(n);
		}
	}

	return dirDOF;
}

vector<int> IBConditions::DiskProblem()
{
	vector<int> dirDOF;
	double minRadius = 1.0E+21;
	double yNode;
	double xNode;
	double radius;

	for (int n = 0; n < MESH->NumberOfNodes; n++)
	{
		yNode = MESH->Nodes[n].Coordinates.y;
		xNode = MESH->Nodes[n].Coordinates.x;
		radius = sqrt(xNode * xNode + yNode * yNode);

		if (radius < minRadius)
		{
			minRadius = sqrt(xNode * xNode + yNode * yNode);
		}
	}

	for (int n = 0; n < MESH->NumberOfNodes; n++)
	{
		yNode = MESH->Nodes[n].Coordinates.y;
		xNode = MESH->Nodes[n].Coordinates.x;
		radius = sqrt(xNode * xNode + yNode * yNode);

		if (radius < 1.01*minRadius)
		{
			dirDOF.push_back(n);
		}
	}

	return dirDOF;
}

vector<int> IBConditions::CornerFreezing()
{
	vector<int> dirDOF;
	for (int n = 0; n < MESH->NumberOfNodes; n++)
	{
		double xNode = MESH->Nodes[n].Coordinates.x;
		double yNode = MESH->Nodes[n].Coordinates.y;
		if (yNode == MESH->MinY)
		{
			dirDOF.push_back(n);
		}
		else if (xNode == MESH->MinX)
		{
			dirDOF.push_back(n);
		}
	}

	return dirDOF;
}

void IBConditions::SetPlotNodes()
{
	for (int n = 0; n < MESH->NumberOfNodes; n++)
	{
		double xNode = MESH->Nodes[n].Coordinates.x;
		double yNode = MESH->Nodes[n].Coordinates.y;
		if (yNode == xNode)
		{
			PlotNodes.push_back(MESH->Nodes[n].n);
		}
	}
}

void IBConditions::SetDirichletBC(vector<int> boundaryNodes, double value)
{
	for (int n = 0; n < boundaryNodes.size(); n++)
	{
		Temp(boundaryNodes[n]) = value;
	}
}

void IBConditions::InputBC(string filePath)
{
	BCInputData.resize(PROPS->Solution.MaxTimestep, 2);
	BCInputData.setZero();
	if (PROPS->Solution.IsInputBC)
	{
		ifstream bcInputFile;
		bcInputFile.open(filePath.c_str());
		if (bcInputFile.fail())
		{
			cout << "Cannot find file " << filePath << endl;
		}
		else
		{
			cout << "BC input file: " << filePath << endl;
			string line;
			int bcTimstep;
			double bcTemperature;
			getline(bcInputFile, line);
			getline(bcInputFile, line);
			istringstream isDataLine(line);
			isDataLine 
				>> DeltaTimeBC 
				>> EndTimestepBC;
			getline(bcInputFile, line);
			int lastIBC = 0;
			for (int ibc = 0; ibc < PROPS->Solution.MaxTimestep; ibc++)
			{
				getline(bcInputFile, line);
				isDataLine.clear();
				isDataLine.str(line);
				isDataLine >> bcTimstep >> bcTemperature;
				if (!(ibc < EndTimestepBC))
				{
					BCInputData(ibc, 0) = BCInputData(ibc - lastIBC, 0);
					BCInputData(ibc, 1) = BCInputData(ibc - lastIBC, 1);
				}
				else
				{
					BCInputData(ibc, 0) = bcTimstep;
					BCInputData(ibc, 1) = bcTemperature;
					lastIBC++;
				}
			}
		}
	}
}

void IBConditions::SetInitialCondition()
{
	Temp_0 = -0.5 * Temp_0.Ones(MESH->NumberOfNodes);
	TempDot_0 = TempDot_0.Zero(MESH->NumberOfNodes);
	Residual = Residual.Zero(MESH->NumberOfNodes);
	Temp = Temp_0;
	TempDot = TempDot_0;
}

void IBConditions::SetBoundaryConditions()
{
	TopBC = FindTopBoundaryNodes();
	BottomBC = FindBottomBoundaryNodes();
	DirichletBC = TopBC;
	//DirichletBC = FindDirichletBoundaryNodes();
	SetBoundaryNodes(DirichletBC);
	//DirichletBC(DirichletDof, 1);
	SetPlotNodes();
}