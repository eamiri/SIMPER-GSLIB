#include "IBConditions.h"

IBConditions::IBConditions (Mesh *m, Properties *p)
{
	MESH = &(*m);
	PROPS = &(*p);
	SetInitialCondition();
	SetBoundaryConditions();
	InputBC("../Inputs/" + PROPS->BCs.BCInputFile);
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
		//if ((yNode == MESH->MinY) || (yNode == MESH->MaxY))
		if (yNode == MESH->MaxY)
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
		double yNode = 10.0 - MESH->Nodes[n].Coordinates.y;
		if (abs(xNode - 0.005) <= 1E-5)
		{
			if    (abs(yNode - 0.00) <= 1E-5
				|| abs(yNode - 0.05) <= 1E-5
				|| abs(yNode - 0.10) <= 1E-5
				|| abs(yNode - 0.15) <= 1E-5
				|| abs(yNode - 0.20) <= 1E-5
				|| abs(yNode - 0.25) <= 1E-5
				|| abs(yNode - 0.30) <= 1E-5
				|| abs(yNode - 0.35) <= 1E-5
				|| abs(yNode - 0.45) <= 1E-5
				|| abs(yNode - 0.55) <= 1E-5
				|| abs(yNode - 0.65) <= 1E-5)
			{
				PlotNodes.push_back(MESH->Nodes[n].n);
			}
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
	BCInputData.resize(PROPS->Solution.MaxTimestep, 3);
	BCInputData.setZero();
	if (PROPS->BCs.isBCInput)
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
			double temperatureSD;
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
				isDataLine >> bcTimstep >> bcTemperature >> temperatureSD;
				if (!(ibc < EndTimestepBC))
				{
					BCInputData(ibc, 0) = BCInputData(ibc - lastIBC, 0);
					BCInputData(ibc, 1) = BCInputData(ibc - lastIBC, 1);
					BCInputData(ibc, 2) = BCInputData(ibc - lastIBC, 2);
				}
				else
				{
					BCInputData(ibc, 0) = bcTimstep;
					BCInputData(ibc, 1) = bcTemperature;
					BCInputData(ibc, 2) = temperatureSD;
					lastIBC++;
				}
			}
		}
	}
}

void IBConditions::SetInitialCondition()
{
	if (PROPS->BCs.isICInput)
	{
		InputIC("../Inputs/" + PROPS->BCs.ICInputFile);
	}
	else
	{
		Temp_0 = PROPS->BCs.UniformIC * Temp_0.Ones(MESH->NumberOfNodes);
		TempDot_0 = TempDot_0.Zero(MESH->NumberOfNodes);
		Residual = Residual.Zero(MESH->NumberOfNodes);
		Temp = Temp_0;
		TempDot = TempDot_0;
	}
}

void IBConditions::InputIC(string filePath)
{
	ifstream bcInputFile;
	bcInputFile.open(filePath.c_str());
	if (bcInputFile.fail())
	{
		cout << "Cannot find file " << filePath << endl;
	}
	else
	{
		cout << "IC input file: " << filePath << endl;
		string line;
		istringstream isDataLine(line);
		getline(bcInputFile, line);
		int iNode;
		double icTemperature;
		Temp_0.resize(MESH->NumberOfNodes);
		int lastIBC = 0;
		for (int iic = 0; iic < MESH->NumberOfNodes; iic++)
		{
			getline(bcInputFile, line);
			isDataLine.clear();
			isDataLine.str(line);
			isDataLine >> iNode >> icTemperature;
			Temp_0(iNode-1) = icTemperature;
			TempDot_0 = TempDot_0.Zero(MESH->NumberOfNodes);
			Residual = Residual.Zero(MESH->NumberOfNodes);
			Temp = Temp_0;
			TempDot = TempDot_0;
		}
	}
}

void IBConditions::SetBoundaryConditions()
{
	TopBC = FindTopBoundaryNodes();
	BottomBC = FindBottomBoundaryNodes();
	//DirichletBC = TopBC;
	DirichletBC = FindDirichletBoundaryNodes();
	SetBoundaryNodes(DirichletBC);
	//DirichletBC(DirichletDof, 1);
	SetPlotNodes();
}