#include "IBConditions.h"

IBConditions::IBConditions (Mesh *m)
{
	MESH = &(*m);
	SetInitialCondition();
	SetBoundaryConditions();
}

void IBConditions::SetBoundaryNodes()
{
	for (int iBC = 0; iBC < DirichletDof.size(); iBC++)
	{
		MESH->Nodes[DirichletDof[iBC]].IsBoundary = true;
	}
}

vector<int> IBConditions::oneDimensionalBC()
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

void IBConditions::DirichletBC(vector<int> boundaryNodes, double value)
{
	for (int n = 0; n < boundaryNodes.size(); n++)
	{
		Temp(boundaryNodes[n]) = value;
	}
}

void IBConditions::SetInitialCondition()
{
	Temp_0 = -3.0 * Temp_0.Ones(MESH->NumberOfNodes);
	TempDot_0 = TempDot_0.Zero(MESH->NumberOfNodes);
	Residual = Residual.Zero(MESH->NumberOfNodes);
	Temp = Temp_0;
	TempDot = TempDot_0;
}

void IBConditions::SetBoundaryConditions()
{
	DirichletDof = oneDimensionalBC();
	SetBoundaryNodes();
	DirichletBC(DirichletDof, 5.0);
	SetPlotNodes();
}