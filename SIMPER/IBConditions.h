#ifndef __IBCONDITIONS__
#define __IBCONDITIONS__

#include <Eigen/Core>
#include <vector>
#include "Mesh.h"
#include "Properties.h"
#include <iomanip>
#include <iostream>
#include <fstream>
#include <string>
#include <strstream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <cstring>

using namespace std;
using namespace Eigen;

class IBConditions
{
public:
	IBConditions(Mesh *mesh, Properties *props);

	VectorXd Temp_0;
	VectorXd TempDot_0;
	VectorXd Residual;
	VectorXd Temp;
	VectorXd TempDot;
	vector<int> TopBC;
	vector<int> BottomBC;
	vector<int> DirichletBC;
	vector<int> PlotNodes;
	MatrixXd BCInputData;
	double DeltaTimeBC;
	int  EndTimestepBC;

private:
	void SetBoundaryNodes(vector<int> boundaryNodes);
	void SetInitialCondition();
	void SetBoundaryConditions();
	void SetPlotNodes();
	vector<int> FindDirichletBoundaryNodes();
	vector<int> FindTopBoundaryNodes();
	vector<int> FindBottomBoundaryNodes();
	vector<int> CornerFreezing();
	vector<int> DiskProblem();
	void SetDirichletBC(vector<int> boundaryNodes, double value);
	void InputBC(string filePath);
	Mesh *MESH;
	Properties *PROPS;
};

#endif 