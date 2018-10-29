#ifndef __IBCONDITIONS__
#define __IBCONDITIONS__

#include <Eigen/Core>
#include <vector>
#include "Mesh.h"

using namespace std;
using namespace Eigen;

class IBConditions
{
public:
	IBConditions(Mesh *mesh);

	VectorXd Temp_0;
	VectorXd TempDot_0;
	VectorXd Residual;
	VectorXd Temp;
	VectorXd TempDot;
	vector<int> DirichletDof;
	vector<int> PlotNodes;

private:
	void SetBoundaryNodes();
	void SetInitialCondition();
	void SetBoundaryConditions();
	void SetPlotNodes();
	vector<int> oneDimensionalBC();
	vector<int> CornerFreezing();
	vector<int> DiskProblem();
	void DirichletBC(vector<int> boundaryNodes, double value);
	Mesh *MESH;
};

#endif 