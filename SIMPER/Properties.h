#include <vector>
#include "PropertyDefinitions.h"

struct Properties
{
	SoilProperties Soil;
	FluidProperties Fluid;
	NonisothermalProperties Nonisothermal;
	BoundaryConditions BCs;
	GaussQuadrature GQ;
	SolutionVariables Solution;
	vector<int> PlotNodes;
	string GSLIBInputFile;
};
