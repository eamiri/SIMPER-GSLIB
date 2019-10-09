#ifndef __PROPERTIES__
#define __PROPERTIES__

#include <vector>
#include "PropertyDefinitions.h"

struct Properties
{
	int nSoilType;
	bool IsLateralFlow;
	vector<SoilProperties> Soil;
	FluidProperties Fluid;
	NonisothermalProperties Nonisothermal;
	BoundaryConditions BCs;
	GaussQuadrature GQ;
	SolutionVariables Solution;
	vector<int> PlotNodes;
	GSLIBProps GSLIB;
};

#endif 