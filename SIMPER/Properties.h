#ifndef __PROPERTIES__
#define __PROPERTIES__

#include <vector>
#include "PropertyDefinitions.h"

struct Properties
{
	SoilProperties Soil;
	FluidProperties Fluid;
	NonisothermalProperties Nonisothermal;
	FenProperties Fen;
	BoundaryConditions BCs;
	GaussQuadrature GQ;
	SolutionVariables Solution;
	vector<int> PlotNodes;
	GSLIBProps GSLIB;
};

#endif 