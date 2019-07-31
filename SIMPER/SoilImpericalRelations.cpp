#include "SoilImpericalRelations.h"

double SoilImpericalRelations::HydraulicConductivity(double bd)
{
	//saturated hydraulic conductivity based on bulk density \cite Liu, H., and Lennartz, B. (2018). https://doi.org/10.1002/hyp.13314
	bd = bd / 1000.0;
	satHydCon = pow(10, (1.935 - 15.802 * bd + 19.552 * bd * bd)) / (3600.0 * 100.0); //[m/s]

	return satHydCon;
}

double SoilImpericalRelations::ThermalConductivity(double bd)
{
	satThermCon = (0.135 * bd + 64.7) / (2700 - 0.947 * bd); //[W/m/K]

	return satThermCon;
}