#include "SoilEmpiricalRelations.h"

double SoilEmpiricalRelations::HydraulicConductivity(double bd)
{
	//saturated hydraulic conductivity based on bulk density \cite Liu, H., and Lennartz, B. (2018). https://doi.org/10.1002/hyp.13314
	if (bd > 700.0) // **** this empirical realtion doesn't work for densities greater than 700 since it's just for peat soil ****
	{
		bd = 700.0; 
	}

	bd = bd / 1000.0;
	satHydCon = pow(10.0, (1.935 - 15.802 * bd + 19.552 * bd * bd)) / (3600.0 * 100.0); //[m/s]

	return satHydCon;
}

double SoilEmpiricalRelations::ThermalConductivity(double bd)
{
	satThermCon = (0.135 * bd + 64.7) / (2700 - 0.947 * bd); //[W/m/K] //CLM manual eqn. 6.62

	return satThermCon;
}

double SoilEmpiricalRelations::ShiftBCgslib(double bc, double Tf, double Tl, double gslibCoeff, double stanDev, double latentHeat, double density, double heatCapacity)
{
	double enthalpy;
	double shiftedTemp = 0.0;
	if (bc > Tl)
	{
		enthalpy = density * heatCapacity * bc + latentHeat * Tl / abs(Tl - Tf);
	}
	else if (bc < Tf)
	{
		enthalpy = density * heatCapacity * bc + latentHeat * Tf / abs(Tl - Tf);
	}
	else
	{
		enthalpy = density * heatCapacity * bc + latentHeat * bc / abs(Tl - Tf);
	}

	enthalpy += density * heatCapacity * stanDev * gslibCoeff;
	if (enthalpy > (density * heatCapacity * Tl + latentHeat * Tl / abs(Tl - Tf)))
	{
		shiftedTemp = (enthalpy - latentHeat * Tl / abs(Tl - Tf)) / (density * latentHeat);
	}
	else if (enthalpy < (density * heatCapacity * Tf + latentHeat * Tf / abs(Tl - Tf)))
	{
		shiftedTemp = (enthalpy - latentHeat * Tf / abs(Tl - Tf)) / (density * latentHeat);
	}
	else
	{
		shiftedTemp = enthalpy / (density * heatCapacity + latentHeat / abs(Tl - Tf));
	}

	return shiftedTemp;
}