#ifndef __SIR__
#define __SIR__

#include <vector>
#include <math.h>

using namespace std;

class SoilEmpiricalRelations
{
public:
	double HydraulicConductivity(double bd);
	double ThermalConductivity(double bd);	
	double ShiftBCgslib(double bc, double Tf, double Tl, double gslibCoeff, double stanDev, double latentHeat, double density, double heatCapacity);

private:
	double satHydCon;
	double satThermCon;

};

#endif