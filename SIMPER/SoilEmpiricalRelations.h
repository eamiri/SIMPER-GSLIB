#ifndef __SIR__
#define __SIR__

#include <vector>
#include <math.h>
#include <stdlib.h>

using namespace std;

class SoilEmpiricalRelations
{
public:
	double HydraulicConductivity(double den, double maxDen, double minDen);
	double ThermalConductivity(double bd);	
	double ShiftBCgslib(double bc, double Tf, double Tl, double gslibCoeff, double stanDev, double latentHeat, double density, double heatCapacity);

private:
	double satHydCon;
	double satThermCon;

};

#endif