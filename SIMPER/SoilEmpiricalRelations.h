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

private:
	double satHydCon;
	double satThermCon;

};

#endif