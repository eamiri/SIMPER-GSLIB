#ifndef __SIR__
#define __SIR__

#include <vector>

using namespace std;

class SoilImpericalRelations
{
public:
	double HydraulicConductivity(double bd);
	double ThermalConductivity(double bd);	

private:
	double satHydCon;
	double satThermCon;

};



#endif