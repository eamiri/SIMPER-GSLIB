#ifndef __SATFUNCTIONs__
#define __SATFUNCTIONS__

class SaturationFunctions
{
public:
	SaturationFunctions(double Tgau, double Tsol, double Tliq, double Sres);

	double Swat;
	double dSwat;
	double ISwat;
	double ISwatxT;
	double IdSwatxT;

	double Sice;
	double dSice;
	double ISice;
	double ISicexT;
	double IdSicexT;

};

#endif 



