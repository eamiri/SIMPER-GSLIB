#include "SimperInclude.h"

SaturationFunctions::SaturationFunctions(double Tgau, double Tsol, double Tliq, double Sres)
{
	double pi = PI;

	Swat =

		Sres - (sin(pi*(Tgau + 1 / 2)) / 2 - 1 / 2)*(sin(pi*(Tgau + 1 / 2)) / 2 + 1 / 2)*(heaviside(Tgau + 1) - heaviside(Tgau));


	dSwat =

		-(pi*cos(pi*(Tgau + 1 / 2))*(sin(pi*(Tgau + 1 / 2)) / 2 - 1 / 2)*(heaviside(Tgau + 1) - heaviside(Tgau))) / 2 - (pi*cos(pi
			*(Tgau + 1 / 2))*(sin(pi*(Tgau + 1 / 2)) / 2 + 1 / 2)*(heaviside(Tgau + 1) - heaviside(Tgau))) / 2;


	ISwat =

		Sres * (Tgau + 1) - heaviside(Tgau)*(Tgau / 8 + sin(2 * pi*(Tgau + 1 / 2)) / (16 * pi)) + heaviside(Tgau + 1)*(Tgau / 8 +
			sin(2 * pi*(Tgau + 1 / 2)) / (16 * pi) + 1 / 8);


	IdSwatxT =

		heaviside(Tgau + 1)*((sin(pi*Tgau) - pi * (Tgau*cos(pi*Tgau) - 1)) / (4 * pi) + (sin(2 * pi*Tgau) / 8 - pi * ((Tgau*cos(2 *
			pi*Tgau)) / 4 + 1 / 4)) / (4 * pi)) - heaviside(Tgau + 1)*((sin(pi*Tgau) - pi * (Tgau*cos(pi*Tgau) - 1)) / (4 * pi) -
			(sin(2 * pi*Tgau) / 8 - pi * ((Tgau*cos(2 * pi*Tgau)) / 4 + 1 / 4)) / (4 * pi)) + heaviside(Tgau)*((sin(pi*Tgau) - Tgau
				* pi*cos(pi*Tgau)) / (4 * pi) - (sin(2 * pi*Tgau) / 2 - Tgau * pi*cos(2 * pi*Tgau)) / (16 * pi)) - heaviside(Tgau)*
				((sin(pi*Tgau) - Tgau * pi*cos(pi*Tgau)) / (4 * pi) + (sin(2 * pi*Tgau) / 2 - Tgau * pi*cos(2 * pi*Tgau)) / (16 * pi));


	ISwatxT =

		heaviside(Tgau + 1)*((sin(pi*Tgau)*sin(pi*Tgau) / 16 - (Tgau*pi*sin(2 * pi*Tgau)) / 16) / (pi*pi) + Tgau * Tgau / 16 - 1 / 16) - Sres / 2 +
		(Sres*Tgau*Tgau) / 2 - heaviside(Tgau)*((sin(pi*Tgau)*sin(pi*Tgau) / 4 - (Tgau*pi*sin(2 * pi*Tgau)) / 4) / (4 * pi*pi) + Tgau * Tgau / 16);


	Sice =

		-heaviside(-Tgau - 1)*(Sres - 1) - (sin(pi*(Tgau + 1 / 2)) / 2 - 1 / 2)*(heaviside(Tgau + 1) - heaviside(Tgau))*(sin(pi*(Tgau + 1 / 2))*
		(Sres / 2 - 1 / 2) - Sres / 2 + 1 / 2);


	dSice =

		-(pi*cos(pi*(Tgau + 1 / 2))*(heaviside(Tgau + 1) - heaviside(Tgau))*(sin(pi*(Tgau + 1 / 2))*(Sres / 2 - 1 / 2) - Sres / 2 + 1 / 2)) / 2
		- pi * cos(pi*(Tgau + 1 / 2))*(Sres / 2 - 1 / 2)*(sin(pi*(Tgau + 1 / 2)) / 2 - 1 / 2)*(heaviside(Tgau + 1) - heaviside(Tgau));


	ISice =

		heaviside(Tgau)*((3 * Tgau*(Sres - 1)) / 8 - ((8 * sin(pi*Tgau) + sin(2 * pi*(Tgau + 1 / 2)))*(Sres - 1)) / (16 * pi)) -
		heaviside(Tgau + 1)*(((6 * Tgau + 6)*(Sres - 1)) / 16 - ((8 * sin(pi*Tgau) + sin(2 * pi*(Tgau + 1 / 2)))*(Sres - 1)) / (16 * pi));


	IdSicexT =

		(heaviside(Tgau)*(Sres - 1)*(8 * sin(pi*Tgau) - sin(2 * pi*Tgau) - 8 * Tgau*pi*cos(pi*Tgau) + 2 * Tgau*pi*cos(2 * pi*Tgau))) / (16 * pi)
		- (heaviside(Tgau + 1)*(Sres - 1)*(10 * pi + 8 * sin(pi*Tgau) - sin(2 * pi*Tgau) - 8 * Tgau*pi*cos(pi*Tgau) + 2 * Tgau*pi*cos(2 * pi*Tgau))) / (16 * pi);


	ISicexT =

		(heaviside(Tgau + 1)*((17 * Sres) / 2 + (Sres*(2 * sin(pi*Tgau)*sin(pi*Tgau) - 1)) / 2 - 8 * Sres*(2 * sin((pi*Tgau) / 2)*sin((pi*Tgau) / 2) - 1) + 3 * Sres*pi*pi - 3 * pi*pi + 3 *
			Tgau*Tgau * pi*pi - sin(pi*Tgau)*sin(pi*Tgau) + 16 * sin((pi*Tgau) / 2)*sin((pi*Tgau) / 2) - 3 * Sres*Tgau*Tgau * pi*pi - 8 * Tgau*pi*sin(pi*Tgau) + Tgau * pi*sin(2 * pi*Tgau) +
			8 * Sres*Tgau*pi*sin(pi*Tgau) - Sres * Tgau*pi*sin(2 * pi*Tgau) - 16)) / (16 * pi*pi) - heaviside(Tgau)*((((sin(pi*Tgau)*sin(pi*Tgau) - 16 * sin((pi*Tgau) / 2)*sin((pi*Tgau) / 2))
				*(Sres - 1)) / 16 + (pi*(Sres - 1)*(8 * Tgau*sin(pi*Tgau) - Tgau * sin(2 * pi*Tgau))) / 16) / pi * pi - (3 * Tgau*Tgau * (Sres - 1)) / 16);

}

//double SATUR(double Tgau, double Tsol, double Tliq, double Sres)
//{
//	double Sw;
//	double pi = PI;
//
//	Sw = Sres - (sin((PI*(Tliq - 2 * Tgau + Tsol)) / (2 * (Tliq - Tsol))) / 2 - 1 / 2)*(sin((PI*(Tliq - 2 * Tgau + Tsol))
//			/ (2 * (Tliq - Tsol))) / 2 + 1 / 2)*(heaviside(Tgau + 1) - heaviside(Tgau));
//
//	return Sw;
//}
//
//double dSATUR(double Tgau, double Tsol, double Tliq, double Sres)
//{
//	double dSw;
//	double pi = PI;
//
//	dSw = +(pi*sin((pi*(Tliq - 2 * Tgau + Tsol)) / (Tliq - Tsol))*(heaviside(Tgau + 1) - heaviside(Tgau))) / (4 * (Tliq - Tsol));
//
//	return dSw;
//}
//
//double ISATUR(double Tgau, double Tsol, double Tliq, double Sres)
//{
//	double ISw;
//	double pi = PI;
//
//	ISw = Sres * Tgau + heaviside(Tgau + 1)*((sin((pi*(Tliq + Tsol + 2)) / (Tliq - Tsol))*(Tliq - Tsol)) / (16 * pi) + 1 / 8)
//			- heaviside(Tgau)*(Tgau / 8 - (sin((pi*(Tliq - 2 * Tgau + Tsol)) / (Tliq - Tsol))*(Tliq - Tsol)) / (16 * pi)) +
//			heaviside(Tgau + 1)*(Tgau / 8 - (sin((pi*(Tliq - 2 * Tgau + Tsol)) / (Tliq - Tsol))*(Tliq - Tsol)) / (16 * pi)) -
//			(heaviside(Tgau)*sin((pi*(Tliq + Tsol)) / (Tliq - Tsol))*(Tliq - Tsol)) / (16 * pi);
//
//	return ISw;
//}
//
//double ISATURxT(double Tgau, double Tsol, double Tliq, double Sres)
//{
//	double ISwxT;
//	double pi = PI;
//
//	ISwxT = (Sres*Tgau*Tgau) / 2 - heaviside(Tgau)*(Tgau*Tgau / 16 + (cos((pi*(Tliq - 2 * Tgau + Tsol)) / (Tliq - Tsol))*
//				(Tliq * Tliq / 4 - (Tliq*Tsol) / 2 + Tsol * Tsol / 4)) / (8 * pi * pi) - (Tgau*sin((pi*(Tliq - 2 * Tgau + Tsol)) /
//				(Tliq - Tsol))*(Tliq - Tsol)) / (16 * pi)) + heaviside(Tgau + 1)*(Tgau*Tgau / 16 + (cos((pi*(Tliq - 2 * Tgau + Tsol))
//				/ (Tliq - Tsol))*(Tliq*Tliq / 4 - (Tliq*Tsol) / 2 + Tsol * Tsol / 4)) / (8 * pi*pi) - (Tgau*sin((pi*(Tliq - 2 * Tgau + Tsol))
//				/ (Tliq - Tsol))*(Tliq - Tsol)) / (16 * pi)) - heaviside(Tgau + 1)*((sin((pi*(Tliq + Tsol + 2)) / (Tliq - Tsol))*
//				(Tliq - Tsol)) / (16 * pi) + (cos((pi*(Tliq + Tsol + 2)) / (Tliq - Tsol))*(Tliq*Tliq / 4 - (Tliq*Tsol) / 2 + Tsol * Tsol / 4))
//				/ (8 * pi*pi) + 1 / 16) + (heaviside(Tgau)*cos((pi*(Tliq + Tsol)) / (Tliq - Tsol))*(Tliq*Tliq / 4 - (Tliq*Tsol) / 2 + Tsol * Tsol / 4)) / (8 * pi*pi);
//
//	return ISwxT;
//}
//
//double IdSATURxT(double Tgau, double Tsol, double Tliq, double Sres)
//{
//	double IdSwxT;
//	double pi = PI;
//
//	IdSwxT = -(2 * pi*heaviside(Tgau + 1) + Tliq * sin((pi*(Tliq - 2 * Tgau + Tsol)) / (Tliq - Tsol))*heaviside(Tgau) -
//				Tsol * sin((pi*(Tliq - 2 * Tgau + Tsol)) / (Tliq - Tsol))*heaviside(Tgau) - Tliq * heaviside(Tgau)*sin((pi*(Tliq + Tsol))
//				/ (Tliq - Tsol)) + Tsol * heaviside(Tgau)*sin((pi*(Tliq + Tsol)) / (Tliq - Tsol)) - Tliq * heaviside(Tgau + 1)*sin((pi*
//				(Tliq - 2 * Tgau + Tsol)) / (Tliq - Tsol)) + Tsol * heaviside(Tgau + 1)*sin((pi*(Tliq - 2 * Tgau + Tsol)) / (Tliq - Tsol)) +
//				Tliq * heaviside(Tgau + 1)*sin((pi*(Tliq + Tsol + 2)) / (Tliq - Tsol)) - Tsol * heaviside(Tgau + 1)*sin((pi*(Tliq + Tsol + 2)) / (Tliq - Tsol))
//				+ 2 * Tgau*pi*cos((pi*(Tliq - 2 * Tgau + Tsol)) / (Tliq - Tsol))*heaviside(Tgau) - 2 * Tgau*pi*heaviside(Tgau + 1)*cos((pi*
//				(Tliq - 2 * Tgau + Tsol)) / (Tliq - Tsol))) / (16 * pi);
//
//	return IdSwxT;
//}
//
//double SATURice(double Tgau, double Tsol, double Tliq, double Sres)
//{
//	double Si;
//	double pi = PI;
//
//	Si = (heaviside(Tgau + 1) - 1)*(Sres - 1) - (sin((pi*(Tliq - 2 * Tgau + Tsol)) / (2 * (Tliq - Tsol))) + 1)*
//			(Sres / 2 - 1 / 2)*(sin((pi*(Tliq - 2 * Tgau + Tsol)) / (2 * (Tliq - Tsol))) / 2 + 1 / 2)*(heaviside(Tgau + 1) - heaviside(Tgau));
//
//	return Si;
//}
//
//double dSATURice(double Tgau, double Tsol, double Tliq, double Sres)
//{
//	double dSi;
//	double pi = PI;
//
//	dSi = -(pi*(heaviside(Tgau + 1) - heaviside(Tgau))*(sin(pi*Tgau) - sin(2 * pi*Tgau) / 2)*(Sres - 1)) / 2;
//
//	return dSi;
//}
//
//double ISATURice(double Tgau, double Tsol, double Tliq, double Sres)
//{
//	double ISi;
//	double pi = PI;
//
//	ISi = heaviside(Tgau + 1)*(((Sres / 8 - 1 / 8)*(Tliq - Tsol)*(4 * cos((pi*(Tliq + Tsol + 2)) / (2 * (Tliq - Tsol))) + sin((pi*(Tliq + Tsol + 2))
//			/ (Tliq - Tsol)) / 2)) / pi - (3 * Sres) / 8 + 3 / 8) + heaviside(Tgau)*(Tgau*((3 * Sres) / 8 - 3 / 8) + ((Sres / 8 - 1 / 8)*(4 * cos((pi*(Tliq - 2 * Tgau + Tsol))
//			/ (2 * (Tliq - Tsol))) + sin((pi*(Tliq - 2 * Tgau + Tsol)) / (Tliq - Tsol)) / 2)*(Tliq - Tsol)) / pi) + heaviside(Tgau + 1)*(Sres - 1) -
//			heaviside(Tgau + 1)*(Tgau*((3 * Sres) / 8 - 3 / 8) + ((Sres / 8 - 1 / 8)*(4 * cos((pi*(Tliq - 2 * Tgau + Tsol)) / (2 * (Tliq - Tsol))) +
//			sin((pi*(Tliq - 2 * Tgau + Tsol)) / (Tliq - Tsol)) / 2)*(Tliq - Tsol)) / pi) + Tgau * (heaviside(Tgau + 1) - 1)*(Sres - 1) -
//			(heaviside(Tgau)*(Sres / 8 - 1 / 8)*(Tliq - Tsol)*(4 * cos((pi*(Tliq + Tsol)) / (2 * (Tliq - Tsol))) + sin((pi*(Tliq + Tsol)) / (Tliq - Tsol)) / 2)) / pi;
//
//	return ISi;
//}
//
//double IdSATURxTice(double Tgau, double Tsol, double Tliq, double Sres)
//{
//	double IdSixT;
//	double pi = PI;
//
//	IdSixT = ((sign(Tgau) / 2 + 1 / 2)*(Sres - 1)*(8 * sin(pi*Tgau) - sin(2 * pi*Tgau) - 6 * pi*Tgau + 16 *
//				Tgau*pi*sin((pi*Tgau) / 2)*sin((pi*Tgau) / 2)*sin((pi*Tgau) / 2)*sin((pi*Tgau) / 2))) / (32 * pi) - (pi*(sign(Tgau + 1) / 2 + 1 / 2)*(((8 * sin(pi*Tgau)
//				- 2 * cos(pi*Tgau)*sin(pi*Tgau))*(Sres - 1)) / (16 * pi*pi) - (Tgau*(Sres - 1)*(4 * cos(pi*Tgau)
//				- 2 * cos(pi*Tgau)*cos(pi*Tgau) + 1)) / (8 * pi))) / 2 - ((5 * Sres) / 16 - 5 / 16)*(sign(Tgau + 1) / 2 + 1 / 2) -
//				((Sres / 2 - 1 / 2)*(sign(Tgau) / 2 + 1 / 2)*((pi*Tgau) / 8 - sin(pi*Tgau) / 2 + (cos(pi*Tgau)*sin(pi*Tgau))
//				/ 8 + (Tgau*pi*cos(pi*Tgau)) / 2 - (Tgau*pi*cos(pi*Tgau)*cos(pi*Tgau)) / 4)) / pi - ((Sres / 2 - 1 / 2)*(sign(Tgau + 1) +
//				1)*(5 * pi + 4 * sin(pi*Tgau) - sin(2 * pi*Tgau) / 2 - pi * Tgau - 4 * Tgau*pi*cos(pi*Tgau) + 2 *
//				Tgau*pi*cos(pi*Tgau)*cos(pi*Tgau))) / (16 * pi);
//
//	return IdSixT;
//}
//
//double ISATURxTice(double Tgau, double Tsol, double Tliq, double Sres)
//{
//	double ISixT;
//	double pi = PI;
//
//	ISixT = heaviside(Tgau + 1)*(Sres / 8 + cos((pi*(Tliq + Tsol + 2)) / (2 * (Tliq - Tsol)))*cos((pi*(Tliq + Tsol + 2)) / (2 * (Tliq - Tsol)))
//				* (Sres / 16 - 1 / 16) + sin((pi*(Tliq + Tsol + 2)) / (2 * (Tliq - Tsol)))*sin((pi*(Tliq + Tsol + 2)) / (2 * (Tliq - Tsol))) * (Sres /
//				16 - 1 / 16) - (cos((pi*(Tliq + Tsol + 2)) / (2 * (Tliq - Tsol)))*(Tliq - Tsol)*(Sres - 1)) / (2 * pi) - (sin((pi*(Tliq + Tsol + 2))
//				/ (Tliq - Tsol))*(Tliq - Tsol)*(Sres - 1)) / (16 * pi) + (sin((pi*(Tliq + Tsol + 2)) / (2 * (Tliq - Tsol)))*(Tliq - Tsol)*(Tliq - Tsol)
//				* (Sres - 1)) / (2 * pi*pi) + (sin((pi*(Tliq + Tsol + 2)) / (2 * (Tliq - Tsol)))*sin((pi*(Tliq + Tsol + 2)) / (2 * (Tliq - Tsol)))
//				* (Tliq - Tsol)*(Tliq - Tsol) * (Sres - 1)) / (16 * pi*pi) - 1 / 8) - heaviside(Tgau + 1)*(Tgau*Tgau * (Sres / 8 - 1 / 8)
//				+ (Tgau*Tgau * cos((pi*(Tliq - 2 * Tgau + Tsol)) / (2 * (Tliq - Tsol)))*cos((pi*(Tliq - 2 * Tgau + Tsol)) / (2 * (Tliq - Tsol))) *
//				(Sres / 4 - 1 / 4)) / 4 + (Tgau*Tgau * sin((pi*(Tliq - 2 * Tgau + Tsol)) / (2 * (Tliq - Tsol)))*sin((pi*(Tliq - 2 * Tgau + Tsol))
//				/ (2 * (Tliq - Tsol))) * (Sres / 4 - 1 / 4)) / 4 + (sin((pi*(Tliq - 2 * Tgau + Tsol)) / (2 * (Tliq - Tsol)))*(Tliq - Tsol)*
//				(Tliq - Tsol) * (Sres - 1)) / (2 * pi*pi) + (sin((pi*(Tliq - 2 * Tgau + Tsol)) / (2 * (Tliq - Tsol)))*sin((pi*(Tliq - 2 * Tgau + Tsol))
//				/ (2 * (Tliq - Tsol))) * (Tliq - Tsol)*(Tliq - Tsol) * (Sres - 1)) / (16 * pi*pi) + (Tgau*cos((pi*(Tliq - 2 * Tgau + Tsol)) /
//				(2 * (Tliq - Tsol)))*(Tliq - Tsol)*(Sres - 1)) / (2 * pi) + (Tgau*sin((pi*(Tliq - 2 * Tgau + Tsol)) / (Tliq - Tsol))*
//				(Tliq - Tsol)*(Sres - 1)) / (16 * pi)) - heaviside(Tgau + 1)*(Sres / 2 - 1 / 2) + heaviside(Tgau)*(Tgau*Tgau *
//				(Sres / 8 - 1 / 8) + (Tgau*Tgau * cos((pi*(Tliq - 2 * Tgau + Tsol)) / (2 * (Tliq - Tsol)))*cos((pi*(Tliq - 2 *
//				Tgau + Tsol)) / (2 * (Tliq - Tsol))) * (Sres / 4 - 1 / 4)) / 4 + (Tgau*Tgau * sin((pi*(Tliq - 2 * Tgau + Tsol))
//				/ (2 * (Tliq - Tsol)))*sin((pi*(Tliq - 2 * Tgau + Tsol)) / (2 * (Tliq - Tsol))) * (Sres / 4 - 1 / 4)) / 4 +
//				(sin((pi*(Tliq - 2 * Tgau + Tsol)) / (2 * (Tliq - Tsol)))*(Tliq - Tsol)*(Tliq - Tsol) * (Sres - 1)) / (2 * pi*pi)
//				+ (sin((pi*(Tliq - 2 * Tgau + Tsol)) / (2 * (Tliq - Tsol))) * sin((pi*(Tliq - 2 * Tgau + Tsol)) / (2 * (Tliq - Tsol)))
//				* (Tliq - Tsol)*(Tliq - Tsol) * (Sres - 1)) / (16 * pi*pi) + (Tgau*cos((pi*(Tliq - 2 * Tgau + Tsol)) /
//				(2 * (Tliq - Tsol)))*(Tliq - Tsol)*(Sres - 1)) / (2 * pi) + (Tgau*sin((pi*(Tliq - 2 * Tgau + Tsol)) /
//				(Tliq - Tsol))*(Tliq - Tsol)*(Sres - 1)) / (16 * pi)) + Tgau * Tgau * (Sres / 2 - 1 / 2)*
//				(heaviside(Tgau + 1) - 1) - (heaviside(Tgau)*sin((pi*(Tliq + Tsol)) / (2 * (Tliq - Tsol)))*(sin((pi*(Tliq + Tsol))
//				/ (2 * (Tliq - Tsol))) + 8)*(Tliq - Tsol)*(Tliq - Tsol) * (Sres - 1)) / (16 * pi*pi);
//
//	return ISixT;
//}

//
//double SATUR(double Tgau, double Tsol, double Tliq, double Sres, double Wpar, double Mpar, int rSFC)
//{
//	double satur, Btem;
//	switch (rSFC)
//	{
//	case 1:
//		if (Tgau <= Tsol)
//		{
//			satur = Sres;
//		}
//		else if (Tgau < Tliq && Tgau > Tsol)
//		{
//			satur = (1.0 - Sres)*exp(-((Tgau - Tliq) / Wpar)*((Tgau - Tliq) / Wpar)) + Sres;
//		}
//		else
//		{
//			satur = 1.0;
//
//		}
//
//		break;
//
//	case 0:
//		Btem = (Sres - 1) / Mpar;
//		if (Tgau >= Btem)
//		{
//			satur = Mpar * Tgau + 1.0;
//			if (Tgau > Tliq)
//			{
//				satur = 1.0;
//			}
//		}
//		else
//		{
//			satur = Sres;
//		}
//
//		break;
//
//	case 2:
//		if (Tgau <= Tsol)
//		{
//			satur = Sres;
//		}
//		else if (Tgau <Tliq && Tgau > Tsol)
//		{
//			satur = 0.5*(1.0 + Sres) + 0.5*(1.0 - Sres)*sin(PI*(Tgau - 0.5*(Tliq + Tsol)) / (Tliq - Tsol));
//		}
//		else
//		{
//			satur = 1.0;
//		}
//
//		break;
//	}
//
//	return satur;
//}
//
//double dSATUR(double Tgau, double Tsol, double Tliq, double Sres, double Wpar, double Mpar, int rSFC)
//{
//	double dSatur, Btem;
//	//
//	rSFC = PROPS.Soil.rSFC;
//	Tsol = PROPS.Nonisothermal.TempSolid;
//	Tliq = PROPS.Nonisothermal.TempLiquid;
//	Sres = PROPS.Soil.ResidualWaterSaturation;
//	Wpar = PROPS.Soil.Wpar;
//	Mpar = PROPS.Soil.Mpar;
//	//
//	switch (rSFC)
//	{
//	case 1:
//		if (Tgau <= Tsol)
//		{
//			dSatur = 0.0;
//		}
//		else if (Tgau < Tliq && Tgau > Tsol)
//		{
//			dSatur = (exp(-(Tgau - Tliq) * (Tgau - Tliq) / Wpar * Wpar)*(Sres - 1.0)*(2.0 * Tgau - 2.0 * Tliq)) / (Wpar * Wpar);
//		}
//		else
//		{
//			dSatur = 0;
//		}
//
//		break;
//
//	case 0:
//		Btem = (Sres - 1.0) / Mpar;
//		if (Tgau >= Btem)
//		{
//			dSatur = Mpar;
//			if (Tgau > Tliq)
//			{
//				dSatur = 0.0;
//			}
//		}
//		else
//		{
//			dSatur = 0.0;
//		}
//
//		break;
//
//	case 2:
//		if (Tgau <= Tsol)
//		{
//			dSatur = Sres;
//		}
//		else if (Tgau < Tliq && Tgau > Tsol)
//		{
//			dSatur = -(PI*cos((PI*(Tliq / 2.0 - Tgau + Tsol / 2.0)) / (Tliq - Tsol))*(Sres / 2.0 - 1.0 / 2.0)) / (Tliq - Tsol);
//		}
//		else
//		{
//			dSatur = 0.0;
//		}
//
//		break;
//	}
//
//	return dSatur;
//}
//
//double ISATUR(double Tgau, double Tsol, double Tliq, double Sres, double Wpar, double Mpar, int rSFC)
//{
//	double iSatur, CSw1, CSw2, CSw1p, Btem;
//	//
//	rSFC = PROPS.Soil.rSFC;
//	Tsol = PROPS.Nonisothermal.TempSolid;
//	Tliq = PROPS.Nonisothermal.TempLiquid;
//	Sres = PROPS.Soil.ResidualWaterSaturation;
//	Wpar = PROPS.Soil.Wpar;
//	Mpar = PROPS.Soil.Mpar;
//	//
//	switch (rSFC)
//	{
//	case 1:
//		if (Tgau <= Tsol)
//		{
//			iSatur = Sres * Tgau;
//		}
//		else if (Tgau < Tliq && Tgau > Tsol)
//		{
//			/*CSw1 = Sres * Tsol - (Sres * Tsol - (sqrt(PI)*erf((Tsol - Tliq)*(-1.0 / (Wpar * Wpar)) ^ (1 / 2))*(Sres - 1)) / (2.0 * (-1.0 / (Wpar * Wpar)) ,(1.0 / 2.0)));
//			iSatur = CSw1 + Sres * Tgau - (sqrt(PI)*erf((Tgau - Tliq)*(-1 / (Wpar * Wpar)) ^ (1 / 2))*(Sres - 1.0)) / (2 * sqrt((-1.0 / Wpar * Wpar)));*/
//		}
//
//		break;
//
//	case 0:
//		Btem = (Sres - 1.0) / Mpar;
//		if (Tgau >= Btem)
//		{
//			CSw1 = Sres * Btem - (Mpar*Btem*Btem / 2 + Btem);
//			iSatur = CSw1 + Mpar * Tgau * Tgau / 2.0 + Tgau;
//			if (Tgau > Tliq)
//			{
//				CSw2 = (CSw1 + Mpar * Tliq * Tliq / 2.0 + Tliq) - Tliq;
//				iSatur = CSw2 + Tgau;
//			}
//		}
//		else
//		{
//			iSatur = Sres * Tgau;
//		}
//
//		break;
//
//	case 2:
//		if (Tgau <= Tsol)
//		{
//			iSatur = Sres * Tgau;
//		}
//		else if (Tgau < Tliq && Tgau > Tsol)
//		{
//			CSw1 = Sres * Tsol - ((Tsol*(Sres + 1.0)) / 2.0 + ((Tliq - Tsol)*(Sres - 1.0)) / (2.0 * PI));
//			iSatur = CSw1 + (Tgau*(Sres + 1.0)) / 2.0 + (cos((PI*(Tliq - 2.0 * Tgau + Tsol)) / (4 * (Tliq - Tsol))) * cos((PI*(Tliq - 2.0 * Tgau + Tsol)) / (4 * (Tliq - Tsol))) * (Tliq - Tsol)*(Sres - 1.0)) / PI;
//		}
//		else
//		{
//			CSw1p = Sres * Tsol - ((Tsol*(Sres + 1.0)) / 2.0 + ((Tliq - Tsol)*(Sres - 1.0)) / (2.0 * PI));
//			CSw1 = CSw1p + (Tliq*(Sres + 1)) / 2 + ((Tliq - Tsol)*(Sres - 1)) / (2.0 * PI);
//			CSw2 = CSw1 - (Tliq);
//			iSatur = CSw2 + Tgau;
//		}
//
//		break;
//	}
//
//	return iSatur;
//}
//
//double ISATURxT(double Tgau, double Tsol, double Tliq, double Sres, double Wpar, double Mpar, int rSFC)
//{
//	double iSatxT, CSw1, CSw2, CSw1p, Btem;
//	//
//	rSFC = PROPS.Soil.rSFC;
//	Tsol = PROPS.Nonisothermal.TempSolid;
//	Tliq = PROPS.Nonisothermal.TempLiquid;
//	Sres = PROPS.Soil.ResidualWaterSaturation;
//	Wpar = PROPS.Soil.Wpar;
//	Mpar = PROPS.Soil.Mpar;
//	//
//	switch (rSFC)
//	{
//	case 0:
//		Btem = (Sres - 1.0) / Mpar;
//		if (Tgau >= Btem)
//		{
//			CSw1 = Sres * Btem * Btem / 2.0 - (Mpar*Btem*Btem*Btem / 3.0 + Btem * Btem / 2.0);
//			iSatxT = CSw1 + Mpar * Tgau * Tgau * Tgau / 3.0 + Tgau * Tgau / 2.0;
//			if (Tgau > Tliq)
//			{
//				CSw2 = (CSw1 + Mpar * Tliq * Tliq * Tliq / 3.0 + Tliq * Tliq / 2.0) - Tliq * Tliq / 2.0;
//				iSatxT = CSw2 + Tgau * Tgau / 2;
//			}
//		}
//		else
//		{
//			iSatxT = Sres * Tgau * Tgau / 2;
//		}
//
//		break;
//
//	case 2:
//		if (Tgau <= Tsol)
//			iSatxT = Sres * Tgau * Tgau / 2;
//		else if (Tgau < Tliq && Tgau > Tsol)
//		{
//			CSw1 = Sres * Tsol * Tsol / 2.0 - (Tsol * Tsol * (Sres / 4.0 + 1.0 / 4.0) + ((Tliq - Tsol) * (Tliq - Tsol) * (Sres - 1.0)) / (2.0 * PI * PI));
//			iSatxT = CSw1 + Tgau * Tgau * (Sres / 4.0 + 1.0 / 4.0) + (sin((PI*(Tliq - 2.0 * Tgau + Tsol)) / (2.0 * (Tliq
//				- Tsol)))*(Tliq - Tsol) * (Tliq - Tsol) * (Sres - 1.0)) / (2.0 * PI * PI) + (Tgau*cos((PI*(Tliq
//					- 2.0 * Tgau + Tsol)) / (2.0 * (Tliq - Tsol)))*(Tliq - Tsol)*(Sres - 1.0)) / (2.0 * PI);
//		}
//		else
//		{
//			CSw1p = Sres * Tsol * Tsol / 2.0 - (Tsol * Tsol * (Sres / 4.0 + 1.0 / 4.0) + ((Tliq - Tsol) * (Tliq - Tsol) * (Sres - 1.0)) / (2.0 * PI * PI));
//			CSw1 = CSw1p + Tliq * Tliq * (Sres / 4.0 + 1.0 / 4.0) - ((Tliq - Tsol) * (Tliq - Tsol) * (Sres - 1.0)) / (2.0 * PI * PI);
//			CSw2 = CSw1 - (Tliq * Tliq / 2.0);
//			iSatxT = CSw2 + Tgau * Tgau / 2.0;
//		}
//
//		break;
//	}
//
//	return iSatxT;
//}
//
//double IdSATURxT(double Tgau, double Tsol, double Tliq, double Sres, double Wpar, double Mpar, int rSFC)
//{
//	double idSaturxT, IdSwxT_0, IdSwxT_0p, Btem;
//	//
//	rSFC = PROPS.Soil.rSFC;
//	Tsol = PROPS.Nonisothermal.TempSolid;
//	Tliq = PROPS.Nonisothermal.TempLiquid;
//	Sres = PROPS.Soil.ResidualWaterSaturation;
//	Wpar = PROPS.Soil.Wpar;
//	Mpar = PROPS.Soil.Mpar;
//	//
//	switch (rSFC)
//	{
//	case 0:
//		Btem = (Sres - 1.0) / Mpar;
//		if (Tgau >= Btem)
//		{
//			IdSwxT_0 = Mpar * Btem * Btem / 2.0;
//			idSaturxT = Mpar * Tgau * Tgau / 2.0 - IdSwxT_0;
//			if (Tgau > Tliq)
//			{
//				IdSwxT_0p = Mpar * Btem * Btem / 2.0;
//				IdSwxT_0 = Mpar * Tliq * Tliq / 2.0 - IdSwxT_0p;
//				idSaturxT = IdSwxT_0;
//			}
//		}
//		else
//		{
//			idSaturxT = 0;
//		}
//
//		break;
//
//	case 2:
//		if (Tgau <= Tsol)
//			idSaturxT = Sres * Tgau * Tgau / 2;
//		else if (Tgau < Tliq && Tgau > Tsol)
//		{
//			IdSwxT_0 = -((cos((PI*(Tliq / 2.0 - Tsol / 2.0)) / (Tliq -
//				Tsol))*(Tliq - Tsol) * (Tliq - Tsol) * (Sres - 1)) / PI -
//				Tsol*sin((PI*(Tliq / 2.0 - Tsol / 2.0)) / (Tliq -
//					Tsol))*(Tliq - Tsol)*(Sres - 1)) / (2.0 * Tliq - 2.0 * Tsol);
//			idSaturxT = (Tgau*sin((PI*(Tliq / 2.0 - Tgau + Tsol / 2.0)) / (Tliq -
//				Tsol))*(Tliq - Tsol)*(Sres - 1.0) -
//				(cos((PI*(Tliq / 2 - Tgau + Tsol / 2)) / (Tliq - Tsol)) * (Tliq - Tsol)
//					* (Tliq - Tsol) * (Sres - 1.0)) / PI) / (2.0 * Tliq - 2.0 * Tsol) - IdSwxT_0;
//		}
//		else
//		{
//			IdSwxT_0p = -((cos((PI*(Tliq / 2.0 - Tsol / 2.0)) / (Tliq -
//				Tsol))*(Tliq - Tsol) * (Tliq - Tsol)  * (Sres - 1.0)) / PI -
//				Tsol*sin((PI*(Tliq / 2.0 - Tsol / 2.0)) / (Tliq -
//					Tsol))*(Tliq - Tsol)*(Sres - 1.0)) / (2.0 * Tliq - 2.0 * Tsol);
//			IdSwxT_0 = -((cos((PI*(Tliq / 2.0 - Tsol / 2.0)) / (Tliq -
//				Tsol))*(Tliq - Tsol) * (Tliq - Tsol) * (Sres - 1.0)) / PI +
//				Tliq*sin((PI*(Tliq / 2.0 - Tsol / 2.0)) / (Tliq -
//					Tsol))*(Tliq - Tsol)*(Sres - 1.0)) / (2.0 * Tliq - 2.0 * Tsol) - IdSwxT_0p;
//			idSaturxT = IdSwxT_0;
//		}
//
//		break;
//	}
//
//	return idSaturxT;
//}