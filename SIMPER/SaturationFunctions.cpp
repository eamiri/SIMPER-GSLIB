#include "SimperInclude.h"

SaturationFunctions::SaturationFunctions(double Tgau, double Tsol, double Tliq, double Sres, bool isSaturated)
{
	double pi = PI;

	if (isSaturated)
	{
		Swat =

			heaviside(Tgau - Tliq) + Sres * heaviside(Tsol - Tgau) - (heaviside(Tgau - Tliq) - heaviside(Tgau - Tsol))*
			(Sres / 2.0 + sin((pi*(Tliq / 2.0 - Tgau + Tsol / 2.0)) / (Tliq - Tsol))*(Sres / 2.0 - 1.0 / 2.0) + 1.0 / 2.0);


		dSwat =

			+(pi*cos((pi*(Tliq / 2.0 - Tgau + Tsol / 2.0)) / (Tliq - Tsol))*(Sres / 2.0 - 1.0 / 2.0)
				*(heaviside(Tgau - Tliq) - heaviside(Tgau - Tsol))) / (Tliq - Tsol);


		ISwat =

			Sres * Tgau - heaviside(Tgau - Tliq)*((Tgau*(Sres + 1.0)) / 2.0 +
			(cos((pi*(Tliq - 2.0 * Tgau + Tsol)) / (4.0 * (Tliq - Tsol)))*cos((pi*(Tliq - 2.0 * Tgau + Tsol)) / (4.0 * (Tliq - Tsol))) *
				(Tliq - Tsol)*(Sres - 1.0)) / pi) + heaviside(Tgau - Tsol)*
				((Tgau*(Sres + 1.0)) / 2.0 + (cos((pi*(Tliq - 2.0 * Tgau + Tsol))
					/ (4.0 * (Tliq - Tsol)))*cos((pi*(Tliq - 2.0 * Tgau + Tsol))
						/ (4.0 * (Tliq - Tsol))) * (Tliq - Tsol)*(Sres - 1.0)) / pi)
			- (Tsol*(Sres + 1.0)) / 2.0 + Tgau * heaviside(Tgau - Tliq) - Tliq *
			heaviside(Tliq - Tsol) - Tliq * (heaviside(Tgau - Tliq) - 1.0) + Tsol
			* (heaviside(Tliq - Tsol) - 1.0) + ((heaviside(Tgau - Tliq) - 1.0)
				*(Tsol - Tliq + pi * Tliq + Sres * Tliq - Sres * Tsol + pi * Sres*Tliq))
			/ (2.0 * pi) - ((heaviside(Tgau - Tsol) - 1.0)*(Tsol - Tliq + pi *
				Tsol + Sres * Tliq - Sres * Tsol + pi * Sres*Tsol)) / (2.0 * pi)
			- ((heaviside(Tliq - Tsol) - 1.0)*(Tsol - Tliq + pi * Tsol + Sres *
				Tliq - Sres * Tsol + pi * Sres*Tsol)) / (2.0 * pi) - ((Tliq - Tsol)
					*(Sres - 1.0)) / (2.0 * pi) + Sres * Tsol*(heaviside(Tgau - Tsol) - 1.0)
			+ (heaviside(Tliq - Tsol)*(Tsol - Tliq + pi *
				Tliq + Sres * Tliq - Sres * Tsol + pi * Sres*Tliq)) /
				(2.0 * pi) - Sres * Tgau*heaviside(Tgau - Tsol);


		ISwatxT =

			heaviside(Tliq - Tsol)*(Tliq*Tliq * (Sres / 4.0 + 1.0 / 4.0) - ((Tliq - Tsol)*(Tliq - Tsol) * (Sres - 1.0)) / (2.0 * pi*pi))
			- Tsol * Tsol * (Sres / 4.0 + 1.0 / 4.0) - (Tliq*Tliq * (heaviside(Tgau - Tliq) - 1.0)) / 2.0 +
			(Tsol*Tsol * (heaviside(Tliq - Tsol) - 1.0)) / 2.0 + (Sres*Tgau*Tgau) / 2.0 + (Tliq*Tliq *
			(Sres / 4.0 + 1.0 / 4.0) - ((Tliq - Tsol)*(Tliq - Tsol) * (Sres - 1.0)) / (2.0 * pi*pi))*(heaviside(Tgau - Tliq) - 1.0) -
				(Tsol*Tsol * (Sres / 4.0 + 1.0 / 4.0) + ((Tliq - Tsol)*(Tliq - Tsol) * (Sres - 1.0)) / (2.0 * pi*pi))*
			(heaviside(Tgau - Tsol) - 1.0) - (Tsol*Tsol * (Sres / 4.0 + 1.0 / 4.0) + ((Tliq - Tsol)*(Tliq - Tsol) *
			(Sres - 1.0)) / (2.0 * pi*pi))*(heaviside(Tliq - Tsol) - 1.0) + (Tgau*Tgau * heaviside(Tgau - Tliq))
			/ 2.0 - (Tliq*Tliq * heaviside(Tliq - Tsol)) / 2.0 -
			heaviside(Tgau - Tliq)*(Tgau*Tgau * (Sres / 4.0 + 1.0 / 4.0) +
			(sin((pi*(Tliq - 2.0 * Tgau + Tsol)) / (2.0 * (Tliq - Tsol)))*(Tliq - Tsol)*(Tliq - Tsol) * (Sres - 1.0)) /
				(2.0 * pi*pi) + (Tgau*cos((pi*(Tliq - 2.0 * Tgau + Tsol)) / (2.0 * (Tliq - Tsol)))*(Tliq - Tsol)
					*(Sres - 1.0)) / (2.0 * pi)) + heaviside(Tgau - Tsol)*(Tgau*Tgau * (Sres / 4.0 + 1.0 / 4.0) +
					(sin((pi*(Tliq - 2.0 * Tgau + Tsol)) / (2.0 * (Tliq - Tsol)))*(Tliq - Tsol)*(Tliq - Tsol) * (Sres - 1.0)) /
						(2.0 * pi*pi) + (Tgau*cos((pi*(Tliq - 2.0 * Tgau + Tsol)) / (2.0 * (Tliq - Tsol)))*
						(Tliq - Tsol)*(Sres - 1.0)) / (2.0 * pi)) - (Sres*Tgau*Tgau * heaviside(Tgau - Tsol))
			/ 2.0 - ((Tliq - Tsol)*(Tliq - Tsol) * (Sres - 1.0)) / (2.0 * pi*pi) + (Sres*Tsol*Tsol * (heaviside(Tgau - Tsol) - 1.0)) / 2.0;


		Sice =

			(heaviside(Tgau - Tsol) - 1.0)*(Sres - 1.0) + (sin((pi*(Tliq - 2.0 * Tgau + Tsol)) / (2.0 * (Tliq - Tsol))) + 1.0)*(Sres / 2.0 - 1.0 / 2.0)*(heaviside(Tgau - Tliq) - heaviside(Tgau - Tsol));


		dSice =

			-(pi*cos((pi*(Tliq - 2.0 * Tgau + Tsol)) / (2.0 * (Tliq - Tsol)))*(Sres / 2.0 - 1.0 / 2.0)*(heaviside(Tgau - Tliq) - heaviside(Tgau - Tsol))) / (Tliq - Tsol);


		ISice =

			((Sres - 1.0)*(pi*Tliq - 2.0 * pi*Tgau + pi * Tsol - Tliq * heaviside(Tgau - Tliq) + Tliq *
				heaviside(Tgau - Tsol) + Tsol * heaviside(Tgau - Tliq) - Tsol * heaviside(Tgau - Tsol)
				+ 2.0 * Tliq*cos((pi*(Tliq - 2.0 * Tgau + Tsol)) / (4.0 * (Tliq - Tsol))) * cos((pi*(Tliq - 2.0 * Tgau + Tsol)) / (4.0 * (Tliq - Tsol))) * heaviside(Tgau - Tliq)
				- 2.0 * Tliq*cos((pi*(Tliq - 2.0 * Tgau + Tsol)) / (4.0 * (Tliq - Tsol)))*cos((pi*(Tliq - 2.0 * Tgau + Tsol)) / (4.0 * (Tliq - Tsol))) * heaviside(Tgau - Tsol)
				- 2.0 * Tsol*cos((pi*(Tliq - 2.0 * Tgau + Tsol)) / (4.0 * (Tliq - Tsol))) * cos((pi*(Tliq - 2.0 * Tgau + Tsol)) / (4.0 * (Tliq - Tsol))) * heaviside(Tgau - Tliq)
				+ 2.0 * Tsol*cos((pi*(Tliq - 2.0 * Tgau + Tsol)) / (4.0 * (Tliq - Tsol)))*cos((pi*(Tliq - 2.0 * Tgau + Tsol)) / (4.0 * (Tliq - Tsol))) * heaviside(Tgau - Tsol)
				+ Tgau * pi*heaviside(Tgau - Tliq) - Tliq * pi*heaviside(Tgau - Tliq) + Tgau *
				pi*heaviside(Tgau - Tsol) - Tliq * pi*heaviside(Tliq - Tsol) - Tsol * pi*heaviside(Tgau - Tsol)
				+ Tsol * pi*heaviside(Tliq - Tsol))) / (2.0 * pi);


		IdSice =

			(Sres / 2.0 - 1.0 / 2.0)*(heaviside(Tgau - Tliq) + heaviside(Tgau - Tsol) + 2.0 * heaviside(Tliq - Tsol)
				+ sin((pi*(Tliq - 2.0 * Tgau + Tsol)) / (2.0 * (Tliq - Tsol)))*heaviside(Tgau - Tliq) -
				sin((pi*(Tliq - 2.0 * Tgau + Tsol)) / (2.0 * (Tliq - Tsol)))*heaviside(Tgau - Tsol) - 2.0);


		IdSicexT =

			((Sres - 1.0)*(Tliq*heaviside(Tliq - Tsol) - Tsol + Tsol * heaviside(Tliq - Tsol)))
			/ 2.0 + (Tsol*(Sres - 1.0)) / 2.0 + Tsol * (heaviside(Tgau - Tsol) - 1.0)*(Sres - 1.0) -
			2.0 * Tsol*(Sres / 2.0 - 1.0 / 2.0)*(heaviside(Tgau - Tsol) - 1.0) - (pi*(Sres / 2.0 - 1.0 / 2.0)*
			(heaviside(Tgau - Tliq)*((cos((pi*(Tliq - 2.0 * Tgau + Tsol)) / (2.0 * (Tliq - Tsol)))*
				(Tliq - Tsol)*(Tliq - Tsol)) / (pi * pi) - (Tgau*sin((pi*(Tliq - 2.0 * Tgau + Tsol)) /
				(2.0 * (Tliq - Tsol)))*(Tliq - Tsol)) / pi) - (Tliq*(heaviside(Tgau - Tliq) - 1.0)
					*(Tliq - Tsol)) / pi)) / (Tliq - Tsol) + (pi*(Sres / 2.0 - 1.0 / 2.0)*(heaviside(Tgau - Tsol)
						*((cos((pi*(Tliq - 2.0 * Tgau + Tsol)) / (2.0 * (Tliq - Tsol)))*(Tliq - Tsol)*(Tliq - Tsol))
							/ (pi *pi) - (Tgau*sin((pi*(Tliq - 2.0 * Tgau + Tsol)) / (2.0 * (Tliq - Tsol)))*
							(Tliq - Tsol)) / pi) + (Tsol*(heaviside(Tgau - Tsol) - 1.0)*(Tliq - Tsol)) / pi))
			/ (Tliq - Tsol);


		ISicexT =

			Tsol * Tsol * (Sres / 4.0 - 1.0 / 4.0) - Tgau * Tgau * (Sres / 2.0 - 1.0 / 2.0) + heaviside(Tliq - Tgau)
			*(Tliq *Tliq * (Sres / 4.0 - 1.0 / 4.0) - ((Tliq - Tsol)*(Tliq - Tsol) * (Sres - 1.0)) / (2.0 * pi*pi))
			- heaviside(Tliq - Tsol)*(Tliq *Tliq * (Sres / 4.0 - 1.0 / 4.0) - ((Tliq - Tsol)*(Tliq - Tsol) * (Sres - 1.0))
				/ (2.0 * pi *pi)) - heaviside(Tsol - Tgau)*(Tsol *Tsol * (Sres / 4.0 - 1.0 / 4.0) + ((Tliq - Tsol)*(Tliq - Tsol) *
				(Sres - 1.0)) / (2.0 * pi *pi)) - heaviside(Tsol - Tliq)*(Tsol *Tsol * (Sres / 4.0 - 1.0 / 4.0) +
					((Tliq - Tsol)*(Tliq - Tsol) * (Sres - 1.0)) / (2.0 * pi*pi)) + heaviside(Tgau - Tliq)*(Tgau *Tgau
						* (Sres / 4.0 - 1.0 / 4.0) + (sin((pi*(Tliq - 2.0 * Tgau + Tsol)) / (2.0 * (Tliq - Tsol)))
							*(Tliq - Tsol)*(Tliq - Tsol) * (Sres - 1.0)) / (2.0 * pi*pi) + (Tgau*cos((pi*(Tliq - 2.0 *
								Tgau + Tsol)) / (2.0 * (Tliq - Tsol)))*(Tliq - Tsol)*(Sres - 1.0)) / (2.0 * pi))
			- heaviside(Tgau - Tsol)*(Tgau*Tgau * (Sres / 4.0 - 1.0 / 4.0) + (sin((pi*(Tliq - 2.0 * Tgau + Tsol)) /
			(2.0 * (Tliq - Tsol)))*(Tliq - Tsol)*(Tliq - Tsol) * (Sres - 1.0)) / (2.0 * pi *pi) + (Tgau*cos((pi*(Tliq - 2.0
				* Tgau + Tsol)) / (2.0 * (Tliq - Tsol)))*(Tliq - Tsol)*(Sres - 1.0)) / (2.0 * pi)) + ((Tliq - Tsol)
					*(Tliq - Tsol) * (Sres - 1.0)) / (2.0 * pi *pi) + Tgau * Tgau * heaviside(Tgau - Tsol)*(Sres / 2.0 - 1.0 / 2.0) +
			Tsol * Tsol * heaviside(Tsol - Tgau)*(Sres / 2.0 - 1.0 / 2.0);


		Sair = 0.0;
		ISair = 0.0;
		ISairxT = 0.0;

		Sair = 0.0;
		ISair = 0.0;
		ISairxT = 0.0;
	}
	else
	{
		Swat =

			Sres + (sin((pi*(Tliq / 2.0 - Tgau + Tsol / 2.0)) / (Tliq - Tsol)) / 2.0 - 1.0 / 2.0)*
			(sin((pi*(Tliq / 2.0 - Tgau + Tsol / 2.0)) / (Tliq - Tsol)) / 2.0 + 1.0 / 2.0)*(heaviside(Tgau - Tliq) - heaviside(Tgau - Tsol));


		dSwat =

			-(pi*sin((pi*(Tliq - 2.0 * Tgau + Tsol)) / (Tliq - Tsol))*(heaviside(Tgau - Tliq) - heaviside(Tgau - Tsol))) / (4.0 * (Tliq - Tsol));


		ISwat =

			Sres * Tgau - Tsol / 8.0 - Sres * Tsol - heaviside(Tgau - Tliq)*(Tgau / 8.0 - (sin((pi*(Tliq - 2.0 * Tgau + Tsol))
				/ (Tliq - Tsol))*(Tliq - Tsol)) / (16.0 * pi)) + heaviside(Tgau - Tsol)*(Tgau / 8.0 - (sin((pi*(Tliq - 2.0 * Tgau + Tsol))
					/ (Tliq - Tsol))*(Tliq - Tsol)) / (16.0 * pi)) + (Tliq*heaviside(Tliq - Tsol)) / 8.0 + (Tliq*(heaviside(Tgau - Tliq) - 1.0))
			/ 8.0 - (Tsol*(heaviside(Tgau - Tsol) - 1.0)) / 8.0 - (Tsol*(heaviside(Tliq - Tsol) - 1.0)) / 8.0;


		ISwatxT =

			(Sres*Tgau*Tgau) / 2.0 - (Sres*Tsol*Tsol) / 2.0 - heaviside(Tgau - Tliq)*(Tgau*Tgau / 16.0 + (cos((pi*(Tliq - 2.0 * Tgau + Tsol)) /
			(Tliq - Tsol))*(Tliq - Tsol)*(Tliq - Tsol)) / (32.0 * pi*pi) - (Tgau*sin((pi*(Tliq - 2.0 * Tgau + Tsol)) / (Tliq - Tsol))*(Tliq - Tsol))
				/ (16.0 * pi)) + heaviside(Tgau - Tsol)*(Tgau*Tgau / 16.0 + (cos((pi*(Tliq - 2.0 * Tgau + Tsol)) / (Tliq - Tsol))*(Tliq - Tsol)*(Tliq - Tsol))
					/ (32.0 * pi*pi) - (Tgau*sin((pi*(Tliq - 2.0 * Tgau + Tsol)) / (Tliq - Tsol))*(Tliq - Tsol)) / (16.0 * pi)) - Tsol * Tsol / 16.0
			+ (Tliq - Tsol)*(Tliq - Tsol) / (32.0 * pi*pi) + (heaviside(Tliq - Tsol)*(2.0 * Tliq*Tsol - Tliq * Tliq - Tsol * Tsol + 2.0 * Tliq*Tliq * pi*pi))
			/ (32.0 * pi*pi) + ((heaviside(Tgau - Tliq) - 1.0)*(2.0 * Tliq*Tsol - Tliq * Tliq - Tsol * Tsol + 2.0 * Tliq*Tliq * pi*pi)) / (32.0 * pi*pi)
			- ((heaviside(Tgau - Tsol) - 1.0)*(2.0 * Tliq*Tsol - Tliq * Tliq - Tsol * Tsol + 2.0 * Tsol*Tsol * pi*pi)) / (32.0 * pi*pi) -
			((heaviside(Tliq - Tsol) - 1.0)*(2.0 * Tliq*Tsol - Tliq * Tliq - Tsol * Tsol + 2.0 * Tsol*Tsol * pi*pi)) / (32.0 * pi*pi);

		Sice =

			(sin((pi*(Tliq / 2.0 - Tgau + Tsol / 2.0)) / (Tliq - Tsol)) / 2.0 + 1.0 / 2.0)*(heaviside(Tgau - Tliq)
				- heaviside(Tgau - Tsol))*(Sres / 2.0 + sin((pi*(Tliq / 2.0 - Tgau + Tsol / 2.0)) /
				(Tliq - Tsol))*(Sres / 2.0 - 1.0 / 2.0) - 1.0 / 2.0) - heaviside(Tsol - Tgau)*(Sres - 1.0);


		dSice =

			-(pi*cos((pi*(Tliq - 2.0 * Tgau + Tsol)) / (2.0 * (Tliq - Tsol)))*(sin((pi*(Tliq - 2.0 * Tgau + Tsol)) / (2.0 * (Tliq - Tsol))) + 1.0)*
			(heaviside(Tgau - Tliq) - heaviside(Tgau - Tsol))*(Sres - 1.0)) / (2.0 * (Tliq - Tsol));


		ISice =

			heaviside(Tgau - Tliq)*(Tgau*((3.0 * Sres) / 8.0 - 3.0 / 8.0) + ((Sres / 8.0 - 1.0 / 8.0)*
			(4.0 * cos((pi*(Tliq - 2.0 * Tgau + Tsol)) / (2.0 * (Tliq - Tsol))) + sin((pi*(Tliq - 2.0 * Tgau + Tsol))
				/ (Tliq - Tsol)) / 2.0)*(Tliq - Tsol)) / pi) - heaviside(Tgau - Tsol)*(Tgau*((3.0 * Sres) / 8.0 - 3.0 / 8.0)
					+ ((Sres / 8.0 - 1.0 / 8.0)*(4.0 * cos((pi*(Tliq - 2.0 * Tgau + Tsol)) / (2.0 * (Tliq - Tsol))) +
						sin((pi*(Tliq - 2.0 * Tgau + Tsol)) / (Tliq - Tsol)) / 2.0)*(Tliq - Tsol)) / pi) + Tsol
			* ((3.0 * Sres) / 8.0 - 3.0 / 8.0) + Tgau * (heaviside(Tgau - Tsol) - 1.0)*(Sres - 1.0) - Tliq *
			heaviside(Tliq - Tsol)*((3.0 * Sres) / 8.0 - 3.0 / 8.0) - Tsol * (heaviside(Tgau - Tsol) - 1.0)*(Sres - 1.0)
			- Tliq * ((3.0 * Sres) / 8.0 - 3.0 / 8.0)*(heaviside(Tgau - Tliq) - 1.0) + Tsol * ((3.0 * Sres) / 8.0 - 3.0 / 8.0)*
			(heaviside(Tgau - Tsol) - 1.0) + Tsol * ((3.0 * Sres) / 8.0 - 3.0 / 8.0)*(heaviside(Tliq - Tsol) - 1.0);


		IdSice =

			(Sres / 4.0 - 1.0 / 4.0)*(heaviside(Tgau - Tliq) + 3.0 * heaviside(Tgau - Tsol) + 4.0 * heaviside(Tliq - Tsol)
				+ sin((pi*(Tliq - 2.0 * Tgau + Tsol)) / (2.0 * (Tliq - Tsol)))*sin((pi*(Tliq - 2.0 * Tgau + Tsol)) / (2.0 * (Tliq - Tsol))) * heaviside(Tgau - Tliq)
				- sin((pi*(Tliq - 2.0 * Tgau + Tsol)) / (2.0 * (Tliq - Tsol)))*sin((pi*(Tliq - 2.0 * Tgau + Tsol)) / (2.0 * (Tliq - Tsol))) * heaviside(Tgau - Tsol) +
				2.0 * sin((pi*(Tliq - 2.0 * Tgau + Tsol)) / (2.0 * (Tliq - Tsol)))*heaviside(Tgau - Tliq) - 2.0
				* sin((pi*(Tliq - 2.0 * Tgau + Tsol)) / (2.0 * (Tliq - Tsol)))*heaviside(Tgau - Tsol) - 4.0);


		IdSicexT =

			((Sres - 1.0)*(8.0 * Tliq*cos((pi*(Tliq - 2.0 * Tgau + Tsol)) / (2.0 * (Tliq - Tsol)))*heaviside(Tgau - Tsol)
				- 10.0 * pi*Tsol - 8.0 * Tliq*cos((pi*(Tliq - 2.0 * Tgau + Tsol)) / (2.0 * (Tliq - Tsol)))*heaviside(Tgau - Tliq)
				- 6.0 * pi*Tliq + 8.0 * Tsol*cos((pi*(Tliq - 2.0 * Tgau + Tsol)) / (2.0 * (Tliq - Tsol)))*heaviside(Tgau - Tliq)
				- 8.0 * Tsol*cos((pi*(Tliq - 2.0 * Tgau + Tsol)) / (2.0 * (Tliq - Tsol)))*heaviside(Tgau - Tsol) - Tliq *
				sin((pi*(Tliq - 2.0 * Tgau + Tsol)) / (Tliq - Tsol))*heaviside(Tgau - Tliq) + Tliq * sin((pi*(Tliq - 2.0
					* Tgau + Tsol)) / (Tliq - Tsol))*heaviside(Tgau - Tsol) + Tsol * sin((pi*(Tliq - 2.0 * Tgau + Tsol))
						/ (Tliq - Tsol))*heaviside(Tgau - Tliq) - Tsol * sin((pi*(Tliq - 2.0 * Tgau + Tsol)) / (Tliq - Tsol))
				*heaviside(Tgau - Tsol) + 6.0 * Tliq*pi*heaviside(Tgau - Tliq) + 6.0 * Tliq*pi*heaviside(Tliq - Tsol) +
				10.0 * Tsol*pi*heaviside(Tgau - Tsol) + 10.0 * Tsol*pi*heaviside(Tliq - Tsol) - 2.0 * Tgau*pi*cos((pi*(Tliq
					- 2.0 * Tgau + Tsol)) / (Tliq - Tsol))*heaviside(Tgau - Tliq) + 2.0 * Tgau*pi*cos((pi*(Tliq - 2.0 *
						Tgau + Tsol)) / (Tliq - Tsol))*heaviside(Tgau - Tsol) + 8.0 * Tgau*pi*sin((pi*(Tliq - 2.0 * Tgau + Tsol))
							/ (2.0 * (Tliq - Tsol)))*heaviside(Tgau - Tliq) - 8.0 * Tgau*pi*sin((pi*(Tliq - 2.0 * Tgau + Tsol)) /
							(2.0 * (Tliq - Tsol)))*heaviside(Tgau - Tsol))) / (16.0 * pi);


		ISicexT =

			(Tsol*Tsol * (Sres / 4.0 - 1.0 / 4.0)) / 4.0 + Tsol * Tsol * (Sres / 8.0 - 1.0 / 8.0) +
			heaviside(Tgau - Tliq)*(Tgau*Tgau * (Sres / 8.0 - 1.0 / 8.0) +
			(Tgau*Tgau * cos((pi*(Tliq - 2.0 * Tgau + Tsol)) /
				(2.0 * (Tliq - Tsol)))*cos((pi*(Tliq - 2.0 * Tgau + Tsol)) /
				(2.0 * (Tliq - Tsol))) * (Sres / 4.0 - 1.0 / 4.0)) / 4.0 + (Tgau*Tgau * sin((pi*(Tliq - 2.0 * Tgau + Tsol))
					/ (2.0 * (Tliq - Tsol)))*sin((pi*(Tliq - 2.0 * Tgau + Tsol))
						/ (2.0 * (Tliq - Tsol))) * (Sres / 4.0 - 1.0 / 4.0)) / 4.0 + (sin((pi*(Tliq - 2.0 * Tgau + Tsol))
							/ (2.0 * (Tliq - Tsol)))*(Tliq - Tsol)*(Tliq - Tsol) * (Sres - 1.0)) / (2.0 * pi*pi) + (sin((pi*(Tliq - 2.0 * Tgau + Tsol))
								/ (2.0 * (Tliq - Tsol)))*sin((pi*(Tliq - 2.0 * Tgau + Tsol))
									/ (2.0 * (Tliq - Tsol))) * (Tliq - Tsol)*(Tliq - Tsol) * (Sres - 1.0)) / (16.0 * pi*pi) +
									(Tgau*cos((pi*(Tliq - 2.0 * Tgau + Tsol)) / (2.0 * (Tliq - Tsol)))*(Tliq - Tsol)*(Sres - 1.0)) /
				(2.0 * pi) + (Tgau*sin((pi*(Tliq - 2.0 * Tgau + Tsol)) / (Tliq - Tsol))*(Tliq - Tsol)*(Sres - 1.0)) / (16.0 * pi))
			- heaviside(Tgau - Tsol)*(Tgau*Tgau * (Sres / 8.0 - 1.0 / 8.0) + (Tgau*Tgau * cos((pi*(Tliq - 2.0 * Tgau + Tsol)) /
			(2.0 * (Tliq - Tsol)))*cos((pi*(Tliq - 2.0 * Tgau + Tsol)) /
				(2.0 * (Tliq - Tsol))) * (Sres / 4.0 - 1.0 / 4.0)) / 4.0 + (Tgau*Tgau * sin((pi*(Tliq - 2.0 * Tgau + Tsol)) /
				(2.0 * (Tliq - Tsol)))*sin((pi*(Tliq - 2.0 * Tgau + Tsol)) /
					(2.0 * (Tliq - Tsol))) * (Sres / 4.0 - 1.0 / 4.0)) / 4.0 + (sin((pi*(Tliq - 2.0 * Tgau + Tsol)) / (2.0 * (Tliq - Tsol)))
						*(Tliq - Tsol)*(Tliq - Tsol) * (Sres - 1.0)) / (2.0 * pi*pi) + (sin((pi*(Tliq - 2.0 * Tgau + Tsol))
							/ (2.0 * (Tliq - Tsol)))*sin((pi*(Tliq - 2.0 * Tgau + Tsol))
								/ (2.0 * (Tliq - Tsol))) * (Tliq - Tsol)*(Tliq - Tsol) * (Sres - 1.0)) / (16.0 * pi*pi) +
								(Tgau*cos((pi*(Tliq - 2.0 * Tgau + Tsol)) / (2.0 * (Tliq - Tsol)))*(Tliq - Tsol)*(Sres - 1.0))
				/ (2.0 * pi) + (Tgau*sin((pi*(Tliq - 2.0 * Tgau + Tsol)) / (Tliq - Tsol))*(Tliq - Tsol)*(Sres - 1.0))
				/ (16.0 * pi)) + Tgau * Tgau * (Sres / 2.0 - 1.0 / 2.0)*(heaviside(Tgau - Tsol) - 1.0) - Tsol * Tsol *
				(Sres / 2.0 - 1.0 / 2.0)*(heaviside(Tgau - Tsol) - 1.0) + (9.0 * (Tliq - Tsol)*(Tliq - Tsol) * (Sres - 1.0)) / (16.0 * pi*pi) -
			((heaviside(Tgau - Tliq) - 1.0)*(Sres - 1.0)*(14.0 * Tliq*Tsol - 7.0 * Tliq*Tliq - 7.0 * Tsol*Tsol + 3.0 * Tliq*Tliq * pi*pi)) /
			(16.0 * pi*pi) + (3.0 * (heaviside(Tgau - Tsol) - 1.0)*(Sres - 1.0)*(3.0 * Tliq*Tliq - 6.0 * Tliq*Tsol +
				3.0 * Tsol*Tsol + Tsol * Tsol * pi*pi)) / (16.0 * pi*pi) + (3.0 * (heaviside(Tliq - Tsol) - 1.0)*(Sres - 1.0)
					*(3.0 * Tliq*Tliq - 6.0 * Tliq*Tsol + 3.0 * Tsol*Tsol + Tsol * Tsol * pi*pi)) / (16.0 * pi*pi) -
					(heaviside(Tliq - Tsol)*(Sres - 1.0)*(14.0 * Tliq*Tsol - 7.0 * Tliq*Tliq - 7.0 * Tsol*Tsol + 3.0 * Tliq*Tliq * pi*pi)) / (16.0 * pi*pi);


		Sair =

			heaviside(Tsol - Tgau)*(Sres - 1.0) - Sres - (sin((pi*(Tliq / 2.0 - Tgau + Tsol / 2.0)) / (Tliq - Tsol)) / 2.0 - 1.0 / 2.0)
			*(sin((pi*(Tliq / 2.0 - Tgau + Tsol / 2.0)) / (Tliq - Tsol)) / 2.0 + 1.0 / 2.0)*(heaviside(Tgau - Tliq) - heaviside(Tgau - Tsol))
			- (sin((pi*(Tliq / 2.0 - Tgau + Tsol / 2.0)) / (Tliq - Tsol)) / 2.0 + 1.0 / 2.0)*(heaviside(Tgau - Tliq) - heaviside(Tgau - Tsol))
			*(Sres / 2.0 + sin((pi*(Tliq / 2.0 - Tgau + Tsol / 2.0)) / (Tliq - Tsol))*(Sres / 2.0 - 1.0 / 2.0) - 1.0 / 2.0) + 1.0;


		ISair =

			Tgau - (7.0 * Tsol) / 8.0 - Sres * Tgau + Sres * Tsol - heaviside(Tgau - Tliq)*(Tgau*((3.0 * Sres) / 8.0 - 3.0 / 8.0)
				+ ((Sres / 8.0 - 1.0 / 8.0)*(4.0 * cos((pi*(Tliq - 2.0 * Tgau + Tsol)) / (2.0 * (Tliq - Tsol))) + sin((pi*(Tliq - 2.0 * Tgau + Tsol))
					/ (Tliq - Tsol)) / 2.0)*(Tliq - Tsol)) / pi) + heaviside(Tgau - Tsol)*(Tgau*((3.0 * Sres) / 8.0 - 3.0 / 8.0) + ((Sres / 8.0 - 1.0 / 8.0)
						*(4.0 * cos((pi*(Tliq - 2.0 * Tgau + Tsol)) / (2.0 * (Tliq - Tsol))) + sin((pi*(Tliq - 2.0 * Tgau + Tsol)) / (Tliq - Tsol)) / 2.0)
						*(Tliq - Tsol)) / pi) + heaviside(Tgau - Tliq)*(Tgau / 8.0 - (sin((pi*(Tliq - 2.0 * Tgau + Tsol)) / (Tliq - Tsol))*(Tliq - Tsol))
							/ (16.0 * pi)) - heaviside(Tgau - Tsol)*(Tgau / 8.0 - (sin((pi*(Tliq - 2.0 * Tgau + Tsol)) / (Tliq - Tsol))*(Tliq - Tsol))
								/ (16.0 * pi)) - (3.0 * Tsol*(Sres - 1.0)) / 8.0 - (Tliq*heaviside(Tliq - Tsol)) / 8.0 - (Tliq*(heaviside(Tgau - Tliq) - 1.0))
			/ 8.0 + (Tsol*(heaviside(Tgau - Tsol) - 1.0)) / 8.0 + (Tsol*(heaviside(Tliq - Tsol) - 1.0)) / 8.0 - Tgau * (heaviside(Tgau - Tsol) - 1.0)*(Sres - 1.0)
			+ Tliq * heaviside(Tliq - Tsol)*((3.0 * Sres) / 8.0 - 3.0 / 8.0) + Tsol * (heaviside(Tgau - Tsol) - 1.0)*(Sres - 1.0) + Tliq * ((3.0 * Sres) / 8.0 - 3.0 / 8.0)
			*(heaviside(Tgau - Tliq) - 1.0) - Tsol * ((3.0 * Sres) / 8.0 - 3.0 / 8.0)*(heaviside(Tgau - Tsol) - 1.0) - Tsol * ((3.0 * Sres) / 8.0 - 3.0 / 8.0)*(heaviside(Tliq - Tsol) - 1.0);


		ISairxT =

			(Tsol*Tsol * (Sres / 4.0 - 1.0 / 4.0)) / 4.0 + Tsol * Tsol * (Sres / 8.0 - 1.0 / 8.0) + heaviside(Tgau - Tliq)*(Tgau*Tgau * (Sres / 8.0 - 1.0 / 8.0) +
			(Tgau*Tgau * cos((pi*(Tliq - 2.0 * Tgau + Tsol)) / (2.0 * (Tliq - Tsol)))*cos((pi*(Tliq - 2.0 * Tgau + Tsol)) / (2.0 * (Tliq - Tsol))) * (Sres / 4.0 - 1.0 / 4.0)) / 4.0
				+ (Tgau*Tgau * sin((pi*(Tliq - 2.0 * Tgau + Tsol))
					/ (2.0 * (Tliq - Tsol)))*sin((pi*(Tliq - 2.0 * Tgau + Tsol))
						/ (2.0 * (Tliq - Tsol))) * (Sres / 4.0 - 1.0 / 4.0)) / 4.0 + (sin((pi*(Tliq - 2.0 * Tgau + Tsol)) / (2.0 * (Tliq - Tsol)))*(Tliq - Tsol)*(Tliq - Tsol) * (Sres - 1.0))
				/ (2.0 * pi*pi) + (sin((pi*(Tliq - 2.0 * Tgau + Tsol)) / (2.0 * (Tliq - Tsol)))*sin((pi*(Tliq - 2.0 * Tgau + Tsol)) / (2.0 * (Tliq - Tsol))) * (Tliq - Tsol)*(Tliq - Tsol) * (Sres - 1.0)) / (16.0 * pi*pi) +
				(Tgau*cos((pi*(Tliq - 2.0 * Tgau + Tsol)) / (2.0 * (Tliq - Tsol)))*(Tliq - Tsol)*(Sres - 1.0)) / (2.0 * pi) + (Tgau*sin((pi*(Tliq - 2.0 * Tgau + Tsol))
					/ (Tliq - Tsol))*(Tliq - Tsol)*(Sres - 1.0)) / (16.0 * pi)) - heaviside(Tgau - Tsol)*(Tgau*Tgau * (Sres / 8.0 - 1.0 / 8.0) + (Tgau*Tgau * cos((pi*(Tliq - 2.0 * Tgau + Tsol))
						/ (2.0 * (Tliq - Tsol)))*cos((pi*(Tliq - 2.0 * Tgau + Tsol))
							/ (2.0 * (Tliq - Tsol))) * (Sres / 4.0 - 1.0 / 4.0)) / 4.0 + (Tgau*Tgau * sin((pi*(Tliq - 2.0 * Tgau + Tsol)) / (2.0 * (Tliq - Tsol)))*sin((pi*(Tliq - 2.0 * Tgau + Tsol)) / (2.0 * (Tliq - Tsol))) * (Sres / 4.0 - 1.0 / 4.0)) / 4.0 +
							(sin((pi*(Tliq - 2.0 * Tgau + Tsol)) / (2.0 * (Tliq - Tsol)))*(Tliq - Tsol)*(Tliq - Tsol) * (Sres - 1.0)) / (2.0 * pi*pi) + (sin((pi*(Tliq - 2.0 * Tgau + Tsol)) / (2.0 * (Tliq - Tsol)))*sin((pi*(Tliq - 2.0 * Tgau + Tsol)) / (2.0 * (Tliq - Tsol))) *
						(Tliq - Tsol)*(Tliq - Tsol) * (Sres - 1.0)) / (16.0 * pi*pi) + (Tgau*cos((pi*(Tliq - 2.0 * Tgau + Tsol)) / (2.0 * (Tliq - Tsol)))*(Tliq - Tsol)*(Sres - 1.0)) / (2.0 * pi) +
								(Tgau*sin((pi*(Tliq - 2.0 * Tgau + Tsol)) / (Tliq - Tsol))*(Tliq - Tsol)*(Sres - 1.0)) / (16.0 * pi)) + Tgau * Tgau * (Sres / 2.0 - 1.0 / 2.0)*(heaviside(Tgau - Tsol) - 1.0)
			- Tsol * Tsol * (Sres / 2.0 - 1.0 / 2.0)*(heaviside(Tgau - Tsol) - 1.0) + (9 * (Tliq - Tsol)*(Tliq - Tsol) * (Sres - 1.0)) / (16.0 * pi*pi) - ((heaviside(Tgau - Tliq) - 1.0)*(Sres - 1.0)
				*(14 * Tliq*Tsol - 7.0 * Tliq*Tliq - 7.0 * Tsol*Tsol + 3.0 * Tliq*Tliq * pi*pi)) / (16.0 * pi*pi) + (3.0 * (heaviside(Tgau - Tsol) - 1.0)*
				(Sres - 1.0)*(3.0 * Tliq*Tliq - 6 * Tliq*Tsol + 3.0 * Tsol*Tsol + Tsol * Tsol * pi*pi)) / (16.0 * pi*pi) + (3.0 * (heaviside(Tliq - Tsol)
					- 1.0)*(Sres - 1.0)*(3.0 * Tliq*Tliq - 6 * Tliq*Tsol + 3.0 * Tsol*Tsol + Tsol * Tsol * pi*pi)) / (16.0 * pi*pi) -
					(heaviside(Tliq - Tsol)*(Sres - 1.0)*(14 * Tliq*Tsol - 7.0 * Tliq*Tliq - 7.0 * Tsol*Tsol + 3.0 * Tliq*Tliq * pi*pi)) / (16.0 * pi*pi);
	}
}