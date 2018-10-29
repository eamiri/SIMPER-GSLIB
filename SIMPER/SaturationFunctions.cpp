#include "SimperInclude.h"

int *rSFC = &PROPS.Soil.rSFC;
double *Tsol = &PROPS.Nonisothermal.TempSolid;
double *Tliq = &PROPS.Nonisothermal.TempLiquid;
double *Sres = &PROPS.Soil.ResidualWaterSaturation;
double *Wpar = &PROPS.Soil.Wpar;
double *Mpar = &PROPS.Soil.Mpar;

double SATUR(double Tgau)
{
	double satur, Btem;
	switch (*rSFC)
	{
	case 1:
		if (Tgau <= *Tsol)
		{
			satur = *Sres;
		}
		else if (Tgau < *Tliq && Tgau > *Tsol)
		{
			satur = (1.0 - *Sres)*exp(-((Tgau - *Tliq) / *Wpar)*((Tgau - *Tliq) / *Wpar)) + *Sres;
		}
		else
		{
			satur = 1.0;

		}

		break;

	case 0:
		Btem = (*Sres - 1) / *Mpar;
		if (Tgau >= Btem)
		{
			satur = *Mpar * Tgau + 1.0;
			if (Tgau > *Tliq)
			{
				satur = 1.0;
			}
		}
		else
		{
			satur = *Sres;
		}

		break;

	case 2:
		if (Tgau <= *Tsol)
		{
			satur = *Sres;
		}
		else if (Tgau <*Tliq && Tgau >* Tsol)
		{
			satur = 0.5*(1.0 + *Sres) + 0.5*(1.0 - *Sres)*sin(PI*(Tgau - 0.5*(*Tliq + *Tsol)) / (*Tliq - *Tsol));
		}
		else
		{
			satur = 1.0;
		}

		break;
	}

	return satur;
}

double dSATUR(double Tgau)
{
	double dSatur, Btem;

	switch (*rSFC)
	{
	case 1:
		if (Tgau <= *Tsol)
		{
			dSatur = 0.0;
		}
		else if (Tgau < *Tliq && Tgau > *Tsol)
		{
			dSatur = (exp(-(Tgau - *Tliq) * (Tgau - *Tliq) / *Wpar * *Wpar)*(*Sres - 1.0)*(2.0 * Tgau - 2.0 * *Tliq)) / (*Wpar * *Wpar);
		}
		else
		{
			dSatur = 0;
		}

		break;

	case 0:
		Btem = (*Sres - 1.0) / *Mpar;
		if (Tgau >= Btem)
		{
			dSatur = *Mpar;
			if (Tgau > *Tliq)
			{
				dSatur = 0.0;
			}
		}
		else
		{
			dSatur = 0.0;
		}

		break;

	case 2:
		if (Tgau <= *Tsol)
		{
			dSatur = *Sres;
		}
		else if (Tgau < *Tliq && Tgau > *Tsol)
		{
			dSatur = -(PI*cos((PI*(*Tliq / 2.0 - Tgau + *Tsol / 2.0)) / (*Tliq - *Tsol))*(*Sres / 2.0 - 1.0 / 2.0)) / (*Tliq - *Tsol);
		}
		else
		{
			dSatur = 0.0;
		}

		break;
	}

	return dSatur;
}

double ISATUR(double Tgau)
{
	double iSatur, CSw1, CSw2, CSw1p, Btem;
	switch (*rSFC)
	{
	case 1:
		if (Tgau <= *Tsol)
		{
			iSatur = *Sres * Tgau;
		}
		else if (Tgau < *Tliq && Tgau > *Tsol)
		{
			/*CSw1 = *Sres * *Tsol - (*Sres * *Tsol - (sqrt(PI)*erf((*Tsol - *Tliq)*(-1.0 / (*Wpar * *Wpar)) ^ (1 / 2))*(*Sres - 1)) / (2.0 * (-1.0 / (*Wpar * *Wpar)) ,(1.0 / 2.0)));
			iSatur = CSw1 + *Sres * Tgau - (sqrt(PI)*erf((Tgau - *Tliq)*(-1 / (*Wpar * *Wpar)) ^ (1 / 2))*(*Sres - 1.0)) / (2 * sqrt((-1.0 / *Wpar * *Wpar)));*/
		}

		break;

	case 0:
		Btem = (*Sres - 1.0) / *Mpar;
		if (Tgau >= Btem)
		{
			CSw1 = *Sres * Btem - (*Mpar*Btem*Btem / 2 + Btem);
			iSatur = CSw1 + *Mpar * Tgau * Tgau / 2.0 + Tgau;
			if (Tgau > *Tliq)
			{
				CSw2 = (CSw1 + *Mpar * *Tliq * *Tliq / 2.0 + *Tliq) - *Tliq;
				iSatur = CSw2 + Tgau;
			}
		}
		else
		{
			iSatur = *Sres * Tgau;
		}

		break;

	case 2:
		if (Tgau <= *Tsol)
		{
			iSatur = *Sres * Tgau;
		}
		else if (Tgau < *Tliq && Tgau > *Tsol)
		{
			CSw1 = *Sres * *Tsol - ((*Tsol*(*Sres + 1.0)) / 2.0 + ((*Tliq - *Tsol)*(*Sres - 1.0)) / (2.0 * PI));
			iSatur = CSw1 + (Tgau*(*Sres + 1.0)) / 2.0 + (cos((PI*(*Tliq - 2.0 * Tgau + *Tsol)) / (4 * (*Tliq - *Tsol))) * cos((PI*(*Tliq - 2.0 * Tgau + *Tsol)) / (4 * (*Tliq - *Tsol))) * (*Tliq - *Tsol)*(*Sres - 1.0)) / PI;
		}
		else
		{
			CSw1p = *Sres * *Tsol - ((*Tsol*(*Sres + 1.0)) / 2.0 + ((*Tliq - *Tsol)*(*Sres - 1.0)) / (2.0 * PI));
			CSw1 = CSw1p + (*Tliq*(*Sres + 1)) / 2 + ((*Tliq - *Tsol)*(*Sres - 1)) / (2.0 * PI);
			CSw2 = CSw1 - (*Tliq);
			iSatur = CSw2 + Tgau;
		}

		break;
	}

	return iSatur;
}

double ISATURxT(double Tgau)
{
	double iSatxT, CSw1, CSw2, CSw1p, Btem;

	switch (*rSFC)
	{
	case 0:
		Btem = (*Sres - 1.0) / *Mpar;
		if (Tgau >= Btem)
		{
			CSw1 = *Sres * Btem * Btem / 2.0 - (*Mpar*Btem*Btem*Btem / 3.0 + Btem * Btem / 2.0);
			iSatxT = CSw1 + *Mpar * Tgau * Tgau * Tgau / 3.0 + Tgau * Tgau / 2.0;
			if (Tgau > *Tliq)
			{
				CSw2 = (CSw1 + *Mpar * *Tliq * *Tliq * *Tliq / 3.0 + *Tliq * *Tliq / 2.0) - *Tliq * *Tliq / 2.0;
				iSatxT = CSw2 + Tgau * Tgau / 2;
			}
		}
		else
		{
			iSatxT = *Sres * Tgau * Tgau / 2;
		}

		break;

	case 2:
		if (Tgau <= *Tsol)
			iSatxT = *Sres * Tgau * Tgau / 2;
		else if (Tgau < *Tliq && Tgau > *Tsol)
		{
			CSw1 = *Sres * *Tsol * *Tsol / 2.0 - (*Tsol * *Tsol * (*Sres / 4.0 + 1.0 / 4.0) + ((*Tliq - *Tsol) * (*Tliq - *Tsol) * (*Sres - 1.0)) / (2.0 * PI * PI));
			iSatxT = CSw1 + Tgau * Tgau * (*Sres / 4.0 + 1.0 / 4.0) + (sin((PI*(*Tliq - 2.0 * Tgau + *Tsol)) / (2.0 * (*Tliq
				- *Tsol)))*(*Tliq - *Tsol) * (*Tliq - *Tsol) * (*Sres - 1.0)) / (2.0 * PI * PI) + (Tgau*cos((PI*(*Tliq
					- 2.0 * Tgau + *Tsol)) / (2.0 * (*Tliq - *Tsol)))*(*Tliq - *Tsol)*(*Sres - 1.0)) / (2.0 * PI);
		}
		else
		{
			CSw1p = *Sres * *Tsol * *Tsol / 2.0 - (*Tsol * *Tsol * (*Sres / 4.0 + 1.0 / 4.0) + ((*Tliq - *Tsol) * (*Tliq - *Tsol) * (*Sres - 1.0)) / (2.0 * PI * PI));
			CSw1 = CSw1p + *Tliq * *Tliq * (*Sres / 4.0 + 1.0 / 4.0) - ((*Tliq - *Tsol) * (*Tliq - *Tsol) * (*Sres - 1.0)) / (2.0 * PI * PI);
			CSw2 = CSw1 - (*Tliq * *Tliq / 2.0);
			iSatxT = CSw2 + Tgau * Tgau / 2.0;
		}

		break;
	}

	return iSatxT;
}

double IdSATURxT(double Tgau)
{
	double idSaturxT, IdSwxT_0, IdSwxT_0p, Btem;
	switch (*rSFC)
	{
	case 0:
		Btem = (*Sres - 1.0) / *Mpar;
		if (Tgau >= Btem)
		{
			IdSwxT_0 = *Mpar * Btem * Btem / 2.0;
			idSaturxT = *Mpar * Tgau * Tgau / 2.0 - IdSwxT_0;
			if (Tgau > *Tliq)
			{
				IdSwxT_0p = *Mpar * Btem * Btem / 2.0;
				IdSwxT_0 = *Mpar * *Tliq * *Tliq / 2.0 - IdSwxT_0p;
				idSaturxT = IdSwxT_0;
			}
		}
		else
		{
			idSaturxT = 0;
		}

		break;

	case 2:
		if (Tgau <= *Tsol)
			idSaturxT = *Sres * Tgau * Tgau / 2;
		else if (Tgau < *Tliq && Tgau > *Tsol)
		{
			IdSwxT_0 = -((cos((PI*(*Tliq / 2.0 - *Tsol / 2.0)) / (*Tliq -
				*Tsol))*(*Tliq - *Tsol) * (*Tliq - *Tsol) * (*Sres - 1)) / PI -
				*Tsol*sin((PI*(*Tliq / 2.0 - *Tsol / 2.0)) / (*Tliq -
					*Tsol))*(*Tliq - *Tsol)*(*Sres - 1)) / (2.0 * *Tliq - 2.0 * *Tsol);
			idSaturxT = (Tgau*sin((PI*(*Tliq / 2.0 - Tgau + *Tsol / 2.0)) / (*Tliq -
				*Tsol))*(*Tliq - *Tsol)*(*Sres - 1.0) -
				(cos((PI*(*Tliq / 2 - Tgau + *Tsol / 2)) / (*Tliq - *Tsol)) * (*Tliq - *Tsol)
					* (*Tliq - *Tsol) * (*Sres - 1.0)) / PI) / (2.0 * *Tliq - 2.0 * *Tsol) - IdSwxT_0;
		}
		else
		{
			IdSwxT_0p = -((cos((PI*(*Tliq / 2.0 - *Tsol / 2.0)) / (*Tliq -
				*Tsol))*(*Tliq - *Tsol) * (*Tliq - *Tsol)  * (*Sres - 1.0)) / PI -
				*Tsol*sin((PI*(*Tliq / 2.0 - *Tsol / 2.0)) / (*Tliq -
					*Tsol))*(*Tliq - *Tsol)*(*Sres - 1.0)) / (2.0 * *Tliq - 2.0 * *Tsol);
			IdSwxT_0 = -((cos((PI*(*Tliq / 2.0 - *Tsol / 2.0)) / (*Tliq -
				*Tsol))*(*Tliq - *Tsol) * (*Tliq - *Tsol) * (*Sres - 1.0)) / PI +
				*Tliq*sin((PI*(*Tliq / 2.0 - *Tsol / 2.0)) / (*Tliq -
					*Tsol))*(*Tliq - *Tsol)*(*Sres - 1.0)) / (2.0 * *Tliq - 2.0 * *Tsol) - IdSwxT_0p;
			idSaturxT = IdSwxT_0;
		}

		break;
	}

	return idSaturxT;
}