#include "SimperInclude.h"

double Det2x2(vector<vector<double>> A)
{
	return A[1][1] * A[2][2] - A[1][2] * A[2][1];
}

vector<vector<double>> Inv2x2(vector<vector<double>> A)
{
	vector<vector<double>> Ap;
	Ap.resize(2);
	for (int i=0; i < 2; i++)
	{
		Ap[i].resize(2);
	}

	double detA = Det2x2(A);
	Ap[1][1] = detA * A[2][2];
	Ap[2][2] = detA * A[1][1];
	Ap[1][2] = -detA * A[1][2];
	Ap[2][1] = -detA * A[2][1];

	return Ap;
}

double heaviside(double x)
{
	double hs = 0.0;
	if (x > 0)
	{
		hs = 1.0;
	}
	else
	{
		hs = 0.0;
	}

	return hs;
}

double sign(double x)
{
	double s;
	if (x > 0)
	{
		s = 1.0;
	}
	else if (x == 0)
	{
		s = 0.0;
	}
	else
	{
		s = -1.0;
	}

	return s;
}