#include "TrustRegion.h"

void TrustRegion::Iterate(SparseMatrix<double> B, VectorXd R, double TRRadius)
{
	B = 0.5 * (SparseMatrix<double>(B.transpose()) + B);
	SimplicialLLT<SparseMatrix<double> > Solver;
	int n = R.size();
	double lambdaLow = 0.0;
	double lambda = 0.0;
	double traceBxBtran = calculateTrace(B * B);
	double lambdaHigh = 1.1*sqrt(traceBxBtran);
	int maxIterations = 40;
	IsFullNewtonRaphson = true;
	IsTRSuccess = true;
	int iIteration = 0;
	//B = (B + B.transpose()) / 2.0;
	deltaTempTest = deltaTempTest.Zero(n);
	Error = 0.0;
	MatrixXd eye, L;
	SparseMatrix<double> shiftedB, shiftedBInv;
	LLT<MatrixXd> fact;
	bool isCholSuccess = true;
	bool increaseLambda = false;
	VectorXd q;
	double consvio, deltaLambda, normDeltaTempTest, normQ;
	double tolError = 1E-3;

	while (true)
	{
		iIteration++;
		if (iIteration > maxIterations)
		{
			IsFullNewtonRaphson = false;
			IsTRSuccess = false;
			break;
		}

		Solver.setShift(lambda);
		Solver.compute(B);

		if (Solver.info() == Success)
		{
			deltaTempTest = -(Solver.solve(R));
			q = Solver.matrixL().solve(deltaTempTest);
			normQ = q.norm();
			normDeltaTempTest = deltaTempTest.norm();
			consvio = normDeltaTempTest - TRRadius;
			deltaLambda = ((consvio) / TRRadius) * (normDeltaTempTest / normQ) * (normDeltaTempTest / normQ);
			Error = consvio / TRRadius;

			if (lambda != 0.0)
			{
				IsFullNewtonRaphson = false;
			}

			if (abs(Error) <= tolError || (consvio <= 0.0 && lambda == 0.0))
			{
				break;
			}

			lambda += deltaLambda;
			increaseLambda = (consvio > 0.0);
		}
		else
		{
			increaseLambda = true;
		}

		if (increaseLambda)
		{
			lambdaLow = lambda;
		}
		else
		{
			lambdaHigh = lambda;
		}

		if (Solver.info() != Success || lambda < lambdaLow || lambda > lambdaHigh)
		{
			lambda = (lambdaLow + lambdaHigh) / 2.0;
		}
	}
}

double TrustRegion::calculateTrace(SparseMatrix<double> sm)
{
	double trace = 0.0;
	for (int k = 0; k < sm.outerSize(); k++)
	{
		trace += sm.coeff(k, k);
	}

	return trace;
}

SparseMatrix<double> TrustRegion::ShiftDiag(SparseMatrix<double> sm, double lambda)
{
	for (int k = 0; k < sm.outerSize(); k++)
	{
		sm.coeffRef(k, k) += lambda;
	}

	return sm;
}