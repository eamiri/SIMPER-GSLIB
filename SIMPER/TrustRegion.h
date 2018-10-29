#ifndef __TRUSTREGION__
#define __TRUSTREGION__

#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <Eigen/SparseCholesky>

using namespace std;
using namespace Eigen;

class TrustRegion
{
public:
	VectorXd deltaTempTest;
	bool IsFullNewtonRaphson;
	bool IsTRSuccess;
	double Error;

	void Iterate(SparseMatrix<double> B, VectorXd R, double TRRadius);

private:
	double calculateTrace(SparseMatrix<double> sm);
	SparseMatrix<double> ShiftDiag(SparseMatrix<double> sm, double lambda);
};

#endif