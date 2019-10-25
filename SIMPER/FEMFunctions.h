#ifndef __FENFUNCTIONS__
#define __FENFUNCTIONS__

#include <Eigen/Core>
#include <Eigen/LU>

using namespace std;
using namespace Eigen;

class FEMFunctions
{
public:
	RowVectorXd SF;
	MatrixXd dSF;
	MatrixXd J;
	MatrixXd Bmat;

	void Calculate (VectorXd Xnode, VectorXd Ynode, double s, double r)
	{
		SF = FourNodedSF(s, r);
		dSF = dFourNodedSF(s, r);
		J = JacobianFourNoded(Xnode, Ynode, s, r);
		Bmat = BMatrixFourNoded(Xnode, Ynode, s, r);
	}

	VectorXd CalculateGlobalCoordinates(VectorXd SF, VectorXd Xnode, VectorXd Ynode)
	{
		VectorXd XYGlob;
		XYGlob.resize(2);
		XYGlob.setZero();

		for (int i = 0; i < SF.size(); i++)
		{
			XYGlob(0) += SF(i) * Xnode(i);
			XYGlob(1) += SF(i) * Ynode(i);
		}

		return XYGlob;
	}

private:
	RowVectorXd FourNodedSF (double s, double r)
	{
		RowVectorXd N(4);
		N << 0.25*(1 - s)*(1 - r), 0.25*(1 + s)*(1 - r), 0.25*(1 + s)*(1 + r), 0.25*(1 - s)*(1 + r);

		return N;
	}

	MatrixXd dFourNodedSF(double s, double r)
	{
		RowVectorXd dNdr(1, 4), dNds(1, 4);
		dNdr << 0.25 * s - 0.25, -0.25 * s - 0.25, 0.25 * s + 0.25, 0.25 - 0.25 * s;
		dNds << 0.25 * r - 0.25, 0.25 - 0.25 * r, 0.25 * r + 0.25, -0.25 * r - 0.25;

		MatrixXd dN(2, 4);
		dN << dNds,
			  dNdr;

		return dN;
	}
	
	MatrixXd JacobianFourNoded(VectorXd Xnode, VectorXd Ynode, double s, double r)
	{
		MatrixXd XYnode(4, 2);
		XYnode << Xnode, Ynode;

		MatrixXd dN(2, 4);
		dN << dFourNodedSF(s, r);
		
		return dN * XYnode;
	}

	MatrixXd BMatrixFourNoded(VectorXd Xnode, VectorXd Ynode, double s, double r)
	{
		MatrixXd J;
		J = JacobianFourNoded(Xnode, Ynode, s, r);

		MatrixXd dN(2, 4);
		dN << dFourNodedSF(s, r);

		return J.inverse() * dN;
	}
};

#endif