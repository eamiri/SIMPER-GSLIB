//#include <Eigen/Core>
//
//using Eigen::MatrixXd;
//using Eigen::RowVectorXd;
//
//using namespace std; 
//
//RowVectorXd FourNoded(double s, double r)
//{
//	RowVectorXd Nmat(1, 4);
//	Nmat << 0.25*(1 - s)*(1 - r), 0.25*(1 + s)*(1 - r), 0.25*(1 + s)*(1 + r), 0.25*(1 - s)*(1 + r);
//
//	return Nmat;
//}
//
//MatrixXd dFourNoded(double s, double r)
//{
//	MatrixXd dNMat(1, 4);
//	RowVectorXd dNdr(1,4), dNds(1,4);
//	dNdr << s / 4 - 1 / 4, -s / 4 - 1 / 4, s / 4 + 1 / 4, 1 / 4 - s / 4 ;
//	dNds << r / 4 - 1 / 4, 1 / 4 - r / 4, r / 4 + 1 / 4, -r / 4 - 1 / 4 ;
//
//	dNMat << dNdr,
//			 dNds;
//
//	return dNMat;
//}
//
////vector<vector<double>> JFourNoded(double s, double r)
////{
////	vector<vector<double>> dN = dFourNoded(s, r);
////	double *dNdr = dN[0];
////	double *d
////	return { dNds, dNds };
////}
//
////vector<vector<double>> BFourNoded(double s, double r)
////{
////	vector<double> dNdr = { s / 4 - 1 / 4, -s / 4 - 1 / 4, s / 4 + 1 / 4, 1 / 4 - s / 4 };
////	vector<double> dNds = { r / 4 - 1 / 4, 1 / 4 - r / 4, r / 4 + 1 / 4, -r / 4 - 1 / 4 };
////	return { dNds, dNds };
////}