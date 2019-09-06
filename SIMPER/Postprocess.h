#ifndef __POSTPROC__
#define __POSTPROC__

#include <Eigen/Core>
#include <vector>
#include "Mesh.h"
#include "Properties.h"
#include "FEMFunctions.h"
#include "SaturationFunctions.h"
#include "GaussPoints.h"
#include <iomanip>
#include <iostream>
#include <fstream>
#include <string>
#include <strstream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <cstring>
#include "OutputFileStruc.h"

using namespace std;
using namespace Eigen;

class Postprocess
{
public:
	VectorXd Temp;
	double SolutionTime;

	Postprocess(Mesh mesh, Properties props, MatrixXd nodalGSLIBCoeffs, OutputFileStruc files);
	void Plot(VectorXd Temp, double SolutionTime, int iRealization);
	void AreaAnalysis(VectorXd Temp, double solutionTime, int iRealization);
	void GetTalikArea(VectorXd minTemp, VectorXd maxTemp, int year, int iRealization);
	void GetPermafrostArea(VectorXd minTemp, VectorXd maxTemp, int year, int iRealization);
	void PlotNodes(VectorXd Temp, double solutionTime, int iRealization);

private:
	VectorXd NodalGSLIBCoeffs;
	Mesh MESH;
	Properties PROPS;
	OutputFileStruc Files;
	int nGP;
	GaussPoints GP;
	int noel;
	int nond;
	int ndoe;

	VectorXd xNodes, yNodes, waterSat, iceSat, distance;
	VectorXd TempNode;
	VectorXd GradTemp;
	VectorXi elementDofs;
	MatrixXd derTemp;
	FEMFunctions femFunctions;
	MatrixXd Bmat;
	RowVectorXd SF;

	double GPi, GPj, Wi, Wj, globalGPx, globalGPy, dNodeGP, xNode, yNode, AreaElement;
	int rSFC;
	double Tsol, Tliq, Sres, Wpar, Mpar, SwSol, SwLiq;

	double frozenArea, thawedArea, slushyArea;
	double talikArea, permafrostArea;
	int iFrozen, iThawed, iSlushy;
	VectorXd minTempNode;
	VectorXd maxTempNode;
};

#endif