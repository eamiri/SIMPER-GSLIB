#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <Eigen/SparseCholesky>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <string>
#include <strstream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <cstdlib>
#include <cstring>
#include <algorithm>
#include <math.h>
#include <limits>
#include <vector>
#include "Mesh.h"
#include "IBConditions.h"
#include "Properties.h"
#include "GaussPoints.h"
#include "CommonFunctions.h"
#include "FEMFunctions.h"
#include "TrustRegion.h"

using namespace std;
using namespace Eigen;

extern int noel;
extern int nond;
extern int ndoe;
extern int nGP;
extern Mesh MESH;
extern VectorXd GSLIBCoeffs;
extern Properties PROPS;
extern GaussPoints GP;
extern VectorXd Temp;
extern VectorXd TempDot;
extern VectorXd Temp_0;
extern VectorXd TempDot_0;
extern VectorXd Residual;
extern vector<int> DirichletDof;
extern MatrixXd BCInputData;
const double  PI = 3.1415926535898; ///< Double approximation of pi

extern FILE *OutputFile, *NodePlotFile;
extern ofstream outputFile;


struct GlobalOptions
{
	string Version;
	string OutputDirectory;
	string OutputFileName;

	bool Pause;
	bool Silent;
};

enum exitcode
{
	BAD_DATA,       ///< For bad input provided by user (requires immediate exit from program)
	BAD_DATA_WARN,  ///< For bad input provided by user (requires shutdown prior to simulation)
	RUNTIME_ERR,    ///< For runtime error (bad programming)
	FILE_OPEN_ERR,  ///< For bad file open (requires immediate exit)
	SIMPER_OPEN_ERR, ///< for bad SimPerErrors.txt file open
	OUT_OF_MEMORY,  ///< When out of memory
	SIMULATION_DONE ///< Upon completion of the simulation
};

///////////////////////////////////////////////////////////////////
/// \brief Converts string parameter to integer type
/// \param *s1 [in] String to be converted to integer
/// \return Integer format of passed string
//
inline int StringToInteger(const char *str) { return (int)atof(str); }