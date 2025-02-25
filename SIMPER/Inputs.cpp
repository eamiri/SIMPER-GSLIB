#ifdef _unix_        
#define linux
#elif defined(_WIN32) || defined(WIN32) 
#define _windows_
#endif

#include "SimperInclude.h"
#include <regex>
#include <experimental/filesystem>

bool GlobalFenIndex = false;
bool Inputs(string propsInputFile, string meshInputFile, int iPReal);
string GetOutputFilePath(string outputName, int iRealization);
string GetGSLIBOutputFilePath(string outputName, int iRealization);
Properties InputProperties(string filePath);

Mesh InputMesh(string filePath);
void SgsimParameterFile();
void BCSgsimParameterFile();
void GSLIBRunSGSIM();
void BCGSLIBRunSGSIM();
void UpscaleGSLIBtoSIMPER();
void AddcoorParameterFile(int nRealization);
void BCAddcoorParameterFile(int nRealization);
void VerticalIntegrationInfo();
void HorizontalIntegrationInfo();

int nRealization;
int noel;
int nond;
int ndoe;
int nGP;
Mesh MESH;
MatrixXd GSLIBCoeffs, BCGSLIBCoeffs;
MatrixXd NodalGSLIBCoeffs;
MatrixXd ElementalGSLIBCoeffs;
MatrixXd GSLIBGrid, BCGSLIBGrid;
Properties PROPS;
GaussPoints GP;
GaussPoints VIntGP;
GaussPoints HIntGP;

bool Inputs(string propsInputFile, string meshInputFile, int iPReal)
{
	iParallelRlzn = iPReal;
	PROPS = InputProperties(propsInputFile);
	MESH = InputMesh(meshInputFile);
	if (PROPS.VInteg.isIntegrate)
	{
		VerticalIntegrationInfo();
		HorizontalIntegrationInfo();
	}

	nRealization = PROPS.GSLIB.NumberOfRealizations;
	if (PROPS.GSLIB.isHeterBC)
	{
		BCSgsimParameterFile();
		BCGSLIBRunSGSIM();
	}
	
	if (PROPS.GSLIB.isGSLIB)
	{
		SgsimParameterFile();
		GSLIBRunSGSIM();
	}

	UpscaleGSLIBtoSIMPER();

	return true;
}

string GetGSLIBOutputFilePath(string outputName, int iRealization)
{
	string rlnNo = to_string(iRealization);
	string path = "GSLIB/outputs/R";
	path += rlnNo;
	path += "_";
	path += outputName;
	return path;
}

Properties InputProperties(string filePath)
{
	Properties props;
	ifstream propsFile;
	propsFile.open(filePath.c_str());
	if (propsFile.fail())
	{
		cout << "Cannot find file " << filePath << endl;
	}
	else
	{
		string line;
		getline(propsFile, line);
		getline(propsFile, line);
		getline(propsFile, line);
		istringstream isDataLine(line);
		isDataLine
			>> props.nSoilType
			>> props.GSLIB.isGSLIB
			>> props.IsLateralFlow;
		props.Soil.resize(props.nSoilType);
		getline(propsFile, line);
		for (int iSoilType = 0; iSoilType < props.nSoilType; iSoilType++)
		{
			getline(propsFile, line);
			istringstream isDataLine(line);
			isDataLine
				>> props.Soil[iSoilType].iType
				>> props.Soil[iSoilType].SoilType
				>> props.Soil[iSoilType].HydraulicConductivity
				>> props.Soil[iSoilType].HeatCapacity
				>> props.Soil[iSoilType].ThermalConductivity
				>> props.Soil[iSoilType].DensityMin
				>> props.Soil[iSoilType].DensityMax
				>> props.Soil[iSoilType].Porosity
				>> props.Soil[iSoilType].ResidualWaterSaturation
				>> props.Soil[iSoilType].rSFC
				>> props.Soil[iSoilType].Wpar
				>> props.Soil[iSoilType].Mpar
				>> props.Soil[iSoilType].IsSaturated
				>> props.Soil[iSoilType].FPmin
				>> props.Soil[iSoilType].FPmax;
		}

		getline(propsFile, line);
		getline(propsFile, line);
		getline(propsFile, line);
		isDataLine.clear();
		isDataLine.str(line);
		isDataLine
			>> props.Fluid.SHeatCapacity
			>> props.Fluid.SThermalConductivity
			>> props.Fluid.SDensity
			>> props.Fluid.LHeatCapacity
			>> props.Fluid.LThermalConductivity
			>> props.Fluid.LDensity
			>> props.Fluid.LatentHeat
			>> props.Fluid.AHeatCapacity
			>> props.Fluid.AThermalConductivity
			>> props.Fluid.ADensity;
		getline(propsFile, line);
		getline(propsFile, line);
		getline(propsFile, line);
		isDataLine.clear();
		isDataLine.str(line);
		isDataLine
			>> props.Nonisothermal.IsIsothermal
			>> props.Nonisothermal.HeatCapacity
			>> props.Nonisothermal.TempTransition
			>> props.Nonisothermal.ThermalConductivity
			>> props.Nonisothermal.TempSolid
			>> props.Nonisothermal.TempLiquid
			>> props.Nonisothermal.SolidSatIndex
			>> props.Nonisothermal.LiquidSatIndex;
		getline(propsFile, line);
		getline(propsFile, line);
		getline(propsFile, line);
		isDataLine.clear();
		isDataLine.str(line);
		isDataLine
			>> props.BCs.isICInput
			>> props.BCs.ICInputFile
			>> props.BCs.UniformIC;
		getline(propsFile, line);
		getline(propsFile, line);
		getline(propsFile, line);
		isDataLine.clear();
		isDataLine.str(line);
		isDataLine
			>> props.BCs.isBCInput
			>> props.BCs.BCInputFile;
		getline(propsFile, line);
		getline(propsFile, line);
		isDataLine.clear();
		isDataLine.str(line);
		isDataLine >> props.GQ.NumberOfPoints;
		GP.GP(props.GQ.NumberOfPoints);
		nGP = props.GQ.NumberOfPoints;
		getline(propsFile, line);
		getline(propsFile, line);
		getline(propsFile, line);
		isDataLine.clear();
		isDataLine.str(line);
		isDataLine
			>> props.Solution.MaxTimestep
			>> props.Solution.DeltaTime
			>> props.Solution.TolPsi
			>> props.Solution.MaxIterations
			>> props.Solution.NewmarkGamma
			>> props.Solution.PlotInterval
			>> props.Solution.IsGMSH;
		getline(propsFile, line);
		getline(propsFile, line);
		getline(propsFile, line);
		isDataLine.clear();
		isDataLine.str(line);
		isDataLine
			>> props.GSLIB.HorCorrelationLength
			>> props.GSLIB.VerCorrelationLength
			>> props.GSLIB.NumberOfCellsX
			>> props.GSLIB.GridSizeX
			>> props.GSLIB.NumberOfCellsY
			>> props.GSLIB.GridSizeY
			>> props.GSLIB.NumberOfRealizations
			>> props.GSLIB.Seed;
		getline(propsFile, line);
		getline(propsFile, line);
		getline(propsFile, line);
		isDataLine.clear();
		isDataLine.str(line);
		isDataLine
			>> props.GSLIB.isHeterLambda
			>> props.GSLIB.isHeterK
			>> props.GSLIB.isHeterBC
			>> props.GSLIB.isHeterFP;
		getline(propsFile, line);
		getline(propsFile, line);
		getline(propsFile, line);
		isDataLine.clear();
		isDataLine.str(line);
		isDataLine
			>> props.GSLIB.BCHeterogeneityParams.CorLen
			>> props.GSLIB.BCHeterogeneityParams.NumberOfCells
			>> props.GSLIB.BCHeterogeneityParams.deltaX;
		getline(propsFile, line);
		getline(propsFile, line);
		getline(propsFile, line);
		isDataLine.clear();
		isDataLine.str(line);
		isDataLine
			>> props.VInteg.isIntegrate
			>> props.VInteg.xResolution
			>> props.VInteg.nGP;
		props.HInteg.isIntegrate = props.VInteg.isIntegrate;
		props.HInteg.xResolution = props.VInteg.xResolution;
		props.HInteg.nGP = props.VInteg.nGP;
	}

	return props;
}

void BCSgsimParameterFile()
{
	srand(static_cast<int>(time(0)));
	long seedGSLIB = 305 * abs(rand());

	/*for (int i = 0; i < 5; ++i) {
		seedGSLIB = abs(seedGSLIB << 15) | (rand() & 0x7FFF);
	}*/

	FILE* inputFileGSLIB = fopen(GetGSLIBOutputFilePath("BCSgsimInput.par", iParallelRlzn).c_str(), "w");

	fprintf(inputFileGSLIB, "Parameters for SGSIM\n\n");
	fprintf(inputFileGSLIB, "START OF PARAMETERS\n");
	fprintf(inputFileGSLIB, "nodata\n");
	fprintf(inputFileGSLIB, "1 2 0 3 5 0                             -columns for X,Y,Z,vr,wt,sec.var.\n");
	fprintf(inputFileGSLIB, "-1E+21 1E+21\n");
	fprintf(inputFileGSLIB, "0                                       -transform the data (0=no, 1=yes)\n");
	fprintf(inputFileGSLIB, "sgsim.trn\n");
	fprintf(inputFileGSLIB, "0                                       -consider ref. dist (0=no, 1=yes)\n");
	fprintf(inputFileGSLIB, "vmodel.var\n");
	fprintf(inputFileGSLIB, "1 2                                     -columns for vr and wt\n");
	fprintf(inputFileGSLIB, "0 15                                    -zmin,zmax (tail extrapolation)\n");
	fprintf(inputFileGSLIB, "1 0                                     -lower tail option\n");
	fprintf(inputFileGSLIB, "1 15                                    -upper tail option\n");
	fprintf(inputFileGSLIB, "0                                       -debug level (0-3)\n");
	fprintf(inputFileGSLIB, "GSLIB/outputs/R%i_BCSGSIM_nodata.dbg\n", iParallelRlzn);
	fprintf(inputFileGSLIB, "GSLIB/outputs/R%i_BCSGSIM_output.out\n", iParallelRlzn);
	fprintf(inputFileGSLIB, "%i                                       -number of realizations to generate\n", nRealization);
	fprintf(inputFileGSLIB, "%i 0 %e                              -nx, xmin, xsize\n", PROPS.GSLIB.BCHeterogeneityParams.NumberOfCells + 1, PROPS.GSLIB.BCHeterogeneityParams.deltaX);
	fprintf(inputFileGSLIB, "%i 0 %e                              -ny, ymin, ysize\n", 1, 0.0);
	fprintf(inputFileGSLIB, "1 0 1                                   -nz, zmin, zsize\n");
	if (!PROPS.GSLIB.Seed)
	{
		fprintf(inputFileGSLIB, "%ld                               -random number seed\n", (abs(seedGSLIB) + iParallelRlzn));
	}
	else
	{
		fprintf(inputFileGSLIB, "%ld                               -random number seed\n", (PROPS.GSLIB.Seed + iParallelRlzn));
	}

	fprintf(inputFileGSLIB, "0 8                                     -Min and max original data for sim\n");
	fprintf(inputFileGSLIB, "12                                      -number of simulated nodes to use\n");
	fprintf(inputFileGSLIB, "1                                       -assign data to nodes (0=no, 1=yes)\n");
	fprintf(inputFileGSLIB, "1 3                                     -multiple grid search (0=no, 1=yes), num\n");
	fprintf(inputFileGSLIB, "0                                       -maximum data per octant (0=not used)\n");
	fprintf(inputFileGSLIB, "10 10 10                                -maximum search radii (hmax, hmin, vert)\n");
	fprintf(inputFileGSLIB, "0 0 0                                   -angles for search ellipsoid\n");
	fprintf(inputFileGSLIB, "51 51 11                                -size of covariance lookup table\n");
	fprintf(inputFileGSLIB, "0 0 1                                   -kType: 0=SK,1=OK,2=LVM,3=EXDR,4=COLC\n");
	fprintf(inputFileGSLIB, "nodata\n");
	fprintf(inputFileGSLIB, "4                                       -column\n");
	fprintf(inputFileGSLIB, "1 0                                     -nst, nugget NOFILE\n");
	fprintf(inputFileGSLIB, "2 1 0 0 0                               -it, cc, ang1, ang2, ang3\n");
	fprintf(inputFileGSLIB, "%e %e 0								 -a_hmax, a_hmin, a_vert\n", PROPS.GSLIB.BCHeterogeneityParams.CorLen, PROPS.GSLIB.BCHeterogeneityParams.CorLen);
	fflush(inputFileGSLIB);
}

void BCAddcoorParameterFile(int nRealization)
{
	FILE* inputFileGSLIB = fopen(GetGSLIBOutputFilePath("BCAddcoorInput.par", iParallelRlzn).c_str(), "w");

	fprintf(inputFileGSLIB, "Parameters for SGSIM\n\n");
	fprintf(inputFileGSLIB, "START OF PARAMETERS\n");
	fprintf(inputFileGSLIB, "GSLIB/outputs/R%i_BCSGSIM_output.out\n", iParallelRlzn);
	fprintf(inputFileGSLIB, "GSLIB/outputs/R%i_BCADDCOOR_output.out\n", iParallelRlzn);
	fprintf(inputFileGSLIB, "%i                                   -realization number to add coordinate\n", nRealization);
	fprintf(inputFileGSLIB, "%i 0 %e                  -nx, xmin, xsize\n", PROPS.GSLIB.BCHeterogeneityParams.NumberOfCells + 1, PROPS.GSLIB.BCHeterogeneityParams.deltaX);
	fprintf(inputFileGSLIB, "%i 0 %e                  -ny, ymin, ysize\n", 1, 1.0);
	fprintf(inputFileGSLIB, "1 0 1                               -nz, zmin, zsize\n");
	fflush(inputFileGSLIB);
}

void BCGSLIBRunSGSIM()
{
#ifdef _windows_
	string sgsimPath = "GSLIB\\bin\\GSLIBSimulation.exe";
	string sgsimInputPath = "\"GSLIB/outputs/R" + to_string(iParallelRlzn) + "_BCSgsimInput.par\"";
	string sgsimArgStr = sgsimPath + " " + sgsimInputPath;
	const char* sgsimArg = sgsimArgStr.c_str();
	string addcoorPath = "GSLIB\\bin\\GSLIBAddCoordinates.exe";
	string addcoorInputPath = "\"GSLIB/outputs/R" + to_string(iParallelRlzn) + "_BCAddcoorInput.par\"";
	string addcoorArgStr = addcoorPath + " " + addcoorInputPath;
	const char* addCoorArg = addcoorArgStr.c_str();
#else
	string sgsimPath = "GSLIB/bin/GSLIBSimulation";
	string sgsimInputPath = "\"GSLIB/outputs/R" + to_string(iParallelRlzn) + "_BCSgsimInput.par\"";
	string sgsimArgStr = sgsimPath + " " + sgsimInputPath;
	const char* sgsimArg = sgsimArgStr.c_str();
	string addcoorPath = "GSLIB/bin/GSLIBAddCoordinates";
	string addcoorInputPath = "\"GSLIB/outputs/R" + to_string(iParallelRlzn) + "_BCAddcoorInput.par\"";
	string addcoorArgStr = addcoorPath + " " + addcoorInputPath;
	const char* addCoorArg = addcoorArgStr.c_str();
#endif  
	cout << endl << "=== GSLIB UNCONDITIONAL SIMULATION ===" << endl;
	int ExecuteGSLIB;
	try
	{
		ExecuteGSLIB = system(sgsimArg); // sequential gaussian simulation
	}
	catch (int e)
	{
		cout << "GSLIB FAILED RUNNING SEQUENTIAL GAUSSIAN SIMULATION! EXCEPTION NO: " << e << endl;
	}

	cout << "=== END GSLIB UNCONDITIONAL SIMULATION ===" << endl;
	BCGSLIBCoeffs.resize(PROPS.GSLIB.BCHeterogeneityParams.NumberOfCells + 1, nRealization);
	BCGSLIBCoeffs.setZero();
	BCGSLIBGrid.resize(PROPS.GSLIB.BCHeterogeneityParams.NumberOfCells + 1, 2);
	BCGSLIBGrid.setZero();
	for (int iRealization = 0; iRealization < nRealization; iRealization++)
	{
		BCAddcoorParameterFile(iRealization + 1);
		try
		{
			ExecuteGSLIB = system(addCoorArg); // adding coordinates
		}
		catch (int e)
		{
			cout << "GSLIB FAILED ADDING COORDINATES! EXCEPTION NO: " << e << endl;
		}

		FILE* inputFileGSLIB = fopen(GetOutputFilePath("BCGSLIB_Simulation.plt", iParallelRlzn).c_str(), "w");
		fprintf(inputFileGSLIB, "TITLE = 'Boundary Condition Spatial Variability GSLIB Output' \n");
		fprintf(inputFileGSLIB, "VARIABLES = \"<i>x </i>(m)\" \"GSLIB Coeff\" \n");
		fprintf(inputFileGSLIB, "Zone T = \"1\" \n");
		string AddcoorOutputFilePathStr = GetGSLIBOutputFilePath("BCADDCOOR_output.out", iParallelRlzn);
		const char* AddcoorOutputFilePath = AddcoorOutputFilePathStr.c_str();
		//system("chmod 755 GSLIB/ADDCOOR_output.out");
		ifstream gslibFile;
		gslibFile.open(AddcoorOutputFilePath);
		string line;
		getline(gslibFile, line);
		getline(gslibFile, line);
		getline(gslibFile, line);
		getline(gslibFile, line);
		getline(gslibFile, line);
		int index = 0;
		double gslibCoeff;
		double xGrid;
		double yGrid;
		double zCoord;
		while (getline(gslibFile, line))
		{
			if (!(gslibFile >> xGrid >> yGrid >> zCoord >> gslibCoeff))
			{
				break;
			}
			BCGSLIBGrid(index, 0) = xGrid;
			BCGSLIBGrid(index, 1) = yGrid;
			BCGSLIBCoeffs(index, iRealization) = gslibCoeff;
			fprintf(inputFileGSLIB, "%e\t%e\n", xGrid, gslibCoeff);
			index++;
		}

		/*double maxGslibCoeff = BCGSLIBCoeffs.col(iRealization).maxCoeff();
		double minGslibCoeff = BCGSLIBCoeffs.col(iRealization).minCoeff();
		double rangeGslibCoeff = maxGslibCoeff - minGslibCoeff;
		for (int i = 0; i < BCGSLIBCoeffs.size(); i++)
		{
			BCGSLIBCoeffs(i, iRealization) -= minGslibCoeff;
			BCGSLIBCoeffs(i, iRealization) = BCGSLIBCoeffs(i, iRealization) / rangeGslibCoeff;
			BCGSLIBCoeffs(i, iRealization) = -1.0 + 2 * BCGSLIBCoeffs(i, iRealization);
		}*/

		IBConditions BC(&MESH, &PROPS);
		fprintf(inputFileGSLIB, "TITLE = 'Normalized Nodal Spatial Variability GSLIB Output' \n");
		fprintf(inputFileGSLIB, "VARIABLES = \"<i>x </i>(m)\" \"GSLIB Coeff\" \n");
		fprintf(inputFileGSLIB, "Zone T = \"2\" \n");
		PROPS.GSLIB.BCHeterogeneityParams.coeffs.resize(BC.TopBC.size(), 2);
		PROPS.GSLIB.BCHeterogeneityParams.coeffs.setZero();
		for (int n = 0; n < BC.TopBC.size(); n++)
		{
			double xNode = MESH.Nodes[BC.TopBC[n]].Coordinates.x;
			double distanceP = INFINITY;
			for (int i = 0; i < BCGSLIBCoeffs.size(); i++)
			{
				double xGrid = BCGSLIBGrid(i, 0);
				double distance = abs(xNode - xGrid);
				if (distance < distanceP)
				{
					distanceP = distance;
					gslibCoeff = BCGSLIBCoeffs(i, iRealization);
					PROPS.GSLIB.BCHeterogeneityParams.coeffs(n, 0) = BC.TopBC[n];
					PROPS.GSLIB.BCHeterogeneityParams.coeffs(n, 1) = gslibCoeff;
				}
			}

			fprintf(inputFileGSLIB, "%e\t%e\n", xNode, PROPS.GSLIB.BCHeterogeneityParams.coeffs(n, 1));
		}

		fflush(inputFileGSLIB);
	}
}

void SgsimParameterFile()
{
	srand(static_cast<int>(time(0)));
	long seedGSLIB = 305 * abs(rand());

	/*for (int i = 0; i < 5; ++i) {
		seedGSLIB = abs(seedGSLIB << 15) | (rand() & 0x7FFF);
	}*/

	FILE *inputFileGSLIB = fopen(GetGSLIBOutputFilePath("SgsimInput.par", iParallelRlzn).c_str(), "w");

	fprintf(inputFileGSLIB, "Parameters for SGSIM\n\n");
	fprintf(inputFileGSLIB, "START OF PARAMETERS\n");
	fprintf(inputFileGSLIB, "nodata\n");
	fprintf(inputFileGSLIB, "1 2 0 3 5 0                             -columns for X,Y,Z,vr,wt,sec.var.\n");
	fprintf(inputFileGSLIB, "-1E+21 1E+21\n");
	fprintf(inputFileGSLIB, "0                                       -transform the data (0=no, 1=yes)\n");
	fprintf(inputFileGSLIB, "sgsim.trn\n");
	fprintf(inputFileGSLIB, "0                                       -consider ref. dist (0=no, 1=yes)\n");
	fprintf(inputFileGSLIB, "vmodel.var\n");
	fprintf(inputFileGSLIB, "1 2                                     -columns for vr and wt\n");
	fprintf(inputFileGSLIB, "0 15                                    -zmin,zmax (tail extrapolation)\n");
	fprintf(inputFileGSLIB, "1 0                                     -lower tail option\n");
	fprintf(inputFileGSLIB, "1 15                                    -upper tail option\n");
	fprintf(inputFileGSLIB, "0                                       -debug level (0-3)\n");
	fprintf(inputFileGSLIB, "GSLIB/outputs/R%i_SGSIM_nodata.dbg\n", iParallelRlzn);
	fprintf(inputFileGSLIB, "GSLIB/outputs/R%i_SGSIM_output.out\n", iParallelRlzn);
	fprintf(inputFileGSLIB, "%i                                       -number of realizations to generate\n", nRealization);
	fprintf(inputFileGSLIB, "%i 0 %e                              -nx, xmin, xsize\n", PROPS.GSLIB.NumberOfCellsX + 1, PROPS.GSLIB.GridSizeX);
	fprintf(inputFileGSLIB, "%i 0 %e                              -ny, ymin, ysize\n", PROPS.GSLIB.NumberOfCellsY + 1, PROPS.GSLIB.GridSizeY);
	fprintf(inputFileGSLIB, "1 0 1                                   -nz, zmin, zsize\n");
	if (!PROPS.GSLIB.Seed)
	{
		fprintf(inputFileGSLIB, "%ld                               -random number seed\n", (abs(seedGSLIB) + iParallelRlzn));
	}
	else
	{
		fprintf(inputFileGSLIB, "%ld                               -random number seed\n", (PROPS.GSLIB.Seed + iParallelRlzn));
	}

	fprintf(inputFileGSLIB, "0 8                                     -Min and max original data for sim\n");
	fprintf(inputFileGSLIB, "12                                      -number of simulated nodes to use\n");
	fprintf(inputFileGSLIB, "1                                       -assign data to nodes (0=no, 1=yes)\n");
	fprintf(inputFileGSLIB, "1 3                                     -multiple grid search (0=no, 1=yes), num\n");
	fprintf(inputFileGSLIB, "0                                       -maximum data per octant (0=not used)\n");
	fprintf(inputFileGSLIB, "10 10 10                                -maximum search radii (hmax, hmin, vert)\n");
	fprintf(inputFileGSLIB, "0 0 0                                   -angles for search ellipsoid\n");
	fprintf(inputFileGSLIB, "51 51 11                                -size of covariance lookup table\n");
	fprintf(inputFileGSLIB, "0 0 1                                   -kType: 0=SK,1=OK,2=LVM,3=EXDR,4=COLC\n");
	fprintf(inputFileGSLIB, "nodata\n");
	fprintf(inputFileGSLIB, "4                                       -column\n");
	fprintf(inputFileGSLIB, "1 0                                     -nst, nugget NOFILE\n");
	fprintf(inputFileGSLIB, "2 1 0 0 0                               -it, cc, ang1, ang2, ang3\n");
	fprintf(inputFileGSLIB, "%e %e 0                           -a_hmax, a_hmin, a_vert\n", PROPS.GSLIB.VerCorrelationLength, PROPS.GSLIB.HorCorrelationLength);
	fflush(inputFileGSLIB);
}

void AddcoorParameterFile(int nRealization)
{
	FILE *inputFileGSLIB = fopen(GetGSLIBOutputFilePath("AddcoorInput.par", iParallelRlzn).c_str(), "w");

	fprintf(inputFileGSLIB, "Parameters for SGSIM\n\n");
	fprintf(inputFileGSLIB, "START OF PARAMETERS\n");
	fprintf(inputFileGSLIB, "GSLIB/outputs/R%i_SGSIM_output.out\n", iParallelRlzn);
	fprintf(inputFileGSLIB, "GSLIB/outputs/R%i_ADDCOOR_output.out\n", iParallelRlzn);
	fprintf(inputFileGSLIB, "%i                                   -realization number to add coordinate\n", nRealization);
	fprintf(inputFileGSLIB, "%i 0 %e                  -nx, xmin, xsize\n", PROPS.GSLIB.NumberOfCellsX + 1, PROPS.GSLIB.GridSizeX);
	fprintf(inputFileGSLIB, "%i 0 %e                  -ny, ymin, ysize\n", PROPS.GSLIB.NumberOfCellsY + 1, PROPS.GSLIB.GridSizeY);
	fprintf(inputFileGSLIB, "1 0 1                               -nz, zmin, zsize\n");
	fflush(inputFileGSLIB);
}

void GSLIBRunSGSIM()
{
	#ifdef _windows_
	string sgsimPath = "GSLIB\\bin\\GSLIBSimulation.exe";
	string sgsimInputPath = "\"GSLIB/outputs/R"+ to_string(iParallelRlzn) +"_SgsimInput.par\"";
	string sgsimArgStr = sgsimPath + " " + sgsimInputPath;
	const char* sgsimArg = sgsimArgStr.c_str();
	string addcoorPath = "GSLIB\\bin\\GSLIBAddCoordinates.exe";
	string addcoorInputPath = "\"GSLIB/outputs/R" + to_string(iParallelRlzn) + "_AddcoorInput.par\"";
	string addcoorArgStr = addcoorPath + " " + addcoorInputPath;
	const char* addCoorArg = addcoorArgStr.c_str();
	#else
	string sgsimPath = "GSLIB/bin/GSLIBSimulation";
	string sgsimInputPath = "\"GSLIB/outputs/R" + to_string(iParallelRlzn) + "_SgsimInput.par\"";
	string sgsimArgStr = sgsimPath + " " + sgsimInputPath;
	const char* sgsimArg = sgsimArgStr.c_str();
	string addcoorPath = "GSLIB/bin/GSLIBAddCoordinates";
	string addcoorInputPath = "\"GSLIB/outputs/R" + to_string(iParallelRlzn) + "_AddcoorInput.par\"";
	string addcoorArgStr = addcoorPath + " " + addcoorInputPath;
	const char* addCoorArg = addcoorArgStr.c_str();
	#endif  
	cout << endl << "=== GSLIB UNCONDITIONAL SIMULATION ===" << endl;
	int ExecuteGSLIB;
	try
	{
		ExecuteGSLIB = system(sgsimArg); // sequential gaussian simulation
	}
	catch (int e)
	{
		cout << "GSLIB FAILED RUNNING SEQUENTIAL GAUSSIAN SIMULATION! EXCEPTION NO: "<< e << endl;
	}
	
	cout << "=== END GSLIB UNCONDITIONAL SIMULATION ===" << endl;
	GSLIBCoeffs.resize((PROPS.GSLIB.NumberOfCellsX + 1) * (PROPS.GSLIB.NumberOfCellsY + 1), nRealization);
	GSLIBCoeffs.setZero();
	GSLIBGrid.resize((PROPS.GSLIB.NumberOfCellsX + 1) * (PROPS.GSLIB.NumberOfCellsY + 1), 2);
	GSLIBGrid.setZero();
	for (int iRealization = 0; iRealization < nRealization; iRealization++)
	{
		AddcoorParameterFile(iRealization + 1);
		try
		{
			ExecuteGSLIB = system(addCoorArg); // adding coordinates
		}
		catch (int e)
		{
			cout << "GSLIB FAILED ADDING COORDINATES! EXCEPTION NO: " << e << endl;
		}

		FILE *inputFileGSLIB = fopen(GetOutputFilePath("GSLIB_Simulation.plt", iParallelRlzn).c_str(), "w");
		string AddcoorOutputFilePathStr = GetGSLIBOutputFilePath("ADDCOOR_output.out", iParallelRlzn);
		const char* AddcoorOutputFilePath = AddcoorOutputFilePathStr.c_str();
		//system("chmod 755 GSLIB/ADDCOOR_output.out");
		ifstream gslibFile;
		gslibFile.open(AddcoorOutputFilePath);
		string line;
		getline(gslibFile, line);
		getline(gslibFile, line);
		getline(gslibFile, line);
		getline(gslibFile, line);
		getline(gslibFile, line);
		int index = 0;
		double gslibCoeff;
		double xGrid;
		double yGrid;
		double zCoord;
		while (getline(gslibFile, line))
		{
			if (!(gslibFile >> xGrid >> yGrid >> zCoord >> gslibCoeff))
			{
				break;
			}
			GSLIBGrid(index, 0) = xGrid;
			GSLIBGrid(index, 1) = yGrid;
			GSLIBCoeffs(index, iRealization) = gslibCoeff;
			fprintf(inputFileGSLIB, "%e\t%e\t%e\n", xGrid, yGrid, gslibCoeff);
			index++;
		}

		fflush(inputFileGSLIB);
	}
}

void UpscaleGSLIBtoSIMPER()
{
	if (PROPS.GSLIB.isGSLIB)
	{
		double GSLIBCoeffE;
		VectorXi elementDofs;
		VectorXd GSLIBCoeffsNode;
		FILE* plotSoilProperties = fopen(GetOutputFilePath("SoilProperties.plt", iParallelRlzn).c_str(), "w");
		SoilEmpiricalRelations SER;
		NodalGSLIBCoeffs.resize(nond, nRealization);
		ElementalGSLIBCoeffs.resize(noel, nRealization);
		NodalGSLIBCoeffs.setZero();
		int gslibNumberOfNodes = (PROPS.GSLIB.NumberOfCellsX + 1) * (PROPS.GSLIB.NumberOfCellsY + 1);
		for (int iReal = 0; iReal < nRealization; iReal++)
		{
			for (int n = 0; n < nond; n++)
			{
				double xNode = MESH.Nodes[n].Coordinates.x;
				double yNode = MESH.Nodes[n].Coordinates.y;
				double gslibCoeff;
				double distanceP = INFINITY;

				for (int i = 0; i < gslibNumberOfNodes; i++)
				{
					double xGrid = GSLIBGrid(i, 0);
					double yGrid = GSLIBGrid(i, 1);
					double distance = sqrt((xNode - xGrid) * (xNode - xGrid) + (yNode - yGrid) * (yNode - yGrid));
					if (distance < distanceP)
					{
						distanceP = distance;
						gslibCoeff = GSLIBCoeffs(i, iReal);
						NodalGSLIBCoeffs(n, iReal) = gslibCoeff;
					}
				}
			}

			for (int e = 0;e < noel;e++)
			{
				elementDofs = MESH.GetElementDofs(e, ndoe);
				GSLIBCoeffsNode = MESH.GetNodalValues(NodalGSLIBCoeffs.col(iReal), elementDofs);
				ElementalGSLIBCoeffs(e, iReal) = 0;
				for (int ig = 0; ig < 4; ig++)
				{
					ElementalGSLIBCoeffs(e, iReal) += GSLIBCoeffsNode(ig);
				}

				ElementalGSLIBCoeffs(e, iReal) = ElementalGSLIBCoeffs(e, iReal) * 0.25;
			}

			double maxGslibCoeff = ElementalGSLIBCoeffs.col(iReal).maxCoeff();
			double minGslibCoeff = ElementalGSLIBCoeffs.col(iReal).minCoeff();
			double rangeGslibCoeff = maxGslibCoeff - minGslibCoeff;
			for (int i = 0; i < noel; i++)
			{
				ElementalGSLIBCoeffs(i, iReal) -= minGslibCoeff;
				ElementalGSLIBCoeffs(i, iReal) = ElementalGSLIBCoeffs(i, iReal) / rangeGslibCoeff;
			}
		}

		for (int iReal = 0; iReal < nRealization; iReal++)
		{
			for (int e = 0; e < noel; e++)
			{
				int iSoilType = MESH.Elements[e].iSoilType;
				GSLIBCoeffE = ElementalGSLIBCoeffs(e, iReal);
				if (PROPS.GSLIB.isHeterBC || PROPS.GSLIB.isHeterFP || PROPS.GSLIB.isHeterK || PROPS.GSLIB.isHeterLambda)
				{
					MESH.Elements[e].SoilDensity = PROPS.Soil[iSoilType].DensityMin + (PROPS.Soil[iSoilType].DensityMax - PROPS.Soil[iSoilType].DensityMin) * GSLIBCoeffE;
				}
				else
				{
					MESH.Elements[e].SoilDensity = 0.5 * (PROPS.Soil[iSoilType].DensityMin + PROPS.Soil[iSoilType].DensityMax);
				}

				MESH.Elements[e].SoilHeatCapacity = PROPS.Soil[iSoilType].HeatCapacity;
				MESH.Elements[e].SoilThermalConductivity = SER.ThermalConductivity(MESH.Elements[e].SoilDensity);

				if (MESH.Elements[e].SoilType == "Fen" || MESH.Elements[e].SoilType == "Bog")
				{
					MESH.Elements[e].SoilHydraulicConductivity = PROPS.Soil[iSoilType].HydraulicConductivity;
				}
				else
				{
					MESH.Elements[e].SoilHydraulicConductivity = SER.HydraulicConductivity(MESH.Elements[e].SoilDensity, PROPS.Soil[iSoilType].DensityMax, PROPS.Soil[iSoilType].DensityMin);
				}

				if (PROPS.GSLIB.isHeterFP) // Check if soil freezing point is heterogeneous
				{
					//MESH.Elements[e].SoilFreezingPoint = PROPS.Nonisothermal.TempLiquid - (PROPS.Nonisothermal.TempLiquid - PROPS.Nonisothermal.TempSolid) * (1.0 + GSLIBCoeffE);
					MESH.Elements[e].SoilFreezingPoint = PROPS.Soil[iSoilType].FPmin - (abs(PROPS.Soil[iSoilType].FPmax - PROPS.Soil[iSoilType].FPmin)) * GSLIBCoeffE;
				}
				else
				{
					MESH.Elements[e].SoilFreezingPoint = PROPS.Nonisothermal.TempSolid;
				}
			}
		}

		for (int e = 0; e < noel; e++)
		{
			fprintf(plotSoilProperties, "variables =\"X\" \"Y\"");
			fprintf(plotSoilProperties, " \"Realization %i <greek>l</greek><sub>soil</sub>\"", iParallelRlzn);
			fprintf(plotSoilProperties, " \"Realization %i <i>c</i><sub>soil</sub>\"", iParallelRlzn);
			fprintf(plotSoilProperties, " \"Realization %i <i>K</i><sub>soil</sub>\"", iParallelRlzn);
			fprintf(plotSoilProperties, " \"Realization %i <i>D</i><sub>soil</sub>\"", iParallelRlzn);
			fprintf(plotSoilProperties, " \"Realization %i Freezing Point\"", iParallelRlzn);
			fprintf(plotSoilProperties, " \"Realization %i Coeff_GSLIB\"", iParallelRlzn);
			fprintf(plotSoilProperties, " \"Realization %i Soil type\"", iParallelRlzn);


			fprintf(plotSoilProperties, "\nZONE N = %5.0d, E = %5.0d, ZONETYPE = FEQuadrilateral, DATAPACKING = POINT\n", ndoe, 1);

			for (int inod = 0; inod < ndoe; inod++)
			{
				fprintf(plotSoilProperties, "%e\t%e\t", MESH.Elements[e].Nodes[inod].Coordinates.x,
					MESH.Elements[e].Nodes[inod].Coordinates.y);
				fprintf(plotSoilProperties, "%e\t", MESH.Elements[e].SoilHydraulicConductivity);
				fprintf(plotSoilProperties, "%e\t", MESH.Elements[e].SoilHeatCapacity);
				fprintf(plotSoilProperties, "%e\t", MESH.Elements[e].SoilThermalConductivity);
				fprintf(plotSoilProperties, "%e\t", MESH.Elements[e].SoilDensity);
				fprintf(plotSoilProperties, "%e\t", MESH.Elements[e].SoilFreezingPoint);
				fprintf(plotSoilProperties, "%e\t", ElementalGSLIBCoeffs(e, 0));
				fprintf(plotSoilProperties, "%i\t", MESH.Elements[e].iSoilType);
				fprintf(plotSoilProperties, "\n");
			}

			fprintf(plotSoilProperties, "1 2 3 4\n");
			fflush(plotSoilProperties);
		}
	}
	else
	{
		cout << "Homogeneous Media" << endl;
		FILE* plotSoilProperties = fopen(GetOutputFilePath("SoilProperties.plt", iParallelRlzn).c_str(), "w");;
		for (int e = 0; e < MESH.NumberOfElements; e++)
		{
			int iSoilType = MESH.Elements[e].iSoilType;
			MESH.Elements[e].SoilHeatCapacity = PROPS.Soil[iSoilType].HeatCapacity;
			MESH.Elements[e].SoilThermalConductivity = PROPS.Soil[iSoilType].ThermalConductivity;
			MESH.Elements[e].SoilDensity = 0.5 * (PROPS.Soil[iSoilType].DensityMax + PROPS.Soil[iSoilType].DensityMin);
			MESH.Elements[e].SoilFreezingPoint = PROPS.Nonisothermal.TempSolid;
			MESH.Elements[e].SoilHydraulicConductivity = PROPS.Soil[iSoilType].HydraulicConductivity;


			fprintf(plotSoilProperties, "variables =\"X\" \"Y\" \"<greek>l</greek><sub>soil</sub>\" \"<i>c</i><sub>soil</sub>\" \"<i>K</i><sub>soil</sub>\" \"<i>D</i><sub>soil</sub>\" \"Freezing Point\" \"Soil type\"\n");
			fprintf(plotSoilProperties, "ZONE N = %5.0d, E = %5.0d, ZONETYPE = FEQuadrilateral, DATAPACKING = POINT\n", ndoe, 1);
			for (int inod = 0; inod < ndoe; inod++)
			{
				int nodeIndex = MESH.Elements[e].Nodes[inod].n - 1;
				fprintf(plotSoilProperties, "%e\t%e\t%e\t%e\t%e\t%e\t%e\t%i\n", MESH.Elements[e].Nodes[inod].Coordinates.x,
					MESH.Elements[e].Nodes[inod].Coordinates.y,
					MESH.Elements[e].SoilHydraulicConductivity,
					MESH.Elements[e].SoilHeatCapacity,
					MESH.Elements[e].SoilThermalConductivity,
					MESH.Elements[e].SoilDensity,
					MESH.Elements[e].SoilFreezingPoint,
					MESH.Elements[e].iSoilType);
			}

			fprintf(plotSoilProperties, "1 2 3 4\n");
			fflush(plotSoilProperties);
		}
	}
}

Mesh InputMesh(string filePath)
{
	Mesh mesh;
	Node node;
	Element element;

	double *maxX = &mesh.MaxX;
	double *minX = &mesh.MinX;
	double *maxY = &mesh.MaxY;
	double *minY = &mesh.MinY;
	double *maxZ = &mesh.MaxZ;
	double *minZ = &mesh.MinZ;

	*maxX = -1.0* INFINITY;
	*minX = INFINITY;
	*maxY = -1.0 * INFINITY;
	*minY = INFINITY;
	*maxZ = -1.0 * INFINITY;
	*minZ = INFINITY;

	if (PROPS.Solution.IsGMSH)
	{
		ifstream meshFile;
		meshFile.open(filePath.c_str());
		if (meshFile.fail())
		{
			cout << "Cannot find file " << filePath << endl;
		}
		else
		{
			string line;

			getline(meshFile, line);
			getline(meshFile, line);
			getline(meshFile, line);
			getline(meshFile, line);
			getline(meshFile, line);

			istringstream isSurfaceNumber(line);
			int nSurface;
			vector<int> surfaceCode;
			int dimensions;
			vector<string> SoilType;

			isSurfaceNumber >> nSurface;
			surfaceCode.resize(nSurface);
			SoilType.resize(nSurface);

			for (int i = 0; i < nSurface; i++)
			{
				getline(meshFile, line);
				istringstream surfaceCodeStr(line);
				if (surfaceCodeStr >> dimensions >> surfaceCode[i] >> SoilType[i])
				{
					surfaceCodeStr >> dimensions >> surfaceCode[i] >> SoilType[i];
					SoilType[i].erase(
						remove(SoilType[i].begin(), SoilType[i].end(), '\"'),
						SoilType[i].end()
					);
				}
			}

			getline(meshFile, line);
			getline(meshFile, line);
			getline(meshFile, line);

			// Nodes Data
			while (getline(meshFile, line))
			{
				istringstream isNodeData(line);
				int n;
				double x, y, z;

				if (!(isNodeData >> n >> x >> y >> z))
				{
					break;
				}

				node.n = n;
				node.Coordinates.x = x;
				node.Coordinates.y = y;
				node.Coordinates.z = z;

				*minX = (x < *minX) ? x : *minX;
				*maxX = (x > *maxX) ? x : *maxX;
				*minY = (y < *minY) ? y : *minY;
				*maxY = (y > *maxY) ? y : *maxY;
				*minZ = (z < *minZ) ? z : *minZ;
				*maxZ = (z > *maxZ) ? z : *maxZ;

				mesh.Nodes.push_back(node);
			}

			getline(meshFile, line);
			getline(meshFile, line);

			//Elements Data
			while (getline(meshFile, line))
			{
				istringstream isElementData(line);
				int e, X1, X2, isoiltype, X4, n1, n2, n3, n4;

				if (!(isElementData >> e >> X1 >> X2 >> isoiltype >> X4 >> n1 >> n2 >> n3 >> n4))
				{
					break;
				}

				element.Nodes.clear();
				element.iSoilType = isoiltype - 1;
				element.SoilType = SoilType[isoiltype - 1];
				if (element.SoilType == "Fen")
				{
					GlobalFenIndex = true;
				}

				element.e = e;
				node.n = n1;
				node.Coordinates = mesh.Nodes[n1 - 1].Coordinates;
				element.Nodes.push_back(node);
				node.n = n2;
				node.Coordinates = mesh.Nodes[n2 - 1].Coordinates;
				element.Nodes.push_back(node);
				node.n = n3;
				node.Coordinates = mesh.Nodes[n3 - 1].Coordinates;
				element.Nodes.push_back(node);
				node.n = n4;
				node.Coordinates = mesh.Nodes[n4 - 1].Coordinates;
				element.Nodes.push_back(node);

				mesh.Elements.push_back(element);
			}

			mesh.NumberOfElements = static_cast<int>(mesh.Elements.size());
			mesh.NumberOfNodes = static_cast<int>(mesh.Nodes.size());
			mesh.ElementNumberOfNodes = 4;
			noel = mesh.NumberOfElements;
			nond = mesh.NumberOfNodes;
			ndoe = mesh.ElementNumberOfNodes;
		}
	}
	else
	{
		ifstream meshFile;
		meshFile.open(filePath.c_str());
		if (meshFile.fail())
		{
			cout << "Cannot find file " << filePath << endl;
		}
		else
		{
			string line;

			getline(meshFile, line);
			getline(meshFile, line);

			// Nodes Data
			while (getline(meshFile, line))
			{
				istringstream isNodeData(line);
				int n;
				double x, y, z;

				if (!(isNodeData >> n >> x >> y >> z))
				{
					break;
				}

				node.n = n;
				node.Coordinates.x = x;
				node.Coordinates.y = y;
				node.Coordinates.z = z;

				*minX = (x < *minX) ? x : *minX;
				*maxX = (x > *maxX) ? x : *maxX;
				*minY = (y < *minY) ? y : *minY;
				*maxY = (y > *maxY) ? y : *maxY;
				*minZ = (z < *minZ) ? z : *minZ;
				*maxZ = (z > *maxZ) ? z : *maxZ;

				mesh.Nodes.push_back(node);
			}

			getline(meshFile, line);
			getline(meshFile, line);

			//Elements Data
			while (getline(meshFile, line))
			{
				istringstream isElementData(line);
				int e, n1, n2, n3, n4;

				if (!(isElementData >> e >> n1 >> n2 >> n3 >> n4))
				{
					break;
				}

				element.Nodes.clear();
				element.e = e;
				node.n = n1;
				node.Coordinates = mesh.Nodes[n1 - 1].Coordinates;
				element.Nodes.push_back(node);
				node.n = n2;
				node.Coordinates = mesh.Nodes[n2 - 1].Coordinates;
				element.Nodes.push_back(node);
				node.n = n3;
				node.Coordinates = mesh.Nodes[n3 - 1].Coordinates;
				element.Nodes.push_back(node);
				node.n = n4;
				node.Coordinates = mesh.Nodes[n4 - 1].Coordinates;
				element.Nodes.push_back(node);

				mesh.Elements.push_back(element);
			}

			mesh.NumberOfElements = static_cast<int>(mesh.Elements.size());
			mesh.NumberOfNodes = static_cast<int>(mesh.Nodes.size());
			mesh.ElementNumberOfNodes = 4;
			noel = mesh.NumberOfElements;
			nond = mesh.NumberOfNodes;
			ndoe = mesh.ElementNumberOfNodes;
		}
	}

	for (int e = 0; e < noel; e++)
	{
		mesh.Elements[e].Area = mesh.GetElementArea(e, ndoe);
	}

	return mesh;
}

//Point ParentCoordinates(vector<Node> nodes, Point p)
//{
//	double xBar = 0.0;
//	double yBar = 0.0;
//	for (int i = 0; i < 4; i++)
//	{
//		xBar += nodes[i].Coordinates.x;
//		yBar += nodes[i].Coordinates.y;
//	}
//
//	xBar = 0.25 * xBar;
//	yBar = 0.25 * yBar;
//
//	double l = sqrt((nodes[1].Coordinates.x - nodes[0].Coordinates.x) * (nodes[1].Coordinates.x - nodes[0].Coordinates.x) +
//		            (nodes[1].Coordinates.y - nodes[0].Coordinates.y) * (nodes[1].Coordinates.y - nodes[0].Coordinates.y));
//
//	double w = sqrt((nodes[2].Coordinates.x - nodes[1].Coordinates.x) * (nodes[2].Coordinates.x - nodes[1].Coordinates.x) +
//		            (nodes[2].Coordinates.y - nodes[1].Coordinates.y) * (nodes[2].Coordinates.y - nodes[1].Coordinates.y));
//	Point point;
//	point.x = 2.0 * (p.x - xBar) / l;
//	point.y = 2.0 * (p.y - yBar) / w;
//
//	return point;
//}

void VerticalIntegrationInfo()
{
	FEMFunctions femFunctions;
	VectorXd xNodes;
	VectorXd yNodes;
	double domainWidth = abs(MESH.MaxX - MESH.MinX);
	double domainDepth = abs(MESH.MaxY - MESH.MinY);
	double deltaX = domainWidth / (double)PROPS.VInteg.xResolution;
	PROPS.VInteg.DomainDepth = domainDepth;
	PROPS.VInteg.DomainWidth = domainWidth;
	for (int i = 0; i <= PROPS.VInteg.xResolution; i++)
	{
		IntGlobalInfo GI;
		IntLocalInfo LI;
		VectorXd SF;
		VectorXd xyGlob;
		double xGlob = MESH.MinX + (double)(i) * deltaX;
		GI.xGlob = xGlob;
		GI.IntSwat = 0.0;
		GI.IntSice = 0.0;
		VIntGP.GP(PROPS.VInteg.nGP);
		GP.GP(nGP);
		for (int iGP = 0; iGP < PROPS.VInteg.nGP; iGP++)
		{
			LI.GPWeight = VIntGP.Weights[iGP];
			double GPi = VIntGP.Points[iGP];
			double yGlob = 0.5 * domainDepth * GPi + 0.5 * domainDepth;
			double distance = 9999.0;
			int iNodeNearest = 1;
			for (int e = 0; e < noel; e++)
			{
				xNodes = MESH.GetNodesXCoordinates(e, ndoe);
				yNodes = MESH.GetNodesYCoordinates(e, ndoe);
				for (int i = 0; i < nGP; i++)
				{
					for (int j = 0; j < nGP; j++)
					{
						double GPiloc = GP.Points[i];
						double GPjloc = GP.Points[j];
						femFunctions.Calculate(xNodes, yNodes, GPiloc, GPjloc);
						SF = femFunctions.SF;
						xyGlob = femFunctions.CalculateGlobalCoordinates(SF, xNodes, yNodes);

						double distance2 = sqrt((xGlob - xyGlob(0)) * (xGlob - xyGlob(0)) +
												(yGlob - xyGlob(1)) * (yGlob - xyGlob(1)));

						if (distance2 < distance)
						{
							distance = distance2;
							LI.iElement = e;
							LI.iGP = GPiloc;
							LI.jGP = GPjloc;
						}
					}
				}
			}
		
			GI.LocalInfo.push_back(LI);
		}

		PROPS.VInteg.GlobalInfo.push_back(GI);
	}
}

void HorizontalIntegrationInfo()
{
	FEMFunctions femFunctions;
	VectorXd xNodes;
	VectorXd yNodes;
	double domainWidth = abs(MESH.MaxX - MESH.MinX);
	double domainDepth = abs(MESH.MaxY - MESH.MinY);
	double deltaY = domainDepth / (double)PROPS.HInteg.xResolution;
	PROPS.HInteg.DomainDepth = domainDepth;
	PROPS.HInteg.DomainWidth = domainWidth;
	for (int i = 0; i <= PROPS.VInteg.xResolution; i++)
	{
		IntGlobalInfo GI;
		IntLocalInfo LI;
		VectorXd SF;
		VectorXd xyGlob;
		double yGlob = MESH.MinY + (double)(i)*deltaY;
		GI.yGlob = yGlob;
		GI.IntSwat = 0.0;
		GI.IntSice = 0.0;
		HIntGP.GP(PROPS.HInteg.nGP);
		GP.GP(nGP);
		for (int iGP = 0; iGP < PROPS.HInteg.nGP; iGP++)
		{
			LI.GPWeight = HIntGP.Weights[iGP];
			double GPi = HIntGP.Points[iGP];
			double xGlob = 0.5 * domainWidth + 0.5 * domainWidth * GPi;
			double distance = 9999.0;
			int iNodeNearest = 1;
			for (int e = 0; e < noel; e++)
			{
				xNodes = MESH.GetNodesXCoordinates(e, ndoe);
				yNodes = MESH.GetNodesYCoordinates(e, ndoe);
				for (int i = 0; i < nGP; i++)
				{
					for (int j = 0; j < nGP; j++)
					{
						double GPiloc = GP.Points[i];
						double GPjloc = GP.Points[j];
						femFunctions.Calculate(xNodes, yNodes, GPiloc, GPjloc);
						SF = femFunctions.SF;
						xyGlob = femFunctions.CalculateGlobalCoordinates(SF, xNodes, yNodes);

						double distance2 = sqrt((xGlob - xyGlob(0)) * (xGlob - xyGlob(0)) +
							(yGlob - xyGlob(1)) * (yGlob - xyGlob(1)));

						if (distance2 < distance)
						{
							distance = distance2;
							LI.iElement = e;
							LI.iGP = GPiloc;
							LI.jGP = GPjloc;
						}
					}
				}
			}

			GI.LocalInfo.push_back(LI);
		}

		PROPS.HInteg.GlobalInfo.push_back(GI);
	}
}

