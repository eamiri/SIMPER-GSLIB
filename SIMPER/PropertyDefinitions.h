using namespace std;

struct SoilProperties
{
	int iType;
	string SoilType;
	double DensityMin;
	double DensityMax;
	double Porosity;
	double ResidualWaterSaturation;
	double HeatCapacity;
	double ThermalConductivity;
	double HydraulicConductivity;
	double Wpar;
	double Mpar;
	int rSFC;
	bool IsSaturated;
	double FPmin;
	double FPmax;
};

struct FluidProperties
{
	double LDensity;
	double LHeatCapacity;
	double LThermalConductivity;
	double SDensity;
	double SHeatCapacity;
	double SThermalConductivity;
	double LatentHeat;
	double ADensity;
	double AHeatCapacity;
	double AThermalConductivity;
};

struct NonisothermalProperties
{
	bool IsIsothermal;
	double HeatCapacity;
	double ThermalConductivity;
	double TempTransition;
	double TempLiquid;
	double TempSolid;
	double SolidSatIndex;
	double LiquidSatIndex;
};

struct BoundaryConditions
{
	bool isBCInput;
	string BCInputFile;
	bool isICInput;
	string ICInputFile;
	double UniformIC;
	vector<int> Node;
	vector<double> Value;
};

struct GaussQuadrature
{
	int NumberOfPoints;
};

struct BCHeterogeneity
{
	double CorLen;
	int NumberOfCells;
	double deltaX;
	MatrixXd coeffs;
};

struct GSLIBProps
{
	bool isGSLIB;
	bool isHeterLambda;
	bool isHeterK;
	bool isHeterBC;
	BCHeterogeneity BCHeterogeneityParams;
	bool isHeterFP;
	double HorCorrelationLength;
	double VerCorrelationLength;
	int NumberOfCellsX;
	int NumberOfCellsY;
	double GridSizeX;
	double GridSizeY;
	double AnisotropyRatio;
	int NumberOfRealizations;
	long Seed;
};

struct SolutionVariables
{
	int MaxTimestep;
	double DeltaTime;
	double TolPsi;
	int MaxIterations;
	double NewmarkGamma;
	int PlotInterval;
	bool IsGMSH;
};

struct IntLocalInfo
{
	int iElement;
	double iGP;
	double jGP;
	double GPWeight;
};

struct IntGlobalInfo
{
	double xGlob;
	double yGlob;
	double IntSwat;
	double IntSice;
	vector<IntLocalInfo> LocalInfo;
};

struct Integration
{
	bool isIntegrate;
	int xResolution;
	int yResolution;
	int nGP;
	double DomainDepth;
	double DomainWidth;
	vector<IntGlobalInfo> GlobalInfo;
};



