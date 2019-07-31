using namespace std;

struct SoilProperties
{
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
	bool IsGSLIB;
	bool IsSaturated;
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

struct FenProperties
{
	bool IsFen;
	double FenDelTemp;
	double FenLambdaRatio;
};

struct BoundaryConditions
{
	string InputFile;
	vector<int> Node;
	vector<double> Value;
};

struct GaussQuadrature
{
	int NumberOfPoints;
};

struct GSLIBProps
{
	bool isHeterLambda;
	bool isHeterK;
	bool isHeterBC;
	bool isHeterFP;
	double CorrelationLength;
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
	bool IsInputBC;
};