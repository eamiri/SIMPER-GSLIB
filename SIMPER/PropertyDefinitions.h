using namespace std;

struct SoilProperties
{
	double Density;
	double Porosity;
	double ResidualWaterSaturation;
	double HeatCapacity;
	double ThermalConductivity;
	double Wpar;
	double Mpar;
	int rSFC;
	bool IsGSLIB;
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
};

struct NonisothermalProperties
{
	bool IsIsothermal;
	double HeatCapacity;
	double ThermalConductivity;
	double TempTransition;
	double TempLiquid;
	double TempSolid;
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
	bool isHeterC;
	bool isHeterK;
	bool isHeterD;
	bool isHeterBC;
	double CorrelationLength;
	int NumberOfCells;
	double GridSize;
	int NumberOfRealizations;
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