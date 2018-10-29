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
	bool isGSLIB;
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
	vector<int> Node;
	vector<double> Value;
};

struct GaussQuadrature
{
	int NumberOfPoints;
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