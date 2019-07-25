#include "SimperInclude.h"

static string SimPerBuildDate(__DATE__);

const string PropertiesInputFile = "../Inputs/input_properties.dat";
const string MeshInputFile = "../Inputs/input_mesh.msh";
GlobalOptions Options;
FILE *OutputFile, *NodePlotFile, *AreaAnalysisFile, *TalikAreaFile, *PermafrostAreaFile;
ofstream outputFile;

bool Inputs(string propsInputFile, string meshInputFile);
void Simulate(int iReal);

VectorXd Temp_0;
VectorXd TempDot_0;
VectorXd Residual;
VectorXd Temp, TempMin, TempMax;
VectorXd TempDot;
vector<int> DirichletBoundary;
vector<int> TopBoundary;
vector<int> BottomBoundary;

MatrixXd BCInputData;

int main(int argc, char* argv[])
{
	//double      t;
	string      filebase;
	clock_t     t0, t1;          //computational time markers

	Options.Pause = true;
	Options.Version = "1.9.0";

	if (!Options.Silent) {
		int year = StringToInteger(SimPerBuildDate.substr(SimPerBuildDate.length() - 4, 4).c_str());
		cout << "============================================================= "  << endl;
		cout << "               SIMPER: Simulate Permafrost					   "  << endl;
		cout << "              a robust permafrost simulator			       "  << endl;
		cout << "    Copyright 2015-" << year << ", the SIMPER Development Team"  << endl;
		cout << "                    Version " << Options.Version                 << endl;
		cout << "                BuildDate " << SimPerBuildDate				      << endl;
		cout << "============================================================= "  << endl;
	}

	//INPUT AND OUTPUT
	if (!Inputs(PropertiesInputFile, MeshInputFile))
	{
		cout << "Input File Error!" << endl;
	}
	else
	{
		cout << "Number of Nodes =        " << MESH.NumberOfNodes << endl;
		cout << "Number of Elements =     " << MESH.NumberOfElements << endl;
		cout << "Number of Gauss Points = " << PROPS.GQ.NumberOfPoints << endl;
	}

	//INITIAL AND BOUNDARY CONDITIONS
	IBConditions Conditions(&MESH, &PROPS);
	Temp = Conditions.Temp;
	TempDot = Conditions.TempDot;
	Temp_0 = Conditions.Temp_0;
	TempDot_0 = Conditions.TempDot_0;
	Residual = Conditions.Residual;
	DirichletBoundary = Conditions.DirichletBC;
	TopBoundary = Conditions.TopBC;
	BottomBoundary = Conditions.BottomBC;
	PROPS.PlotNodes = Conditions.PlotNodes;
	BCInputData.resize(PROPS.Solution.MaxTimestep, 2);
	BCInputData = Conditions.BCInputData;

	//PLOT FILES
	NodePlotFile = fopen("../Results/NodePlot.plt", "w");
	fprintf(NodePlotFile, "TITLE = \"Temperature profile of the nodes\"\n");
	fprintf(NodePlotFile, "VARIABLES = \"<i>t</i>(sec)\"");
	for (int n = 0; n < PROPS.PlotNodes.size(); n++)
	{
		fprintf(NodePlotFile, ", \"Node %i(<sup>o</sup>C)\"", n + 1);
	}
	
	OutputFile = fopen("../Results/OutPut.plt", "w"); 
	AreaAnalysisFile = fopen("../Results/AreaAnalysisTimeSeries.csv", "w");
	TalikAreaFile = fopen("../Results/TalikArea.csv", "w");
	PermafrostAreaFile = fopen("../Results/PermafrostArea.csv", "w");
	fprintf(AreaAnalysisFile, "SolutionTime,FrozenArea,ThawedArea,SlushyArea\n");
	fprintf(TalikAreaFile, "Year,TalikArea\n");
	fprintf(PermafrostAreaFile, "Year,PermafrostArea\n");
	//SOLUTION
	for (int iReal = 0; iReal < PROPS.GSLIB.NumberOfRealizations; iReal++)
	{
		Simulate(iReal);
	}
	

	fclose(OutputFile);
	fclose(NodePlotFile);
	fclose(AreaAnalysisFile);
	fclose(TalikAreaFile);
	fclose(PermafrostAreaFile);

	cin.get();
	
	return 0;
}

