#include "SimperInclude.h"

static string SimPerBuildDate(__DATE__);

const string PropertiesInputFile = "../Inputs/input_properties.dat";
const string MeshInputFile = "../Inputs/input_mesh.msh";
GlobalOptions Options;

bool Inputs(string propsInputFile, string meshInputFile, int iParallelRealization);
void Simulate(int iReal);
const char* GetOutputFilePath(string outputName, int iRealization);

VectorXd Temp_0;
VectorXd TempDot_0;
VectorXd Residual;
VectorXd Temp, TempMin, TempMax;
VectorXd TempDot;
vector<int> DirichletBoundary;
vector<int> TopBoundary;
vector<int> BottomBoundary;

MatrixXd BCInputData;
int iParallel;

const char * GetOutputFilePath(string outputName, int iRealization)
{
	const char* ch1 = "../Results/R";
	string ch2 = to_string(iRealization);
	string path = ch1;
	path += ch2;
	path += "_";
	path += outputName;
	
	return path.c_str();
}

int main(int argc, char* argv[])
{
	//double      t;
	string      filebase;
	//clock_t     t0, t1;          //computational time markers

	Options.Pause = true;
	Options.Version = "1.5.0";

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

	//PARALLEL REALIZATIONS
	if (argc > 1)
	{

		iParallel = atoi(argv[1]);
	}
	else
	{
		iParallel = 1;
	}

	//INPUT AND OUTPUT
	if (!Inputs(PropertiesInputFile, MeshInputFile, iParallel))
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
	OutputFiles.NodePlotFile = fopen(GetOutputFilePath("NodePlot.plt", iParallel), "w");
	fprintf(OutputFiles.NodePlotFile, "TITLE = \"Temperature profile of the nodes\"\n");
	fprintf(OutputFiles.NodePlotFile, "VARIABLES = \"<i>t</i>(sec)\"");
	for (int n = 0; n < PROPS.PlotNodes.size(); n++)
	{
		fprintf(OutputFiles.NodePlotFile, ", \"Node %i(<sup>o</sup>C)\"", n + 1);
	}

	OutputFiles.OutputFile = fopen(GetOutputFilePath("Output.plt", iParallel), "w");
	OutputFiles.AreaAnalysisFile = fopen(GetOutputFilePath("AreaAnalysisTimeSeries.csv", iParallel), "w");
	OutputFiles.TalikAreaFile = fopen(GetOutputFilePath("TalikArea.csv", iParallel), "w");
	OutputFiles.AnnualMinTemperatures = fopen(GetOutputFilePath("AnnualMinTemperatures.csv", iParallel), "w");
	OutputFiles.AnnualMaxTemperatures = fopen(GetOutputFilePath("AnnualMaxTemperatures.csv", iParallel), "w");
	OutputFiles.PermafrostAreaFile = fopen(GetOutputFilePath("PermafrostAreaFile.csv", iParallel), "w");
	OutputFiles.BiannualMinTemperatures = fopen(GetOutputFilePath("BiannualMinTemperatures.csv", iParallel), "w");
	OutputFiles.BiannualMaxTemperatures = fopen(GetOutputFilePath("BiannualMaxTemperatures.csv", iParallel), "w");

	fprintf(OutputFiles.AreaAnalysisFile, "Realization,SolutionTime,FrozenArea,ThawedArea,SlushyArea\n");
	fprintf(OutputFiles.TalikAreaFile, "Realization,Year,TalikArea\n");
	fprintf(OutputFiles.PermafrostAreaFile, "Realization,Year,PermafrostArea\n");
	fprintf(OutputFiles.AnnualMinTemperatures, "Annual minimum nodal temperatures:\n Realization no., Year, <Nodal Temperatures>\n\n");
	fprintf(OutputFiles.AnnualMaxTemperatures, "Annual maximum nodal temperatures:\n Realization no., Year, <Nodal Temperatures>\n\n");
	fprintf(OutputFiles.BiannualMinTemperatures, "Biannual minimum nodal temperatures:\n Realization no., Year, <Nodal Temperatures>\n\n");
	fprintf(OutputFiles.BiannualMaxTemperatures, "Biannual maximum nodal temperatures:\n Realization no., Year, <Nodal Temperatures>\n\n");

	//SOLUTION
	for (int iReal = 0; iReal < PROPS.GSLIB.NumberOfRealizations; iReal++)
	{
		Simulate(iReal);
	}	

	fclose(OutputFiles.OutputFile);
	fclose(OutputFiles.NodePlotFile);
	fclose(OutputFiles.AreaAnalysisFile);
	fclose(OutputFiles.TalikAreaFile);
	fclose(OutputFiles.PermafrostAreaFile);
	cin.get();
	
	return 0;
}

