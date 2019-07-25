#ifndef __OUTPUTFILES__
#define __OUTPUTFILES__

#include <stdio.h>
#include <iomanip>
#include <iostream>
using namespace std;
struct OutputFileStruc
{
	FILE* OutputFile;
	FILE* NodePlotFile;
	FILE* AreaAnalysisFile;
	FILE* TalikAreaFile;
	FILE* PermafrostAreaFile;
	FILE* AnnualMinTemperatures;
	FILE* AnnualMaxTemperatures;
	FILE* BiannualMinTemperatures;
	FILE* BiannualMaxTemperatures;
};

#endif 