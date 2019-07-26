#ifdef _unix_        
#define linux
#elif defined(_WIN32) || defined(WIN32) 
#define _windows_
#endif
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

using namespace std;

struct inputs{
    string account;
    string memory;
    string email;
    string email_case;
    string runtime;
    string job_name;
    string program;
	int nRealization;
};

string GetOutputFilePath(string outputName, int iRealization)
{
	string rlnNo = to_string(iRealization);
	string path = "bashFiles/R";
	path += rlnNo;
	path += "_";
	path += outputName;
	return path;
}

int main(int argc, char* argv[])
{
    inputs input;
    ifstream inputFile;
    inputFile.open("inputBash.dat");
    string line;
    getline(inputFile, line);
    inputFile >> input.account;
    getline(inputFile, line);
    getline(inputFile, line);
	getline(inputFile, line);
    inputFile >> input.memory;
	getline(inputFile, line);
	getline(inputFile, line);
	getline(inputFile, line);
	inputFile >> input.email;
	getline(inputFile, line);
	getline(inputFile, line);
	getline(inputFile, line);
	inputFile >> input.email_case;
	getline(inputFile, line);
	getline(inputFile, line);
	getline(inputFile, line);
	getline(inputFile, line);
	string DD, HH, MM;
	inputFile >> DD >> HH >> MM;
	input.runtime = DD + "-" + HH + ":" + MM;
	getline(inputFile, line);
	getline(inputFile, line);
	getline(inputFile, line);
	inputFile >> input.job_name;
	getline(inputFile, line);
	getline(inputFile, line);
	getline(inputFile, line);
	inputFile >> input.program;
	getline(inputFile, line);
	getline(inputFile, line);
	getline(inputFile, line);
	inputFile >> input.nRealization;

	#ifdef linux
	str = "chmod 755 bashFiles";
	system(str.c_str());
	#endif

	FILE* submitFile;
	submitFile = fopen("bashFiles/submit-to-cluster.sh", "w");
	fprintf(submitFile, "#!/bin/bash\n");
	FILE* bashFile;
	string str;
	int permissionChange;
	for (int iReal = 1; iReal <= input.nRealization; iReal++)
	{
		bashFile = fopen(GetOutputFilePath("BashFile.sh", iReal).c_str(), "w");
		fprintf(bashFile, "#!/bin/bash\n");
		str = "SBATCH --account=" + input.account + "\n";
		fprintf(bashFile, str.c_str());
		str = "SBATCH --mem-per-cpu=" + input.memory + "\n";
		fprintf(bashFile, str.c_str());
		str = "SBATCH --mail-user=" + input.email + "\n";
		fprintf(bashFile, str.c_str());
		str = "SBATCH --mail-type=" + input.email_case + "\n";
		fprintf(bashFile, str.c_str());
		str = "SBATCH --time=" + input.runtime + "\n";
		fprintf(bashFile, str.c_str());
		str = "SBATCH --job-name=" + input.job_name + "_" + to_string(iReal) + "\n\n";
		fprintf(bashFile, str.c_str());
		str = "./" + input.program + " " + to_string(iReal);
		fprintf(bashFile, str.c_str());
		fclose(bashFile);
		str = "chmod 755 " + GetOutputFilePath("BashFile.sh", iReal);

#		ifdef _windows_
		
		#else
		str = "chmod 755 " + GetOutputFilePath("BashFile.sh", iReal);
		permissionChange = system(str.c_str());
		#endif
		
		str = "sbatch R" + to_string(iReal) + "_BashFile.sh\n";
		fprintf(submitFile, str.c_str());
	}

	#ifdef _windows_

	#else
	str = "chmod 755 bashFiles/submit-to-cluster.sh";
	permissionChange = system(str.c_str());
	#endif

    return 0;
}