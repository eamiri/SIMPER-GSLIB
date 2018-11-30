#include "SimperInclude.h"

void POSTPROCESS(VectorXd Temp, double solutionTime)
{
	VectorXd xNodes, yNodes;
	VectorXd TempNode;
	VectorXd GradTemp;
	VectorXi elementDofs;
	MatrixXd derTemp;
	FEMFunctions femFunctions;
	MatrixXd Bmat;
	RowVectorXd SF;
	
	double waterSat, iceSat;
	double GPi, GPj, Wi, Wj, globalGPx, globalGPy, dNodeGP, xNode, yNode;
	double distance;

	int rSFC;
	double Tsol, Tliq, Sres, Wpar, Mpar;

	//
	rSFC = PROPS.Soil.rSFC;
	Tsol = PROPS.Nonisothermal.TempSolid;
	Tliq = PROPS.Nonisothermal.TempLiquid;
	Sres = PROPS.Soil.ResidualWaterSaturation;
	Wpar = PROPS.Soil.Wpar;
	Mpar = PROPS.Soil.Mpar;
	//

	derTemp = derTemp.Zero(nond, 2);

	for (int e = 0; e < MESH.NumberOfElements; e++)
	{
		xNodes = MESH.GetNodesXCoordinates(e, MESH.ElementNumberOfNodes);
		yNodes = MESH.GetNodesYCoordinates(e, MESH.ElementNumberOfNodes);
		//
		elementDofs = MESH.GetElementDofs(e, ndoe);
		//
		TempNode = MESH.GetNodalValues(Temp, elementDofs);
		//
		distance = 0.0;

		for (int i = 0; i < nGP; i++)
		{
			for (int j = 0; j < nGP; j++)
			{
				GPi = GP.Points[i];
				GPj = GP.Points[j];
				Wi = GP.Weights[i];
				Wj = GP.Weights[j];
				//
				femFunctions.Calculate(xNodes, yNodes, GPi, GPj);
				SF = femFunctions.SF;
				Bmat = femFunctions.Bmat;
				GradTemp = Bmat * TempNode;
				//
				globalGPx = SF * xNodes;
				globalGPy = SF * yNodes;
				//
				
				for (int n = 0; n < ndoe; n++)
				{
					xNode = MESH.Elements[e].Nodes[n].Coordinates.x;
					yNode = MESH.Elements[e].Nodes[n].Coordinates.y;
					dNodeGP = sqrt((globalGPx - xNode)*(globalGPx - xNode) + (globalGPy - yNode)*(globalGPy - yNode));

					derTemp(elementDofs(n), 0) +=  dNodeGP * GradTemp(0);
					derTemp(elementDofs(n), 1) += dNodeGP * GradTemp(1);

					distance += dNodeGP;
				}
			}
		}

		for (int n = 0; n < ndoe; n++)
		{
			derTemp(elementDofs(n), 0) /= distance;
			derTemp(elementDofs(n), 1) /= distance;
		}
	}

	if (PROPS.PlotNodes.size())
	{
		fprintf(NodePlotFile, "\n%e\t", solutionTime);
		for (int n = 0; n < PROPS.PlotNodes.size(); n++)
		{
			fprintf(NodePlotFile, "%e\t", Temp(PROPS.PlotNodes[n] - 1));
		}

		fflush(NodePlotFile);
	}	

	if (PROPS.Soil.IsGSLIB)
	{
		fprintf(OutputFile, "variables =\"X\" \"Y\" \"T\" \"Tx\" \"Ty\" \"Sw\" \"Si\" \"GSLIB Coeff\"\n");
		fprintf(OutputFile, "ZONE N = %5.0d, E = %5.0d, ZONETYPE = FEQuadrilateral, DATAPACKING = POINT\n", nond, noel);
		for (int n = 0; n < nond; n++)
		{
			SaturationFunctions SATFUNCS(Temp(n), Tsol, Tliq, Sres, PROPS.Soil.IsSaturated);
			waterSat = SATFUNCS.Swat;
			iceSat = SATFUNCS.Sice;
			fprintf(OutputFile, "%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n", MESH.Nodes[n].Coordinates.x, MESH.Nodes[n].Coordinates.y, Temp(n), derTemp(n, 0), derTemp(n, 1), waterSat, iceSat, NodalGSLIBCoeffs(n));
		}
	}
	else
	{
		fprintf(OutputFile, "variables =\"X\" \"Y\" \"T\" \"Tx\" \"Ty\" \"Sw\" \"Si\" \n");
		fprintf(OutputFile, "ZONE N = %5.0d, E = %5.0d, ZONETYPE = FEQuadrilateral, DATAPACKING = POINT\n", nond, noel);
		for (int n = 0; n < nond; n++)
		{
			SaturationFunctions SATFUNCS(Temp(n), Tsol, Tliq, Sres, PROPS.Soil.IsSaturated);
			waterSat = SATFUNCS.Swat;
			iceSat = SATFUNCS.Sice;
			fprintf(OutputFile, "%e\t%e\t%e\t%e\t%e\t%e\t%e\n", MESH.Nodes[n].Coordinates.x, MESH.Nodes[n].Coordinates.y, Temp(n), derTemp(n, 0), derTemp(n, 1), waterSat, iceSat);
		}
	}
	

	for (int e = 0; e < noel; e++)
	{
		fprintf(OutputFile, "%i\t%i\t%i\t%i\t\n", MESH.Elements[e].Nodes[0].n, MESH.Elements[e].Nodes[1].n, MESH.Elements[e].Nodes[2].n, MESH.Elements[e].Nodes[3].n);
	}

	fflush(OutputFile);
}