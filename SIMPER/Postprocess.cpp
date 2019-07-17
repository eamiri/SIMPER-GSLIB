#include "Postprocess.h"

Postprocess::Postprocess(Mesh m, Properties p, FILE *plot, FILE *node, FILE *area, VectorXd nodalGSLIBCoeffs)
{
	MESH = m;
	PROPS = p;
	OutputFile = plot;
	AreaAnalysisFile = area;
	NodePlotFile = node;
	nGP = PROPS.GQ.NumberOfPoints;

	GP.GP(nGP);
	noel = MESH.NumberOfElements;
	nond = MESH.NumberOfNodes;
	ndoe = MESH.ElementNumberOfNodes;
	NodalGSLIBCoeffs = nodalGSLIBCoeffs;
}

void Postprocess::AreaAnalysis(VectorXd Temp, double solutionTime)
{
	rSFC = PROPS.Soil.rSFC;
	Tliq = PROPS.Nonisothermal.TempLiquid;
	Sres = PROPS.Soil.ResidualWaterSaturation;

	for (int e = 0; e < MESH.NumberOfElements; e++)
	{
		Tsol = MESH.Elements[e].SoilFreezingPoint;
		//
		xNodes = MESH.GetNodesXCoordinates(e, MESH.ElementNumberOfNodes);
		yNodes = MESH.GetNodesYCoordinates(e, MESH.ElementNumberOfNodes);
		AreaElement = MESH.GetElementArea(e, ndoe);
		//
		elementDofs = MESH.GetElementDofs(e, ndoe);
		//
		TempNode = MESH.GetNodalValues(Temp, elementDofs);
		//
		iFrozen = 0; iThawed = 0; iSlushy = 0;
		for (int n = 0; n < ndoe; n++)
		{
			if (TempNode(n) <= PROPS.Nonisothermal.TempSolid)
			{
				iFrozen++;
			}
			else if (TempNode(n) > PROPS.Nonisothermal.TempSolid && TempNode(n) < PROPS.Nonisothermal.TempLiquid)
			{
				iSlushy++;
			}
			else
			{
				iThawed++;
			}
		}

		frozenArea += iFrozen * AreaElement / ndoe;
		slushyArea += iSlushy * AreaElement / ndoe;
		thawedArea += iThawed * AreaElement / ndoe;
	}

	fprintf(AreaAnalysisFile, "%e, %e, %e, %e\n", solutionTime, frozenArea, thawedArea, slushyArea);
	fflush(AreaAnalysisFile);
}

void Postprocess::Plot(VectorXd Temp, double solutionTime)
{
	//
	rSFC = PROPS.Soil.rSFC;
	Tliq = PROPS.Nonisothermal.TempLiquid;
	Sres = PROPS.Soil.ResidualWaterSaturation;
	Wpar = PROPS.Soil.Wpar;
	Mpar = PROPS.Soil.Mpar;
	//
	derTemp = derTemp.Zero(nond, 2);
	waterSat = waterSat.Zero(nond);
	iceSat = iceSat.Zero(nond);
	distance = distance.Zero(nond);

	for (int e = 0; e < MESH.NumberOfElements; e++)
	{
		Tsol = MESH.Elements[e].SoilFreezingPoint;
		//
		xNodes = MESH.GetNodesXCoordinates(e, MESH.ElementNumberOfNodes);
		yNodes = MESH.GetNodesYCoordinates(e, MESH.ElementNumberOfNodes);
		AreaElement = MESH.GetElementArea(e, ndoe);
		//
		elementDofs = MESH.GetElementDofs(e, ndoe);
		//
		TempNode = MESH.GetNodalValues(Temp, elementDofs);
		//
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

					SaturationFunctions SATFUNCS(Temp(elementDofs(n)), Tsol, Tliq, Sres, PROPS.Soil.IsSaturated);
					waterSat(elementDofs(n)) += dNodeGP * SATFUNCS.Swat;
					iceSat(elementDofs(n)) += dNodeGP * SATFUNCS.Sice;

					derTemp(elementDofs(n), 0) +=  dNodeGP * GradTemp(0);
					derTemp(elementDofs(n), 1) += dNodeGP * GradTemp(1);

					distance(elementDofs(n)) += dNodeGP;
				}
			}
		}
	}

	for (int n = 0; n < nond; n++)
	{
		derTemp(n) /= distance(n);
		derTemp(n) /= distance(n);
		waterSat(n) /= distance(n);
		iceSat(n) /= distance(n);
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
			fprintf(OutputFile, "%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n", MESH.Nodes[n].Coordinates.x, MESH.Nodes[n].Coordinates.y, Temp(n), derTemp(n, 0), derTemp(n, 1), waterSat(n), iceSat(n), NodalGSLIBCoeffs(n));
		}
	}
	else
	{
		fprintf(OutputFile, "variables =\"X\" \"Y\" \"T\" \"Tx\" \"Ty\" \"Sw\" \"Si\" \n");
		fprintf(OutputFile, "ZONE N = %5.0d, E = %5.0d, ZONETYPE = FEQuadrilateral, DATAPACKING = POINT\n", nond, noel);
		for (int n = 0; n < nond; n++)
		{
			fprintf(OutputFile, "%e\t%e\t%e\t%e\t%e\t%e\t%e\n", MESH.Nodes[n].Coordinates.x, MESH.Nodes[n].Coordinates.y, Temp(n), derTemp(n, 0), derTemp(n, 1), waterSat(n), iceSat(n));
		}
	}

	for (int e = 0; e < noel; e++)
	{
		fprintf(OutputFile, "%i\t%i\t%i\t%i\t\n", MESH.Elements[e].Nodes[0].n, MESH.Elements[e].Nodes[1].n, MESH.Elements[e].Nodes[2].n, MESH.Elements[e].Nodes[3].n);
	}

	fflush(OutputFile);
}