#include "Postprocess.h"

Postprocess::Postprocess(Mesh mesh, Properties props, MatrixXd nodalGSLIBCoeffs, OutputFileStruc files)
{
	MESH = mesh;
	PROPS = props;
	nGP = PROPS.GQ.NumberOfPoints;
	GP.GP(nGP);
	noel = MESH.NumberOfElements;
	nond = MESH.NumberOfNodes;
	ndoe = MESH.ElementNumberOfNodes;
	if (PROPS.GSLIB.isGSLIB)
	{
		NodalGSLIBCoeffs = nodalGSLIBCoeffs;
	}

	Files = files;
}

void Postprocess::GetTalikArea(VectorXd minTemp, VectorXd maxTemp, int year, int iRealization)
{
	SwLiq = PROPS.Nonisothermal.LiquidSatIndex;
	talikArea = 0.0;
	for (int e = 0; e < MESH.NumberOfElements; e++)
	{	
		int iSoilType = MESH.Elements[e].iSoilType;
		//
		rSFC = PROPS.Soil[iSoilType].rSFC;
		Tliq = PROPS.Nonisothermal.TempLiquid;
		Sres = PROPS.Soil[iSoilType].ResidualWaterSaturation;
		Wpar = PROPS.Soil[iSoilType].Wpar;
		Mpar = PROPS.Soil[iSoilType].Mpar;
		Tsol = MESH.Elements[e].SoilFreezingPoint;
		//
		xNodes = MESH.GetNodesXCoordinates(e, MESH.ElementNumberOfNodes);
		yNodes = MESH.GetNodesYCoordinates(e, MESH.ElementNumberOfNodes);
		AreaElement = MESH.GetElementArea(e, ndoe);
		//
		elementDofs = MESH.GetElementDofs(e, ndoe);
		//
		minTempNode = MESH.GetNodalValues(minTemp, elementDofs);
		maxTempNode = MESH.GetNodalValues(maxTemp, elementDofs);
		//
		iThawed = 0;
		for (int n = 0; n < ndoe; n++)
		{
			SaturationFunctions SATFUNCS(minTempNode(n), Tsol, Tliq, Sres, PROPS.Soil[iSoilType].IsSaturated);
			if (SATFUNCS.Swat >= SwLiq)
			{
				iThawed++;
			}
		}

		talikArea += iThawed * AreaElement / ndoe;
	}

	fprintf(Files.AnnualMinTemperatures, "%i,%i,", iRealization, year);
	fprintf(Files.AnnualMaxTemperatures, "%i,%i,", iRealization, year);
	for (int n = 0; n < nond; n++)
	{
		fprintf(Files.AnnualMinTemperatures, "%e,", minTemp(n));
		fprintf(Files.AnnualMaxTemperatures, "%e,", maxTemp(n));
	}

	fprintf(Files.AnnualMinTemperatures, "\n");
	fprintf(Files.AnnualMaxTemperatures, "\n");

	fprintf(Files.TalikAreaFile, "%i,%i,%e\n", iRealization, year, talikArea);
	fflush(Files.TalikAreaFile);
	fflush(Files.AnnualMinTemperatures);
	fflush(Files.AnnualMaxTemperatures);
}

void Postprocess::GetPermafrostArea(VectorXd minTemp, VectorXd maxTemp, int year, int iRealization)
{
	
	SwSol = PROPS.Nonisothermal.SolidSatIndex;
	if (SwSol < Sres)
	{
		SwSol = Sres;
	}

	permafrostArea = 0.0;
	for (int e = 0; e < MESH.NumberOfElements; e++)
	{
		int iSoilType = MESH.Elements[e].iSoilType;
		//
		rSFC = PROPS.Soil[iSoilType].rSFC;
		Tliq = PROPS.Nonisothermal.TempLiquid;
		Sres = PROPS.Soil[iSoilType].ResidualWaterSaturation;
		Wpar = PROPS.Soil[iSoilType].Wpar;
		Mpar = PROPS.Soil[iSoilType].Mpar;
		Tsol = MESH.Elements[e].SoilFreezingPoint;
		//
		xNodes = MESH.GetNodesXCoordinates(e, MESH.ElementNumberOfNodes);
		yNodes = MESH.GetNodesYCoordinates(e, MESH.ElementNumberOfNodes);
		AreaElement = MESH.GetElementArea(e, ndoe);
		//
		elementDofs = MESH.GetElementDofs(e, ndoe);
		//
		minTempNode = MESH.GetNodalValues(minTemp, elementDofs);
		maxTempNode = MESH.GetNodalValues(maxTemp, elementDofs);
		//
		iFrozen = 0;
		for (int n = 0; n < ndoe; n++)
		{
			SaturationFunctions SATFUNCS(maxTempNode(n), Tsol, Tliq, Sres, PROPS.Soil[iSoilType].IsSaturated);
			if (SATFUNCS.Swat <= SwSol)
			{
				iFrozen++;
			}
		}

		permafrostArea += iFrozen * AreaElement / ndoe;
	}
	fprintf(Files.BiannualMinTemperatures, "%i,%i,", iRealization, year);
	fprintf(Files.BiannualMaxTemperatures, "%i,%i,", iRealization, year);
	for (int n = 0; n < nond; n++)
	{
		fprintf(Files.BiannualMinTemperatures, "%e,", minTemp(n));
		fprintf(Files.BiannualMaxTemperatures, "%e,", maxTemp(n));
	}

	fprintf(Files.BiannualMinTemperatures, "\n");
	fprintf(Files.BiannualMaxTemperatures, "\n");

	fprintf(Files.PermafrostAreaFile, "%i,%i,%e\n", iRealization, year, permafrostArea);
	fflush(Files.PermafrostAreaFile);
	fflush(Files.BiannualMinTemperatures);
	fflush(Files.BiannualMaxTemperatures);
}

void Postprocess::AreaAnalysis(VectorXd Temp, double solutionTime, int iRealization)
{
	SwSol = PROPS.Nonisothermal.SolidSatIndex;
	if (SwSol < Sres)
	{
		SwSol = Sres;
	}

	SwLiq = PROPS.Nonisothermal.LiquidSatIndex;
	//
	frozenArea = 0.0;
	thawedArea = 0.0;
	slushyArea = 0.0;
	for (int e = 0; e < MESH.NumberOfElements; e++)
	{	
		int iSoilType = MESH.Elements[e].iSoilType;
		//
		rSFC = PROPS.Soil[iSoilType].rSFC;
		Tliq = PROPS.Nonisothermal.TempLiquid;
		Sres = PROPS.Soil[iSoilType].ResidualWaterSaturation;
		Wpar = PROPS.Soil[iSoilType].Wpar;
		Mpar = PROPS.Soil[iSoilType].Mpar;
		//
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
			SaturationFunctions SATFUNCS(TempNode(n), Tsol, Tliq, Sres, PROPS.Soil[iSoilType].IsSaturated);
			if (SATFUNCS.Swat <= SwSol)
			{
				iFrozen++;
			}
			else if (SATFUNCS.Swat > SwSol && SATFUNCS.Swat < SwLiq)
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

	fprintf(Files.AreaAnalysisFile, "%i,%e,%e,%e,%e\n", iRealization, solutionTime, frozenArea, thawedArea, slushyArea);
	fflush(Files.AreaAnalysisFile);
}

void Postprocess::PlotNodes(VectorXd Temp, double solutionTime, int iRealization)
{
	fprintf(Files.NodePlotFile, "\n%e\t", solutionTime);
	for (int n = 0; n < PROPS.PlotNodes.size(); n++)
	{
		fprintf(Files.NodePlotFile, "%e\t", Temp(PROPS.PlotNodes[n] - 1));
	}

	fflush(Files.NodePlotFile);
}

void Postprocess::Plot(VectorXd Temp, double solutionTime, int iRealization)
{
	
	derTemp = derTemp.Zero(nond, 2);
	waterSat = waterSat.Zero(nond);
	iceSat = iceSat.Zero(nond);
	distance = distance.Zero(nond);

	for (int e = 0; e < MESH.NumberOfElements; e++)
	{
		int iSoilType = MESH.Elements[e].iSoilType;
		//
		rSFC = PROPS.Soil[iSoilType].rSFC;
		Tliq = PROPS.Nonisothermal.TempLiquid;
		Sres = PROPS.Soil[iSoilType].ResidualWaterSaturation;
		Wpar = PROPS.Soil[iSoilType].Wpar;
		Mpar = PROPS.Soil[iSoilType].Mpar;
		//
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

					SaturationFunctions SATFUNCS(Temp(elementDofs(n)), Tsol, Tliq, Sres, PROPS.Soil[iSoilType].IsSaturated);
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

	if (PROPS.GSLIB.isGSLIB)
	{
		fprintf(Files.OutputFile, "variables =\"X\" \"Y\" \"T\" \"Tx\" \"Ty\" \"Sw\" \"Si\" \"GSLIB Coeff\"\n");
		fprintf(Files.OutputFile, "ZONE N = %5.0d, E = %5.0d, ZONETYPE = FEQuadrilateral, DATAPACKING = POINT\n", nond, noel);
		for (int n = 0; n < nond; n++)
		{
			fprintf(Files.OutputFile, "%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n", MESH.Nodes[n].Coordinates.x, MESH.Nodes[n].Coordinates.y, Temp(n), derTemp(n, 0), derTemp(n, 1), waterSat(n), iceSat(n), NodalGSLIBCoeffs(n, 0));
		}
	}
	else
	{
		fprintf(Files.OutputFile, "variables =\"X\" \"Y\" \"T\" \"Tx\" \"Ty\" \"Sw\" \"Si\" \n");
		fprintf(Files.OutputFile, "ZONE N = %5.0d, E = %5.0d, ZONETYPE = FEQuadrilateral, DATAPACKING = POINT\n", nond, noel);
		for (int n = 0; n < nond; n++)
		{
			fprintf(Files.OutputFile, "%e\t%e\t%e\t%e\t%e\t%e\t%e\n", MESH.Nodes[n].Coordinates.x, MESH.Nodes[n].Coordinates.y, Temp(n), derTemp(n, 0), derTemp(n, 1), waterSat(n), iceSat(n));
		}
	}

	for (int e = 0; e < noel; e++)
	{
		fprintf(Files.OutputFile, "%i\t%i\t%i\t%i\t\n", MESH.Elements[e].Nodes[0].n, MESH.Elements[e].Nodes[1].n, MESH.Elements[e].Nodes[2].n, MESH.Elements[e].Nodes[3].n);
	}

	fflush(Files.OutputFile);
}