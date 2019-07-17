#include "SimperInclude.h"

using namespace std;
using namespace Eigen;

SimplicialLLT<SparseMatrix<double>> Solver;

int iTimestep, iIteration, iPlot;
int *maxTimestep = &PROPS.Solution.MaxTimestep;
int *maxIterations = &PROPS.Solution.MaxIterations;
int *plotInterval = &PROPS.Solution.PlotInterval;
int talikYear, permafrostYear;

double *TolResidual = &PROPS.Solution.TolPsi;

void InitializeSolution();
void ComputePotentialStar();
Postprocess POSTPROCESS(MESH, PROPS, OutputFile, NodePlotFile, AreaAnalysisFile, nGP, NodalGSLIBCoeffs);
void UpdateTopBC(int iTimestep);
void Simulate();

VectorXd delTempStar, TempStar, TempHat, TempDotStar;
VectorXd dResidual;
VectorXd xNodes, yNodes;
VectorXd TempNode, TempNodeDot, TempNodeHat, TempNode_0, GSLIBCoeffsNode;
VectorXd GradTemp;
VectorXi elementDofs;
VectorXd bcResidual;

SparseMatrix<double> StiffSparse;

MatrixXd dSTIFF;
MatrixXd dMmat, Mmat, Fmat;
MatrixXd Jmat, dJmat;
MatrixXd xyNodes;
FEMFunctions femFunctions;

RowVectorXd SF;
MatrixXd dSF;
MatrixXd Jacob;
MatrixXd Bmat;

TrustRegion TR;

double *gammaNewmark = &PROPS.Solution.NewmarkGamma;
double *deltaTime = &PROPS.Solution.DeltaTime;
double solutionTime;
double trustRegionRadius, maxTrustRegionRadius;
double Potential, PotentialStar, Potential_0, PotetialP, PI1, PI2, PI3, PI4, NormResidual, NormResidualP, NormResidual_0;
double GPi, GPj, Wi, Wj;
double Sair, ISair, ISairxT, Swat, ISwat, dSwat, ISwatxT, IdSwatxT, Sice, ISice, IdSice, ISicexT, dSice, IdSicexT;
double CPar, ICparxT, ICpar, Kpar, Cpar;
double Hfun, Gfun;
double TempG, TempGHat, TempGDot, GradTempGradTemp, GSLIBCoeffG;
double detJacob;
double ErrorResidual, ErrorTemp, ErrorPotential, Reaction;
double m_Temp, m_TempStar, TR_Ratio, actualReduce, modelReduce;
double FenFlux;

bool IsSaturated;

int rSFC;
double Cair, Kair, Dair, Cwat, Kwat, Dwat, Cice, Kice, Dice, Lhea, npor, HydCon, Csol, Ksol, Dsol, Sres, Tsol, Tliq, Wpar, Mpar;

void InitializeSolution()
{
	StiffSparse.resize(nond, nond);
	Residual = Residual.Zero(nond);
	bcResidual.resize(DirichletBoundary.size());
	//delTempStar = delTempStar.Zero(nond);
	TempStar = Temp;
	TempHat = Temp_0 + (1.0 - *gammaNewmark) * (*deltaTime) * TempDot_0;
	iIteration = 0;
	TempMin = Temp;
	//TempMax = Temp;
}
 
void ComputePotentialStar()
{
	PotentialStar = 0.0;
	TempDotStar = (TempStar - TempHat) / (*gammaNewmark * *deltaTime);
	for (int e = 0; e < noel; e++)
	{
		Tsol = MESH.Elements[e].SoilFreezingPoint;
		//
		dMmat = dMmat.Zero(ndoe, ndoe);
		Mmat = Mmat.Zero(ndoe, 1);
		//
		PI1 = 0.0;
		PI2 = 0.0;
		PI3 = 0.0;
		PI4 = 0.0;
		//
		Jmat = Jmat.Zero(ndoe, ndoe);
		dJmat = dJmat.Zero(ndoe, ndoe);
		//
		xNodes = MESH.GetNodesXCoordinates(e, ndoe);
		yNodes = MESH.GetNodesYCoordinates(e, ndoe);
		//
		elementDofs = MESH.GetElementDofs(e, ndoe);
		//
		TempNode = MESH.GetNodalValues(TempStar, elementDofs);
		TempNodeDot = MESH.GetNodalValues(TempDotStar, elementDofs);
		TempNodeHat = MESH.GetNodalValues(TempHat, elementDofs);
		TempNode_0 = MESH.GetNodalValues(Temp_0, elementDofs);
		// GSLIB
		Csol = MESH.Elements[e].SoilHeatCapacity;
		Dsol = MESH.Elements[e].SoilDensity;
		Ksol = MESH.Elements[e].SoilThermalConductivity;

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
				dSF = femFunctions.dSF;
				Jacob = femFunctions.J;
				Bmat = femFunctions.Bmat;
				//
				TempG = SF * TempNode;
				TempGHat = SF * TempNodeHat;
				TempGDot = SF * TempNodeDot;
				GradTemp = Bmat * TempNode;
				//
				SaturationFunctions SATFUNCS(TempG, Tsol, Tliq, Sres, IsSaturated);
				Swat = SATFUNCS.Swat;
				ISwat = SATFUNCS.ISwat;
				ISwatxT = SATFUNCS.ISwatxT;
				dSwat = SATFUNCS.dSwat;
				Sice = SATFUNCS.Sice;
				ISice = SATFUNCS.ISice;
				IdSice = SATFUNCS.IdSice;
				ISicexT = SATFUNCS.ISicexT;
				IdSicexT = SATFUNCS.IdSicexT;
				dSice = SATFUNCS.dSice;
				dSice = SATFUNCS.dSice;
				Sair = SATFUNCS.Sair;
				ISairxT = SATFUNCS.ISairxT;
				ISair = SATFUNCS.ISair;
				//
				Cpar = npor * (Swat * Dwat * Cwat + Sice * Dice * Cice + Sair * Dair * Cair) + (1.0 - npor)* Dsol* Csol + npor * Dice * Lhea * (-dSice);
				ICparxT = npor * (ISwatxT * Dwat * Cwat + ISicexT * Dice * Cice + ISairxT * Dair * Cair) + (1 - npor) * Dsol * Csol * (TempG * TempG / 2.0 - Tsol * Tsol / 2.0) + npor * Dice * Lhea * (-IdSicexT);
				//ICpar = npor * (ISwat * Dwat * Cwat + ISice * Dice * Cice) + (1.0 - npor) * Dsol * Csol * TempG + npor * Dice * Lhea*(Swat - Sres);
				ICpar = npor * (ISwat * Dwat * Cwat + ISice * Dice * Cice + ISair * Dair * Cair) + (1.0 - npor) * Dsol * Csol * (TempG - Tsol) + npor * Dice * Lhea*(-IdSice);
				Kpar = pow(Kwat, (Swat * npor)) * pow(Kice, (Sice * npor))* pow(Ksol, (1.0 - npor));
				//
				detJacob = Jacob.determinant();
				//
				Hfun = ICpar;
				Gfun = ICparxT;
				PI1 += 1.0 / (*deltaTime * *gammaNewmark) * Wi * Wj * (Gfun) * detJacob;
				PI2 += TempGHat / (*deltaTime * *gammaNewmark) * Wi * Wj * Hfun * detJacob;
				GradTempGradTemp = (GradTemp.transpose()) * GradTemp;
				PI3 += 0.5 * Wi * Wj * Kpar * GradTempGradTemp * detJacob;
				//PI4 += Swat * Wi * Wj * FenFlux * TempG * detJacob;
			}
		}

		PotentialStar += PI1 - PI2 + PI3 + PI4;
	}
}

void UpdateTopBC(int iTimestep)
{
	for (int i = 0; i < TopBoundary.size(); i++)
	{
		if (PROPS.GSLIB.isHeterBC && PROPS.Soil.IsGSLIB)
		{
			Temp(TopBoundary[i]) = NodalGSLIBCoeffs(TopBoundary[i]) * BCInputData(iTimestep, 1);
		}
		else
		{
			Temp(TopBoundary[i]) = BCInputData(iTimestep, 1);
		}
	}
}

void Simulate()
{
	IsSaturated = PROPS.Soil.IsSaturated;

	double trRatioParameter = abs(PROPS.Nonisothermal.TempLiquid - PROPS.Nonisothermal.TempSolid);
	maxTrustRegionRadius = 5.0E+3 * trRatioParameter;
	trustRegionRadius = 1000 * trRatioParameter;
	iPlot = 0;

	//
	Cair = PROPS.Fluid.AHeatCapacity;
	Kair = PROPS.Fluid.AThermalConductivity;
	Dair = PROPS.Fluid.ADensity;
	Cwat = PROPS.Fluid.LHeatCapacity;
	Kwat = PROPS.Fluid.LThermalConductivity;
	Dwat = PROPS.Fluid.LDensity;
	Cice = PROPS.Fluid.SHeatCapacity;
	Kice = PROPS.Fluid.SThermalConductivity;
	Dice = PROPS.Fluid.SDensity;
	Lhea = PROPS.Fluid.LatentHeat;
	npor = PROPS.Soil.Porosity;
	Csol = PROPS.Soil.HeatCapacity;
	Ksol = PROPS.Soil.ThermalConductivity;
	Dsol = PROPS.Soil.Density;
	//
	rSFC = PROPS.Soil.rSFC;
	Tliq = PROPS.Nonisothermal.TempLiquid;
	Sres = PROPS.Soil.ResidualWaterSaturation;
	Wpar = PROPS.Soil.Wpar;
	Mpar = PROPS.Soil.Mpar;
	//

	for (iTimestep = 0; iTimestep < *maxTimestep; iTimestep++)
	{
		//Second to Year, Month, Day Conversion
		solutionTime = (iTimestep + 1) * *deltaTime;
		double Year = floor(solutionTime / (3600.0 * 24.0 * 365.0));
		double Month = floor(solutionTime / (3600.0 * 24.0 * 31.0));
		double Week = floor(solutionTime / (3600.0 * 24.0 * 7));
		double Day = floor(solutionTime / (3600.0 * 24.0));
		//
		printf("======================================================================================================================================================================");
		// update boundary conditions
		if (PROPS.Solution.IsInputBC)
		{
			UpdateTopBC(iTimestep);
		}
		else if (!PROPS.Solution.IsInputBC && iTimestep == 0)
		{
			for (int i = 0; i < TopBoundary.size(); i++)
			{
				if (PROPS.GSLIB.isHeterBC && PROPS.Soil.IsGSLIB)
				{
					Temp(TopBoundary[i]) = NodalGSLIBCoeffs(TopBoundary[i]) * PROPS.BCs.Value[0];
				}
				else
				{
					Temp(TopBoundary[i]) = PROPS.BCs.Value[0];
				}
			}
		}
		
		// initialize solution vectors and variables
		InitializeSolution();
		bool isDelTempStarApproved = true;
		iIteration = 0;
		while (true)
		{
			iIteration++;
			printf("\nSTEP= %6i \tITRN= %3i", iTimestep + 1, iIteration);
			if (isDelTempStarApproved)
			{
				Residual = Residual.Zero(MESH.NumberOfNodes);
				TempDot = (Temp - TempHat) / ((*gammaNewmark) * (*deltaTime));
				Potential = 0.0;
				StiffSparse.setZero();

				for (int e = 0; e < noel; e++)
				{
					Tsol = MESH.Elements[e].SoilFreezingPoint;
					//
					dMmat = dMmat.Zero(ndoe, ndoe);
					Mmat = Mmat.Zero(ndoe, 1);
					Fmat = Fmat.Zero(ndoe, 1);
					//
					PI1 = 0.0;
					PI2 = 0.0;
					PI3 = 0.0;
					PI4 = 0.0;
					//
					Jmat = Jmat.Zero(ndoe, ndoe);
					dJmat = dJmat.Zero(ndoe, ndoe);
					//
					xNodes = MESH.GetNodesXCoordinates(e, ndoe);
					yNodes = MESH.GetNodesYCoordinates(e, ndoe);
					//
					elementDofs = MESH.GetElementDofs(e, ndoe);
					//
					TempNode = MESH.GetNodalValues(Temp, elementDofs);
					TempNodeDot = MESH.GetNodalValues(TempDot, elementDofs);
					TempNodeHat = MESH.GetNodalValues(TempHat, elementDofs);
					TempNode_0 = MESH.GetNodalValues(Temp_0, elementDofs);
					// GSLIB
					Csol = MESH.Elements[e].SoilHeatCapacity;
					Dsol = MESH.Elements[e].SoilDensity;
					Ksol = MESH.Elements[e].SoilThermalConductivity;
					HydCon = MESH.Elements[e].SoilHydraulicConductivity;

					for (int i = 0; i < nGP; i++)
					{
						for (int j = 0; j < nGP; j++)
						{
							GP.GP(nGP);
							GPi = GP.Points[i];
							GPj = GP.Points[j];
							Wi = GP.Weights[i];
							Wj = GP.Weights[j];
							//
							femFunctions.Calculate(xNodes, yNodes, GPi, GPj);
							SF = femFunctions.SF;
							dSF = femFunctions.dSF;
							Jacob = femFunctions.J;
							Bmat = femFunctions.Bmat;
							//
							TempG = SF * TempNode;
							TempGHat = SF * TempNodeHat;
							TempGDot = SF * TempNodeDot;
							GradTemp = Bmat * TempNode;							
							//
							SaturationFunctions SATFUNCS(TempG, Tsol, Tliq, Sres, IsSaturated);
							Swat = SATFUNCS.Swat;
							ISwat = SATFUNCS.ISwat;
							ISwatxT = SATFUNCS.ISwatxT;
							dSwat = SATFUNCS.dSwat;
							Sice = SATFUNCS.Sice;
							ISice = SATFUNCS.ISice;
							IdSice = SATFUNCS.IdSice;
							ISicexT = SATFUNCS.ISicexT;
							IdSicexT = SATFUNCS.IdSicexT;
							dSice = SATFUNCS.dSice;
							Sair = SATFUNCS.Sair;
							ISairxT = SATFUNCS.ISairxT;
							ISair = SATFUNCS.ISair;
							//
							Cpar = npor * (Swat * Dwat * Cwat + Sice * Dice * Cice + Sair * Dair * Cair) + (1.0 - npor)* Dsol* Csol + npor * Dice * Lhea * (-dSice);
							ICparxT = npor * (ISwatxT * Dwat * Cwat + ISicexT * Dice * Cice + ISairxT * Dair * Cair) + (1 - npor) * Dsol * Csol * (TempG * TempG / 2.0 - Tsol * Tsol / 2.0) + npor * Dice * Lhea * (-IdSicexT);
							//ICpar = npor * (ISwat * Dwat * Cwat + ISice * Dice * Cice) + (1.0 - npor) * Dsol * Csol * TempG + npor * Dice * Lhea*(Swat - Sres);
							ICpar = npor * (ISwat * Dwat * Cwat + ISice * Dice * Cice + ISair * Dair * Cair) + (1.0 - npor) * Dsol * Csol * (TempG - Tsol) + npor * Dice * Lhea*(-IdSice);
							Kpar = pow(Kwat, (Swat * npor)) * pow(Kice, (Sice * npor))* pow(Ksol, (1.0 - npor));
							//
							SaturationFunctions FenSATFUNCS(TempG, -0.1, 0, 0, IsSaturated);
							FenFlux = 0.0;
							if (BCInputData(iTimestep, 1) < 0)
							{
								if (MESH.Elements[e].SoilType == "Fen")
								{
									FenFlux = Dwat * Cwat * (0.2) * (-HydCon * FenSATFUNCS.Swat * 100.0 * 0.01);
								}
								else
								{
									FenFlux = Dwat * Cwat * (0.2) * (-HydCon * FenSATFUNCS.Swat * 0.01);
								}
							}
							
							//
							detJacob = Jacob.determinant();
							Jmat = Jmat + Wi * Wj*Bmat.transpose()*Kpar*Bmat*detJacob;
							dMmat = dMmat + Wi * Wj*SF.transpose()*Cpar*SF*detJacob;
							Mmat = Mmat + Wi * Wj*SF.transpose()*Cpar*TempGDot*detJacob;
							Fmat = Fmat + Wi * Wj*SF.transpose()*FenFlux*detJacob;
							//
							Hfun = ICpar;
							Gfun = ICparxT;
							PI1 += 1.0 / (*deltaTime * *gammaNewmark)*Wi*Wj*(Gfun)*detJacob;
							PI2 += TempGHat / (*deltaTime * *gammaNewmark)*Wi*Wj*Hfun*detJacob;
							GradTempGradTemp = (GradTemp.transpose()) * GradTemp;
							PI3 += 0.5 * Wi * Wj * Kpar * GradTempGradTemp * detJacob;
							//PI4 += Wi * Wj * FenFlux * TempG * detJacob;
						}
					}

					Potential += PI1 - PI2 + PI3 + PI4;
					dSTIFF = Jmat + dMmat / (*gammaNewmark * *deltaTime);
					dResidual = Jmat * TempNode + Mmat + Fmat;

					for (int m = 0; m < ndoe; m++)
					{
						Residual(elementDofs(m)) += dResidual(m);
						for (int n = 0; n < ndoe; n++)
						{
							bool isBoudary = (MESH.Nodes[elementDofs(m)].IsBoundary || MESH.Nodes[elementDofs(n)].IsBoundary);
							if (dSTIFF(m, n) != 0 && !isBoudary)
							{
								StiffSparse.coeffRef(elementDofs(m), elementDofs(n)) += dSTIFF(m, n);
							}
						}
					}
				}

				if (iIteration == 1)
				{
					Potential_0 = Potential;
					PotetialP = 0.0;
					NormResidualP = 0.0;
					NormResidual_0 = Residual.norm();
				}

				for (int iBC = 0; iBC < DirichletBoundary.size(); iBC++)
				{
					bcResidual(iBC) = Residual(DirichletBoundary[iBC]);
					Residual(DirichletBoundary[iBC]) = 0.0;
					StiffSparse.coeffRef(DirichletBoundary[iBC], DirichletBoundary[iBC]) = 1.0;
				}

				Solver.compute(StiffSparse);
				Solver.determinant();
				if (Solver.info() != Success)
				{
					cout << "Sparse Failed!";
				}

				Reaction = bcResidual.norm();
				NormResidual = Residual.norm();
				ErrorResidual = NormResidual / NormResidual_0;
				NormResidualP = NormResidual;

				ErrorPotential = abs((Potential - PotetialP) / Potential_0);
				PotetialP = Potential;

				StiffSparse.makeCompressed();

			}// end of isDelTempStarApproved condition

			printf("\tERR_RESI= %.3e \tERR_POTEN= %.3e", ErrorResidual, ErrorPotential);

			TR.Iterate(StiffSparse, Residual, trustRegionRadius);

			ErrorTemp = TR.deltaTempTest.norm();
			printf("\tERR_TEMP= %.3e", ErrorTemp);			

			if ((ErrorResidual < *TolResidual && ErrorPotential < *TolResidual) || iIteration > *maxIterations)
			{
				TempDot_0 = TempDot;
				Temp_0 = Temp;
				cout << endl;
				break;
			}

			printf("\tTR_SUCCESS= %s", TR.IsTRSuccess ? "TRUE" : "FALSE");
			printf("\n\t\t\t\tTR_ERR=  %.3e", abs(TR.Error));

			if (!TR.IsTRSuccess && TR.Error > 1.0E+3)
			{
				TempDot_0 = TempDot;
				Temp_0 = Temp;
				cout << endl;
				break;
			}

			delTempStar = TR.deltaTempTest;
			TempStar = Temp + delTempStar;

			ComputePotentialStar();

			m_Temp = Potential;
			m_TempStar = Potential + (Residual.transpose() * delTempStar) + 0.5 * (delTempStar.transpose()) * StiffSparse * delTempStar;

			actualReduce = PotentialStar - Potential;
			modelReduce = m_TempStar - m_Temp;
			TR_Ratio = actualReduce / modelReduce;

			if (abs(actualReduce) / abs(Potential) < 1E-8 || modelReduce == 0)
			{
				TR_Ratio = 1.0;
			}

			printf("\tIS_FNR= %s", TR.IsFullNewtonRaphson ? "TRUE" : "FALSE");
			if (TR_Ratio > 0)
			{
				printf("\tTR_RATIO=  %.3e", TR_Ratio);
			}
			else
			{
				printf("\tTR_RATIO= %.3e", TR_Ratio);
			}

			//update solution estimate
			if (TR_Ratio < 0.25 || TR_Ratio >= 2.0)
			{
				trustRegionRadius = 0.5 * trustRegionRadius;
			}
			else if (TR_Ratio > 3 / 4 && !TR.IsFullNewtonRaphson)
			{
				trustRegionRadius = min(2 * trustRegionRadius, maxTrustRegionRadius);
			}

			// update trust region radius
			if (TR_Ratio > 0.25 && TR_Ratio < 2.0)
			{
				Temp = TempStar;
				isDelTempStarApproved = true;
			}
			else
			{
				isDelTempStarApproved = false;
			}			
		}

		// Analyzing permafrost degradation and talik formations based on min and max temperatures (Permafrost: every two years
		// and talik: every year - based on the official definition)

		//for (int n = 0; n < nond; n++)
		//{
		//	if (Temp(n) < TempMin(n))
		//	{
		//		TempMin(n) = Temp(n);
		//	}

		//	if (Temp(n) > TempMax(n))
		//	{
		//		TempMax(n) = Temp(n);
		//	}
		//}

		if (iTimestep == iPlot * *plotInterval)
		{
			iPlot++;
			//POSTPROCESS(Temp, solutionTime);
		}
	}
}