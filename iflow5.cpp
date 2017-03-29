// iflow.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "param.h"
#include "geom.h"
#include "immersedboundary.h"
#include <iostream>

using namespace std;



bool CylinderMask(int i, int j, int k, int L, int M, int K);
bool InitialConcentrationDistribution(int i, int j, int k, int L, int M, int K);
bool InitialViscosityDistribution(int i, int j, int k, int L, int M, int K);
long double GetBorderConditionU(int lI, int lJ, int lK, int lN);
long double GetBorderConditionV(int lI, int lJ, int lK, int lN);
long double GetBorderConditionW(int lI, int lJ, int lK, int lN);

int _tmain(int argc, _TCHAR* argv[])
{
	// файл параметров должен быть аргументом командной строки
	// но это пока не получилось
	// поэтому

	
	char inpfile[20] = "input.txt";          // файл входных параметров
	cInputParameters inpParams(inpfile);     // заполнение входных параметров из файла
	

	cGeometry geometry(inpParams.domainParams); //основная сетка
	inpParams.timeParams.UpdateSeparationValue(); //  временные параметры

	cImmersedBoundary immersedBoundary(&inpParams); 
	
	 

	//Vectors and operators global pointers defining (Begin)

	//Current time step (Begin)
	//Already defined as class fields
	//Current time step (End)

	//Previous time step (Begin)
	/*
	TRnRelease3dSpace* geometry.U1->fieldPtr = NULL;
	TRnRelease3dSpace* geometry.V1->fieldPtr = NULL;
	TRnRelease3dSpace* geometry.W1->fieldPtr = NULL;
	//Previous time step (End)

	//Intermediate (Begin)
	TRnRelease3dSpace* u_Ptr = NULL;
	TRnRelease3dSpace* v_Ptr = NULL;
	TRnRelease3dSpace* w_Ptr = NULL;
	//Intermediate (End)

	//For control stating (Begin)
	TRnRelease3dSpace* uSPtr = NULL;
	TRnRelease3dSpace* vSPtr = NULL;
	TRnRelease3dSpace* wSPtr = NULL;
	//For control stating (End)

	//Poisson equation right part (Begin)
	TRnRelease3dSpace* rightPartPtr = NULL;
	//Poisson equation right part (End)
	
	//Concentration (Begin)
	TRnRelease3dSpace* C_1Ptr = NULL;
	TRnRelease3dSpace* CPtr = NULL;
	T3dNormalGrid* gridCPtr = NULL;
	long double visValue = 1.5 / fRe;
	//Concentration (End)

	TRnRelease3dSpace* densityPtr = NULL;
	T3dNormalGrid* gridDensityPtr = NULL;
	long double soluteDensityValue = 1.5;
	long double envDensityValue = 1.0;

	//Viscosity (Begin)
	T3dNormalGrid* gridVisPtr = NULL;
	TRnRelease3dSpace* visPtr = NULL;
	//Viscosity (Begin)
	

	T3dVariableDensityOperator* APtr = NULL;
	*/

	/*
	int spaceModelNameRn = 1054;
	T1dPackOperator* ARnPtr = NULL;
	*/

	//Vectors and operators global pointers defining (End)




	//initial conditions (Begin)
	/*
	(*fU) |= 0;
	(*fV) |= 0;
	(*fW) |= 0;
	(*fP) |= 0;
	*/
	//geometry.setFToZero();
	//cField rightPart(geometry.NxP, geometry.NyP, geometry.NzP);// Наверно нужно пересоздавать для U, V ,W !!!!!!!!!!!!!!!
	cMask *maskPtr;
	cField* vpPrevPtr;
	cField* vpCurPtr;

	//initial conditions (End)

	/*
	long double tau = fTimeSeparator->SeparationValue;
	int timeStepsCount = fTimeSeparator->EndIndex;


	bool flagRestruct = true;


	clock_t start = clock();
	clock_t end = clock();
	double diff = 0;

	*/
	
	geometry.setFToZero();
	geometry.SetInitCondVis();
	geometry.SetInitCondDensity();
	//geometry.SetInitCondRho(inpParams);
	geometry.SetInitCondC(inpParams);
	geometry.prepareMesh();
	
	

	
	
	int timeStepNumber;
	long double currentTime,beginTime = inpParams.timeParams.startValue, endTime = inpParams.timeParams.endValue, dTime = inpParams.timeParams.separationValue;

	for (timeStepNumber = 1, currentTime = beginTime; currentTime <= endTime; timeStepNumber++, currentTime += dTime)
	{
		
		//ChangeStifness(timeStepNumber);
		printf("\n\n");
		printf("*****************************************************\n");
		printf("Time step = %d: Starting...\n", timeStepNumber);

		

		


			
			
				
			



		

		
			//SolveTransferEquation_(cnc, cnc_1, uTr, vTr, wTr, rpTr, gridC, 0.0);
			/*
			//Recalculation (Begin)
			for (int i = 0; i<N + 1; i++)
			{
				for (int j = 0; j<L + 1; j++)
				{
					for (int k = 0; k<M + 1; k++)
					{
						vis[i][j][k] = (1 - cnc[i][j][k])*(1 / fRe) + (cnc[i][j][k])*visValue;
						density[i][j][k] = cnc[i][j][k] * (soluteDensityValue - envDensityValue) + envDensityValue;
					}
				}
			}

			//Recalculation (End)

			delete &uTr;
			delete &vTr;
			delete &wTr;
			delete &rpTr;


			end = clock();
			diff = (double)(end - start);
			diff = diff / CLOCKS_PER_SEC;
			printf("Concentration equation: OK. (time = %lf s.)\n", diff);
		}
		//Concentration (End)
		*/

		//Viscosity (Begin)
		/*
		*    if (turbulenceMode==true)
		*    {
		*    printf("\n");
		*    printf("Turbulence model: Starting...\n");
		*    start = clock();
		*
		*    int N = gridVis.GetSeparator1().EndIndex;
		*    int L = gridVis.GetSeparator2().EndIndex;
		*    int M = gridVis.GetSeparator3().EndIndex;
		*    const T3dNumberMask& mask = gridVis.Mask;
		*
		*    const long double* hx = gridVis.GetSeparator1().Separation;
		*    const long double* hy = gridVis.GetSeparator2().Separation;
		*    const long double* hz = gridVis.GetSeparator3().Separation;
		*
		*
		*    //Velocity projecting (Begin)
		*    TRnRelease3dSpace& uTr  = dynamic_cast<TRnRelease3dSpace &>(vis.CreateInstance());
		*    TRnRelease3dSpace& vTr  = dynamic_cast<TRnRelease3dSpace &>(vis.CreateInstance());
		*    TRnRelease3dSpace& wTr  = dynamic_cast<TRnRelease3dSpace &>(vis.CreateInstance());
		*
		*
		*    for (int i=0;i<N+1;i++)
		*    {
		*      for (int j=0;j<L+1;j++)
		*      {
		*        for (int k=0;k<M+1;k++)
		*        {
		*           int mijk = mask[i][j][k];
		*           if (mijk==TFictivePoint)
		*           {
		*                  uTr[i][j][k] = 0; vTr[i][j][k] = 0; wTr[i][j][k] = 0;
		*           }
		*           else if (mijk==TActualPoint)
		*           {
		*                 long double uu = ( u[i-1][j][k] + u[i][j][k])/2;
		*                 uTr[i][j][k] = uu;
		*
		*                 long double vv = ( v[i][j-1][k] + v[i][j][k])/2;
		*                 vTr[i][j][k] = vv;
		*
		*                 long double ww = ( w[i][j][k-1] + w[i][j][k])/2;
		*                 wTr[i][j][k] = ww;
		*           }
		*           else if ( (mijk==TPreNormalBorderPoint) || (mijk==TPreDefinedBorderPoint) || (mijk==TDefinedBorderPoint) || (mijk==TNormalBorderPoint) || (mijk==TEquationBorderPoint) || (mijk==TPreEquationBorderPoint) )
		*           {
		*                throw "Error in TVPProblem::SolveBySplitting: Method is not implemented for this point border type of viscosity.";
		*           }
		*           else
		*           {
		*                throw "Error in TVPProblem::SolveBySplitting: Unknown viscosity point mask type.";
		*           }
		*        } //by k
		*      } //by j
		*    } //by i
		*    //Velocity projecting (End)
		*
		*
		*
		*    //Smogarinsky (Begin)
		*    for (int i=0;i<N+1;i++)
		*    {
		*      for (int j=0;j<L+1;j++)
		*      {
		*        for (int k=0;k<M+1;k++)
		*        {
		*           int mijk = mask[i][j][k];
		*           if (mijk==TFictivePoint)
		*           {
		*                 vis[i][j][k] = 0;
		*           }
		*           else if (mijk==TActualPoint)
		*           {
		*                long double Cs = TProblem::GlobalTurbulenceParameter;
		*
		*                long double delta = powl(hx[i-1]*hy[j-1]*hz[k-1],1.0/3.0);
		*
		*                long double Sxx = (uTr[i+1][j][k] - uTr[i-1][j][k])/(hx[i-1]+hx[i]);
		*                long double Syy = (vTr[i][j+1][k] - vTr[i][j-1][k])/(hy[j-1]+hy[j]);
		*                long double Szz = (wTr[i][j][k+1] - wTr[i][j][k-1])/(hz[k-1]+hz[k]);
		*                long double Sxy = ( (vTr[i+1][j][k] - vTr[i-1][j][k])/(hx[i-1]+hx[i]) + (uTr[i][j+1][k] - uTr[i][j-1][k])/(hy[j-1]+hy[j]) )/2;
		*                long double Sxz = ( (wTr[i+1][j][k] - wTr[i-1][j][k])/(hx[i-1]+hx[i]) + (uTr[i][j][k+1] - uTr[i][j][k-1])/(hz[k-1]+hz[k]) )/2;
		*                long double Syz = ( (vTr[i][j][k+1] - vTr[i][j][k-1])/(hz[k-1]+hz[k]) + (wTr[i][j+1][k] - wTr[i][j-1][k])/(hy[j-1]+hy[j]) )/2;
		*
		*                long double Sc = Sxx*Sxx + Syy*Syy + Szz*Szz + 2*Sxy*Sxy + 2*Sxz*Sxz + 2*Syz*Syz;
		*
		*                long double InvS = powl(2*Sc,1.0/2.0);
		*
		*                long double vSGS = powl(Cs*delta,2)*InvS;
		*
		*                vis[i][j][k] = 1/fRe + vSGS;
		*           }
		*           else if ( (mijk==TPreNormalBorderPoint) || (mijk==TPreDefinedBorderPoint) || (mijk==TDefinedBorderPoint) || (mijk==TNormalBorderPoint) || (mijk==TEquationBorderPoint) || (mijk==TPreEquationBorderPoint) )
		*           {
		*                throw "Error in TVPProblem::SolveBySplitting: Method is not implemented for this point border type of viscosity.";
		*           }
		*           else
		*           {
		*                throw "Error in TVPProblem::SolveBySplitting: Unknown viscosity point mask type.";
		*           }
		*        } //by k
		*      } //by j
		*    } //by i
		*    //Smogarinsky (End)
		*
		*
		*    delete &uTr;
		*    delete &vTr;
		*    delete &wTr;
		*
		*
		*    end = clock();
		*    diff = (double)(end - start);
		*    diff = diff/CLOCKS_PER_SEC;
		*    printf("Turbulence model: OK. (time = %lf s.)\n",diff);
		*    }
		*/
		//Viscosity (End)





		//Saving results (Begin)
		//{
         /*
		if (timeStepNumber % 1 == 0) {
			printf("\n");
			printf("Saving results: Starting...\n");

			start = clock();

			TRnRelease3dSpace& uR = u;
			TRnRelease3dSpace& vR = v;
			TRnRelease3dSpace& wR = w;

			TRnRelease3dSpace& rpTrXR = rpTrX;
			TRnRelease3dSpace& rpTrYR = rpTrY;
			TRnRelease3dSpace& rpTrZR = rpTrZ;


			/*
			TRnRelease3dSpace& uR = u_;
			TRnRelease3dSpace& vR = v_;
			TRnRelease3dSpace& wR = w_;
			*/


			//int len = strlen(TProblem::GlobalCatalog);
			//char *zonesName = new char [len+50];
			//string_copy(zonesName, TProblem::GlobalCatalog, len+50);


			//------------------------------------------------------------------
			// Two format output
			//------------------------------------------------------------------
			//if (timeStepNumber % fBoundary-> == 0 || timeStepNumber == 1)
            /*
			{

				char fileName[40];
				string_print(fileName, 40, "flow%d.dat", 1000 + timeStepNumber);
				//string_concat(zonesName, fileName, len+50);

				//FILE* f = file_open(zonesName, "w");

				FILE* ftp = file_open(fileName, "w");
				//Concentration (Begin) 
				fprintf(ftp, "TITLE = CVP\n");
				fprintf(ftp, "VARIABLES = X,Y,Z,U,V,W,P,Rho,CC\n");
				//Concentration (End) 

				//Viscosity (Begin) 	  	    
				//fprintf(f,"TITLE = VisVP\n");
				//fprintf(f,"VARIABLES = X,Y,Z,Re,U,V,W,P\n");	  
				//Viscosity (End) 	  	    

				T3dNormalGrid& grid = dynamic_cast<T3dNormalGrid&>(*fGrid);
				int N = grid.GetSeparator1().EndIndex;
				int L = grid.GetSeparator2().EndIndex;
				int M = grid.GetSeparator3().EndIndex;
				const long double* x = grid.GetSeparator1().Dimension;
				const long double* y = grid.GetSeparator2().Dimension;
				const long double* z = grid.GetSeparator3().Dimension;
				const T3dNumberMask& mask = grid.Mask;

#if CNPY_VISUALIZATION
				long double *U = new long double[(M + 1)*(L + 1)*(N + 1)];
				long double *V = new long double[(M + 1)*(L + 1)*(N + 1)];
				long double *W = new long double[(M + 1)*(L + 1)*(N + 1)];
				long double *P = new long double[(M + 1)*(L + 1)*(N + 1)];

				long double *FX = new long double[(M + 1)*(L + 1)*(N + 1)];
				long double *FY = new long double[(M + 1)*(L + 1)*(N + 1)];
				long double *FZ = new long double[(M + 1)*(L + 1)*(N + 1)];

				long double *CNC = new long double[(M + 1)*(L + 1)*(N + 1)];
				long double *VIS = new long double[(M + 1)*(L + 1)*(N + 1)];
				long double *DNST = new long double[(M + 1)*(L + 1)*(N + 1)];
#endif

				T3dNormalGrid& gridP = dynamic_cast<T3dNormalGrid&>(*fGridP);
				const T3dNumberMask& maskP = gridP.Mask;

				fprintf(ftp, "ZONE I=%d, J=%d, K=%d, F=POINT\n", N + 1, L + 1, M + 1);
				for (int k = 0; k < M + 1; k++)
				{
					for (int j = 0; j < L + 1; j++)
					{
						for (int i = 0; i < N + 1; i++)
						{
							fprintf(ftp, "%f  %f  %f ", (float)x[i] * (TProblem::GlobalLength), (float)y[j] * (TProblem::GlobalLength), (float)z[k] * (TProblem::GlobalLength));

							long double uu;
							long double vv;
							long double ww;
							long double pp;

							long double fx;
							long double fy;
							long double fz;

							int mijk = mask[i][j][k];
							if (mijk == TFictivePoint) { uu = 0; vv = 0; ww = 0; pp = 0; fx = 0; fy = 0; fz = 0; }
							else
							{
								uu = (uR[i][j][k] + uR[i][j + 1][k] + uR[i][j][k + 1] + uR[i][j + 1][k + 1]) / 4;
								vv = (vR[i][j][k] + vR[i + 1][j][k] + vR[i][j][k + 1] + vR[i + 1][j][k + 1]) / 4;
								ww = (wR[i][j][k] + wR[i + 1][j][k] + wR[i][j + 1][k] + wR[i + 1][j + 1][k]) / 4;

								fx = (rpTrXR[i][j][k] + rpTrXR[i][j + 1][k] + rpTrXR[i][j][k + 1] + rpTrXR[i][j + 1][k + 1]) / 4;
								fy = (rpTrYR[i][j][k] + rpTrYR[i + 1][j][k] + rpTrYR[i][j][k + 1] + rpTrYR[i + 1][j][k + 1]) / 4;
								fz = (rpTrZR[i][j][k] + rpTrZR[i + 1][j][k] + rpTrZR[i][j + 1][k] + rpTrZR[i + 1][j + 1][k]) / 4;

								bool avr1 = ((maskP[i + 1][j][k]) != TFictivePoint);
								bool avr2 = ((maskP[i + 1][j + 1][k]) != TFictivePoint);
								bool avr3 = ((maskP[i + 1][j][k + 1]) != TFictivePoint);
								bool avr4 = ((maskP[i + 1][j + 1][k + 1]) != TFictivePoint);
								bool avr5 = ((maskP[i][j][k]) != TFictivePoint);
								bool avr6 = ((maskP[i][j + 1][k]) != TFictivePoint);
								bool avr7 = ((maskP[i][j][k + 1]) != TFictivePoint);
								bool avr8 = ((maskP[i][j + 1][k + 1]) != TFictivePoint);
								int avr = (int)avr1 + (int)avr2 + (int)avr3 + (int)avr4 + (int)avr5 + (int)avr6 + (int)avr7 + (int)avr8;

								if (avr == 0)
								{
									throw "Error in TVPProblem::SolveBySplitting: Grid is incorrect while saving pressure results.";
								}

								double p1 = p[i + 1][j][k];
								double p2 = p[i + 1][j + 1][k];
								double p3 = p[i + 1][j][k + 1];
								double p4 = p[i + 1][j + 1][k + 1];
								double p5 = p[i][j][k];
								double p6 = p[i][j + 1][k];
								double p7 = p[i][j][k + 1];
								double p8 = p[i][j + 1][k + 1];

								pp = (p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8) / avr;
							}


							//Concentration (Begin)
							const T3dNumberMask& maskC = gridC.Mask;
							long double cc;
							if (mijk == TFictivePoint) { cc = 0; }
							else
							{
								bool avr1 = ((maskC[i + 1][j][k]) != TFictivePoint);
								bool avr2 = ((maskC[i + 1][j + 1][k]) != TFictivePoint);
								bool avr3 = ((maskC[i + 1][j][k + 1]) != TFictivePoint);
								bool avr4 = ((maskC[i + 1][j + 1][k + 1]) != TFictivePoint);
								bool avr5 = ((maskC[i][j][k]) != TFictivePoint);
								bool avr6 = ((maskC[i][j + 1][k]) != TFictivePoint);
								bool avr7 = ((maskC[i][j][k + 1]) != TFictivePoint);
								bool avr8 = ((maskC[i][j + 1][k + 1]) != TFictivePoint);
								int avr = (int)avr1 + (int)avr2 + (int)avr3 + (int)avr4 + (int)avr5 + (int)avr6 + (int)avr7 + (int)avr8;

								if (avr == 0)
								{
									throw "Error in TVPProblem::SolveBySplitting: Grid is incorrect while saving concentration results.";
								}

								double c1 = cnc[i + 1][j][k];
								double c2 = cnc[i + 1][j + 1][k];
								double c3 = cnc[i + 1][j][k + 1];
								double c4 = cnc[i + 1][j + 1][k + 1];
								double c5 = cnc[i][j][k];
								double c6 = cnc[i][j + 1][k];
								double c7 = cnc[i][j][k + 1];
								double c8 = cnc[i][j + 1][k + 1];

								cc = (c1 + c2 + c3 + c4 + c5 + c6 + c7 + c8) / avr;
							}
							//fprintf(f,"%f,",cc);
							//Concentration (End)

							const T3dNumberMask& maskD = gridDensity.Mask;
							long double dnst;
							if (mijk == TFictivePoint) { dnst = 0; }
							else
							{
								bool avr1 = ((maskD[i + 1][j][k]) != TFictivePoint);
								bool avr2 = ((maskD[i + 1][j + 1][k]) != TFictivePoint);
								bool avr3 = ((maskD[i + 1][j][k + 1]) != TFictivePoint);
								bool avr4 = ((maskD[i + 1][j + 1][k + 1]) != TFictivePoint);
								bool avr5 = ((maskD[i][j][k]) != TFictivePoint);
								bool avr6 = ((maskD[i][j + 1][k]) != TFictivePoint);
								bool avr7 = ((maskD[i][j][k + 1]) != TFictivePoint);
								bool avr8 = ((maskD[i][j + 1][k + 1]) != TFictivePoint);
								int avr = (int)avr1 + (int)avr2 + (int)avr3 + (int)avr4 + (int)avr5 + (int)avr6 + (int)avr7 + (int)avr8;

								if (avr == 0)
								{
									throw "Error in TVPProblem::SolveBySplitting: Grid is incorrect while saving concentration results.";
								}

								double d1 = density[i + 1][j][k];
								double d2 = density[i + 1][j + 1][k];
								double d3 = density[i + 1][j][k + 1];
								double d4 = density[i + 1][j + 1][k + 1];
								double d5 = density[i][j][k];
								double d6 = density[i][j + 1][k];
								double d7 = density[i][j][k + 1];
								double d8 = density[i][j + 1][k + 1];

								dnst = (d1 + d2 + d3 + d4 + d5 + d6 + d7 + d8) / avr;
							}
							//fprintf(f,"%f,",dnst);
							fprintf(ftp, "%lf ", (double)uu*(TProblem::GlobalVelocity));
							fprintf(ftp, "%lf ", (double)vv*(TProblem::GlobalVelocity));
							fprintf(ftp, "%lf ", (double)ww*(TProblem::GlobalVelocity));
							fprintf(ftp, "%lf ", (double)pp);
							fprintf(ftp, "%lf ", (double)dnst);
							fprintf(ftp, "%lf\n", (double)cc);

							//Viscosity (Begin)			
							//const T3dNumberMask& maskVis = gridVis.Mask;			
							//long double vs;
							//if (mijk==TFictivePoint) {vs=0;}
							//else
							//{					
							//bool avr1 = ((maskVis[i+1][j][k])!=TFictivePoint);
							//bool avr2 = ((maskVis[i+1][j+1][k])!=TFictivePoint);
							//bool avr3 = ((maskVis[i+1][j][k+1])!=TFictivePoint);
							//bool avr4 = ((maskVis[i+1][j+1][k+1])!=TFictivePoint);
							//bool avr5 = ((maskVis[i][j][k])!=TFictivePoint);
							//bool avr6 = ((maskVis[i][j+1][k])!=TFictivePoint);
							//bool avr7 = ((maskVis[i][j][k+1])!=TFictivePoint);
							//bool avr8 = ((maskVis[i][j+1][k+1])!=TFictivePoint);
							//int avr = (int)avr1 + (int)avr2 + (int)avr3 + (int)avr4 + (int)avr5 + (int)avr6 + (int)avr7 + (int)avr8;

							//if (avr==0)
							//{
							//throw "Error in TVPProblem::SolveBySplitting: Grid is incorrect while saving viscosity results.";  
							//}

							//double vs1 = vis[i+1][j][k];
							//double vs2 = vis[i+1][j+1][k];
							//double vs3 = vis[i+1][j][k+1];
							//double vs4 = vis[i+1][j+1][k+1];
							//double vs5 = vis[i][j][k];
							//double vs6 = vis[i][j+1][k];
							//double vs7 = vis[i][j][k+1];
							//double vs8 = vis[i][j+1][k+1];

							//vs = ( vs1 + vs2 + vs3 + vs4 + vs5  + vs6 + vs7 + vs8 )/avr;	 				
							//int stop = 0;

							//}
							//fprintf(f,"%f,",vs);	        			
							//Viscosity (End)


#if CNPY_VISUALIZATION
							U[k*(L + 1)*(N + 1) + j*(N + 1) + i] = uu*(TProblem::GlobalVelocity);
							V[k*(L + 1)*(N + 1) + j*(N + 1) + i] = vv*(TProblem::GlobalVelocity);
							W[k*(L + 1)*(N + 1) + j*(N + 1) + i] = ww*(TProblem::GlobalVelocity);
							P[k*(L + 1)*(N + 1) + j*(N + 1) + i] = pp*(TProblem::GlobalPressure);

							FX[k*(L + 1)*(N + 1) + j*(N + 1) + i] = fx*(TProblem::GlobalVelocity);
							FY[k*(L + 1)*(N + 1) + j*(N + 1) + i] = fy*(TProblem::GlobalVelocity);
							FZ[k*(L + 1)*(N + 1) + j*(N + 1) + i] = fz*(TProblem::GlobalVelocity);

							CNC[k*(L + 1)*(N + 1) + j*(N + 1) + i] = cc;
							//VIS[k*(L+1)*(N+1)+j*(N+1)+i] = vs;
							DNST[k*(L + 1)*(N + 1) + j*(N + 1) + i] = dnst;
#endif

							/*
							fprintf(f,"%f, ",(float) uu*(TProblem::GlobalVelocity));
							fprintf(f,"%f, ", (float) vv*(TProblem::GlobalVelocity));
							fprintf(f,"%f, ", (float) ww*(TProblem::GlobalVelocity));
							fprintf(f,"%f, ", (float) pp*(TProblem::GlobalPressure));
							*/
                /*
						}
					}
				}
				fclose(ftp);

				//---------------------------------------------------------------
				// Output to VTK
				//---------------------------------------------------------------
				/*
				char vtk_fileName[40];
				if (timeStepNumber == 1)
					string_print(vtk_fileName, 40, "flow%d.vtk", 1000);
				else
					string_print(vtk_fileName, 40, "flow%d.vtk", 1000 + timeStepNumber);
				FILE* vtk_ftp = file_open(vtk_fileName, "w");

				fprintf(vtk_ftp, "# vtk DataFile Version 3.0\n");
				fprintf(vtk_ftp, "Flow\nASCII\n\n");
				fprintf(vtk_ftp, "DATASET STRUCTURED_GRID\n");

				T3dNormalGrid& vtk_grid = dynamic_cast<T3dNormalGrid&>(*fGrid);
				int vtk_N = grid.GetSeparator1().EndIndex;
				int vtk_L = grid.GetSeparator2().EndIndex;
				int vtk_M = grid.GetSeparator3().EndIndex;
				int vtk_NML = (N + 1)*(L + 1)*(M + 1);
				fprintf(vtk_ftp, "DIMENSIONS %d %d %d\n", vtk_N + 1, vtk_L + 1, vtk_M + 1);
				fprintf(vtk_ftp, "POINTS %d double\n", vtk_NML);
				const T3dNumberMask& vtk_mask = grid.Mask;

				const long double* vtk_x = vtk_grid.GetSeparator1().Dimension;
				const long double* vtk_y = vtk_grid.GetSeparator2().Dimension;
				const long double* vtk_z = vtk_grid.GetSeparator3().Dimension;
				//const T3dNumberMask& vtk_mask = grid.Mask;
				for (int k = 0; k < M + 1; k++)
				for (int j = 0; j < L + 1; j++)
				for (int i = 0; i < N + 1; i++)
					fprintf(vtk_ftp, "%f  %f  %f\n", (float)vtk_x[i] * (TProblem::GlobalLength), (float)vtk_y[j] * (TProblem::GlobalLength), (float)vtk_z[k] * (TProblem::GlobalLength));

				T3dNormalGrid& vtk_gridP = dynamic_cast<T3dNormalGrid&>(*fGridP);
				const T3dNumberMask& vtk_maskP = gridP.Mask;
				fprintf(vtk_ftp, "POINT_DATA %d\n\nVECTORS Velocity double\n", vtk_NML);
				for (int k = 0; k < M + 1; k++)
				{
					for (int j = 0; j < L + 1; j++)
					{
						for (int i = 0; i < N + 1; i++)
						{

							long double uu;
							long double vv;
							long double ww;

							int mijk = vtk_mask[i][j][k];
							if (mijk == TFictivePoint) { uu = 0; vv = 0; ww = 0; }
							else
							{
								uu = (uR[i][j][k] + uR[i][j + 1][k] + uR[i][j][k + 1] + uR[i][j + 1][k + 1]) / 4;
								vv = (vR[i][j][k] + vR[i + 1][j][k] + vR[i][j][k + 1] + vR[i + 1][j][k + 1]) / 4;
								ww = (wR[i][j][k] + wR[i + 1][j][k] + wR[i][j + 1][k] + wR[i + 1][j + 1][k]) / 4;
							}

							//fprintf(f,"%f,",dnst);
							fprintf(vtk_ftp, "%lf ", (double)uu*(TProblem::GlobalVelocity));
							fprintf(vtk_ftp, "%lf ", (double)vv*(TProblem::GlobalVelocity));
							fprintf(vtk_ftp, "%lf\n ", (double)ww*(TProblem::GlobalVelocity));

						}
					}
				}

				fprintf(vtk_ftp, "\nSCALARS Pressure float\nLOOKUP_TABLE default\n");
				for (int k = 0; k < M + 1; k++)
				{
					for (int j = 0; j < L + 1; j++)
					{
						for (int i = 0; i < N + 1; i++)
						{
							long double pp;
							int mijk = vtk_maskP[i][j][k];
							if (mijk == TFictivePoint) { pp = 0; }
							else
							{
								bool avr1 = ((vtk_maskP[i + 1][j][k]) != TFictivePoint);
								bool avr2 = ((vtk_maskP[i + 1][j + 1][k]) != TFictivePoint);
								bool avr3 = ((vtk_maskP[i + 1][j][k + 1]) != TFictivePoint);
								bool avr4 = ((vtk_maskP[i + 1][j + 1][k + 1]) != TFictivePoint);
								bool avr5 = ((vtk_maskP[i][j][k]) != TFictivePoint);
								bool avr6 = ((vtk_maskP[i][j + 1][k]) != TFictivePoint);
								bool avr7 = ((vtk_maskP[i][j][k + 1]) != TFictivePoint);
								bool avr8 = ((vtk_maskP[i][j + 1][k + 1]) != TFictivePoint);
								int avr = (int)avr1 + (int)avr2 + (int)avr3 + (int)avr4 + (int)avr5 + (int)avr6 + (int)avr7 + (int)avr8;

								if (avr == 0)
								{
									throw "Error in TVPProblem::SolveBySplitting: Grid is incorrect while saving pressure results.";
								}

								double p1 = p[i + 1][j][k];
								double p2 = p[i + 1][j + 1][k];
								double p3 = p[i + 1][j][k + 1];
								double p4 = p[i + 1][j + 1][k + 1];
								double p5 = p[i][j][k];
								double p6 = p[i][j + 1][k];
								double p7 = p[i][j][k + 1];
								double p8 = p[i][j + 1][k + 1];

								pp = (p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8) / avr;
							}
							fprintf(vtk_ftp, "%lf\n ", (double)pp);
						}
					}
				}

				fprintf(vtk_ftp, "\nSCALARS Density float\nLOOKUP_TABLE default\n");
				const T3dNumberMask& vtk_maskD = gridDensity.Mask;
				for (int k = 0; k < M + 1; k++)
				{
					for (int j = 0; j < L + 1; j++)
					{
						for (int i = 0; i < N + 1; i++)
						{
							long double dnst;
							int mijk = vtk_maskD[i][j][k];
							if (mijk == TFictivePoint) { dnst = 0; }
							else
							{
								bool avr1 = ((vtk_maskD[i + 1][j][k]) != TFictivePoint);
								bool avr2 = ((vtk_maskD[i + 1][j + 1][k]) != TFictivePoint);
								bool avr3 = ((vtk_maskD[i + 1][j][k + 1]) != TFictivePoint);
								bool avr4 = ((vtk_maskD[i + 1][j + 1][k + 1]) != TFictivePoint);
								bool avr5 = ((vtk_maskD[i][j][k]) != TFictivePoint);
								bool avr6 = ((vtk_maskD[i][j + 1][k]) != TFictivePoint);
								bool avr7 = ((vtk_maskD[i][j][k + 1]) != TFictivePoint);
								bool avr8 = ((vtk_maskD[i][j + 1][k + 1]) != TFictivePoint);
								int avr = (int)avr1 + (int)avr2 + (int)avr3 + (int)avr4 + (int)avr5 + (int)avr6 + (int)avr7 + (int)avr8;

								if (avr == 0)
								{
									throw "Error in TVPProblem::SolveBySplitting: Grid is incorrect while saving concentration results.";
								}

								double d1 = density[i + 1][j][k];
								double d2 = density[i + 1][j + 1][k];
								double d3 = density[i + 1][j][k + 1];
								double d4 = density[i + 1][j + 1][k + 1];
								double d5 = density[i][j][k];
								double d6 = density[i][j + 1][k];
								double d7 = density[i][j][k + 1];
								double d8 = density[i][j + 1][k + 1];

								dnst = (d1 + d2 + d3 + d4 + d5 + d6 + d7 + d8) / avr;
							}
							fprintf(vtk_ftp, "%lf\n", (double)dnst);

						}
					}
				}

				fprintf(vtk_ftp, "\nSCALARS Concentration float\nLOOKUP_TABLE default\n");
				const T3dNumberMask& vtk_maskC = gridC.Mask;
				//	  fprintf(ftp,"ZONE I=%d, J=%d, K=%d, F=POINT\n", N+1, L+1, M+1);
				for (int k = 0; k < M + 1; k++)
				{
					for (int j = 0; j < L + 1; j++)
					{
						for (int i = 0; i < N + 1; i++)
						{
							long double cc;
							int mijk = vtk_maskC[i][j][k];
							if (mijk == TFictivePoint) { cc = 0; }
							else
							{
								bool avr1 = ((vtk_maskC[i + 1][j][k]) != TFictivePoint);
								bool avr2 = ((vtk_maskC[i + 1][j + 1][k]) != TFictivePoint);
								bool avr3 = ((vtk_maskC[i + 1][j][k + 1]) != TFictivePoint);
								bool avr4 = ((vtk_maskC[i + 1][j + 1][k + 1]) != TFictivePoint);
								bool avr5 = ((vtk_maskC[i][j][k]) != TFictivePoint);
								bool avr6 = ((vtk_maskC[i][j + 1][k]) != TFictivePoint);
								bool avr7 = ((vtk_maskC[i][j][k + 1]) != TFictivePoint);
								bool avr8 = ((vtk_maskC[i][j + 1][k + 1]) != TFictivePoint);
								int avr = (int)avr1 + (int)avr2 + (int)avr3 + (int)avr4 + (int)avr5 + (int)avr6 + (int)avr7 + (int)avr8;

								if (avr == 0)
								{
									throw "Error in TVPProblem::SolveBySplitting: Grid is incorrect while saving concentration results.";
								}

								double c1 = cnc[i + 1][j][k];
								double c2 = cnc[i + 1][j + 1][k];
								double c3 = cnc[i + 1][j][k + 1];
								double c4 = cnc[i + 1][j + 1][k + 1];
								double c5 = cnc[i][j][k];
								double c6 = cnc[i][j + 1][k];
								double c7 = cnc[i][j][k + 1];
								double c8 = cnc[i][j + 1][k + 1];

								cc = (c1 + c2 + c3 + c4 + c5 + c6 + c7 + c8) / avr;
							}

							fprintf(vtk_ftp, "%lf\n", (double)cc);
						}
					}
				}
				fclose(vtk_ftp);
			}

			//------------------------------------------------------------------
			// Two format output end
			//------------------------------------------------------------------


			CountBoundaryP(p);
#if CNPY_VISUALIZATION
			const unsigned int shape[] = { M + 1, L + 1, N + 1 };

			string_print(fileName, 40, "u_velocity_%03d.npy", timeStepNumber);
			cnpy::npy_save(fileName, U, shape, 3, "w");

			string_print(fileName, 40, "v_velocity_%03d.npy", timeStepNumber);
			cnpy::npy_save(fileName, V, shape, 3, "w");

			string_print(fileName, 40, "w_velocity_%03d.npy", timeStepNumber);
			cnpy::npy_save(fileName, W, shape, 3, "w");

			string_print(fileName, 40, "fx_%03d.npy", timeStepNumber);
			cnpy::npy_save(fileName, FX, shape, 3, "w");

			string_print(fileName, 40, "fy_%03d.npy", timeStepNumber);
			cnpy::npy_save(fileName, FY, shape, 3, "w");

			string_print(fileName, 40, "fz_%03d.npy", timeStepNumber);
			cnpy::npy_save(fileName, FZ, shape, 3, "w");


			string_print(fileName, 40, "pressure_%03d.npy", timeStepNumber);
			cnpy::npy_save(fileName, P, shape, 3, "w");

			string_print(fileName, 40, "concentration_%03d.npy", timeStepNumber);
			cnpy::npy_save(fileName, CNC, shape, 3, "w");

			//snprintf(fileName, 40, "viscosity_%03d.npy", timeStepNumber);	  
			//cnpy::npy_save(fileName, VIS, shape, 3, "w");

			string_print(fileName, 40, "density_%03d.npy", timeStepNumber);
			cnpy::npy_save(fileName, DNST, shape, 3, "w");

			const unsigned int shape_x[] = { N + 1 };
			const unsigned int shape_y[] = { L + 1 };
			const unsigned int shape_z[] = { M + 1 };

			string_print(fileName, 40, "x_coord.npy", timeStepNumber);
			cnpy::npy_save(fileName, x, shape_x, 1, "w");

			string_print(fileName, 40, "y_coord.npy", timeStepNumber);
			cnpy::npy_save(fileName, y, shape_y, 1, "w");

			string_print(fileName, 40, "z_coord.npy", timeStepNumber);
			cnpy::npy_save(fileName, z, shape_z, 1, "w");

			delete U;
			delete V;
			delete W;
			delete P;

			delete FX;
			delete FY;
			delete FZ;
#endif


			OutputBoundary(timeStepNumber);
			system("./replace.sh");

			delete &rpTrXR;
			delete &rpTrYR;
			delete &rpTrZR;

			end = clock();
			diff = (double)(end - start);
			diff = diff / CLOCKS_PER_SEC;
			printf("Saving rerults: OK. (time = %lf s.)\n", diff);
		}
		//Saving results (End)


		//Reassigning (Begin)
		{
			printf("\n");
			printf("Reassigning: Starting...\n");
			start = clock();

			u_.Assign(u);
			v_.Assign(v);
			w_.Assign(w);

			geometry.U1->field.Assign(u);
			geometry.V1->field.Assign(v);
			geometry.W1->field.Assign(w);

			//Concentration (Begin)
			cnc_1.Assign(cnc);
			//Concentration (End)

			end = clock();
			diff = (double)(end - start);
			diff = diff / CLOCKS_PER_SEC;
			printf("Reassigning: OK. (time = %lf s.)\n", diff);
		}
		//Reassigning (End)	


		flagRestruct = RestructStatement(fStatement, fU, fV, fW, fP);

		printf("\n");
		printf("Time step = %d: OK.\n", timeStepNumber);
		printf("*****************************************************\n");
	} // by time 


	//Deleting objects (Begin)
	printf("\n\n");
	printf("Deleting temporary objects: Starting...\n");

	delete geometry.U1->fieldPtr;
	delete geometry.V1->fieldPtr;
	delete geometry.W1->fieldPtr;

	delete u_Ptr;
	delete v_Ptr;
	delete w_Ptr;

	delete uSPtr;
	delete vSPtr;
	delete wSPtr;


	delete rightPartPtr;


	delete APtr;
	if (ARnPtr != NULL) { delete ARnPtr; }


	//Concentration (Begin)
	delete C_1Ptr;
	delete CPtr;
	delete gridCPtr;
	//Concentration (End)

	//Viscosity (Begin)  
	delete gridVisPtr;
	delete visPtr;
	//Viscosity (End)


	printf("Deleting temporary objects: OK.\n");
	//Deleting objects (End)

	printf("\n\n");
	printf("Splitting method: OK.\n");*/
	}
	

	

	return 0;
}


