// iflow.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "param.h"
#include "geom.h"
#include "service.h"
#include "immersedboundary.h"
#include <iostream>

using namespace std;


bool InputHole(int i, int j, int k, int N, int L, int M);
bool CylinderMask(int i, int j, int k, int L, int M, int K);
bool ConcentrationInletMask(int i, int j, int k, int L, int M, int K);
bool InitialConcentrationDistribution(int i, int j, int k, int L, int M, int K);
bool InitialViscosityDistribution(int i, int j, int k, int L, int M, int K);
bool initialDensityDistribution(int i, int j, int k, int N, int L, int M);

int _tmain(int argc, _TCHAR* argv[])
{
	// файл параметров должен быть аргументом командной строки
	// но это пока не получилось
	// поэтому

	
	char inpfile[20] = "input.txt";          // файл входных параметров
	cInputParameters inpParams(inpfile);     // заполнение входных параметров из файла
	

	cGeometry problem(inpParams.domainParams); //основна€ сетка
	inpParams.timeParams.UpdateSeparationValue(); //  временные параметры

	cImmersedBoundary immersedBoundary(&inpParams); 
	
	 

	//Vectors and operators global pointers defining (Begin)

	//Current time step (Begin)
	//Already defined as class fields
	//Current time step (End)

	//Previous time step (Begin)
	/*
	TRnRelease3dSpace* u_1Ptr = NULL;
	TRnRelease3dSpace* v_1Ptr = NULL;
	TRnRelease3dSpace* w_1Ptr = NULL;
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
	problem.setFToZero();
	cField rightPart(problem.NxP, problem.NyP, problem.NzP);// Ќаверно нужно пересоздавать дл€ U, V ,W !!!!!!!!!!!!!!!
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
	int timeStepNumber;
	long double currentTime,beginTime = inpParams.timeParams.startValue, endTime = inpParams.timeParams.endValue, dTime = inpParams.timeParams.separationValue;

	for (timeStepNumber = 1, currentTime = beginTime; currentTime <= endTime; timeStepNumber++, currentTime += dTime)
	{
		
		//ChangeStifness(timeStepNumber);
		printf("\n\n");
		printf("*****************************************************\n");
		printf("Time step = %d: Starting...\n", timeStepNumber);

		

		problem.Uold = problem.Fx;
		problem.Vold = problem.Fy;
		problem.Wold = problem.Fz;
		problem.Pold = problem.Fp;

			//Previous time step (Begin)
			problem.U1 = problem.Fx;
			problem.V1 = problem.Fy;
			problem.W1 = problem.Fz;
			//Previous time step (End)

			//Intermediate (Begin)
			problem.U = problem.Fx;
			problem.V = problem.Fy;
			problem.W = problem.Fz;
			//Intermediate (End)

			//For control stating (Begin)
			problem.Us = problem.Fx;
			problem.Vs = problem.Fy;
			problem.Ws = problem.Fz;
			//For control stating (End)
			rightPart = 0.0;




		
	

		
			problem.MaskRho = problem.MaskP; //нужно добавить переприсвоение нормалей
			for (int i = 0; i < problem.MaskRho.Nx; i++)
			{
				for (int j = 0; j < problem.MaskRho.Ny; j++)
				{
					for (int k = 0; k < problem.MaskRho.Nz; k++)
					{
						int mijk = problem.MaskRho.param[i][j][k].type;
						
						if (mijk == TPreDefinedBorderPoint)
						{
							if (i == 0 || j == 0 || k == 0)
							{
								problem.MaskRho.param[i][j][k].type = mijk;
							}
							else if (i == problem.MaskRho.Nx || j == problem.MaskRho.Ny || k == problem.MaskRho.Nz)
							{
								problem.MaskRho.param[i][j][k].type = TPreNormalBorderPoint;
							}
							else
							{
								//printf("N = %d, L = %d, M = %d\n", N, L, M);
								printf("i = %d, j = %d, k = %d\n", i, j, k);
								throw "Error in TVPProblem::SolveBySplitting: Unknown concentration PreDefined point mask type.";
							}
						}
						else if (mijk == TFictivePoint)
						{
							if (((j == 0) || (j == problem.MaskRho.Ny)) && (i != 0) && (i != problem.MaskRho.Nx) && (k != 0) && (k != problem.MaskRho.Nz))
							{
								problem.MaskRho.param[i][j][k].type = TPreNormalBorderPoint;
							}
							else if (((k == 0) || (k == problem.MaskRho.Nz)) && (i != 0) && (i != problem.MaskRho.Nx) && (j != 0) && (j != problem.MaskRho.Ny))
							{
								problem.MaskRho.param[i][j][k].type = TPreNormalBorderPoint;
							}
							else if (((i == 0) || (i == problem.MaskRho.Nx)) && (k != 0) && (k != problem.MaskRho.Nz) && (j != 0) && (j != problem.MaskRho.Ny))
							{
								problem.MaskRho.param[i][j][k].type = TPreNormalBorderPoint;
							}
							else
							{
								problem.MaskRho.param[i][j][k].type = mijk;
							}
						}
						else
						{
							throw "Error in TVPProblem::SolveBySplitting: Unknown density point mask type.";
						}
					}
				}
			}

			for (int i = 0; i < problem.MaskRho.Nx; i++)
			{
				for (int j = 0; j <problem.MaskRho.Ny; j++)
				{
					for (int k = 0; k < problem.MaskRho.Nz; k++)
					{
						int mijk = problem.MaskRho.param[i][j][k].type;
						if ((mijk == TActualPoint) || (mijk == TFictivePoint))
						{
							problem.MaskRho.param[i][j][k].type = mijk;
						}
						else if (mijk == TPreDefinedBorderPoint)
						{
							problem.MaskRho.param[i][j][k].type = TFictivePoint;
						}
						else if (mijk == TPreNormalBorderPoint)
						{
							problem.MaskRho.param[i][j][k].type = TActualPoint;
						}
						else
						{
							throw "Error in TVPProblem::SolveBySplitting: Unknown density point mask type.";
						}
					}
				}
			}


				//initial conditions (Begin)
			for (int i = 0; i < problem.MaskRho.Nx; i++)
				{
				for (int j = 0; j <problem.MaskRho.Ny; j++)
					{
					for (int k = 0; k < problem.MaskRho.Nz; k++)
						{
						problem.Rho.field[i][j][k]=inpParams.physicalParams.envDensityValue;
						}
					}
				}

			for (int i = 0; i < problem.MaskRho.Nx; i++)  
				{
				for (int j = 0; j < problem.MaskRho.Ny; j++)
					{
					for (int k = 0; k < problem.MaskRho.Nz; k++)
						{
						if (initialDensityDistribution(i, j, k, problem.MaskRho.Nx, problem.MaskRho.Ny, problem.MaskRho.Nz))
							{
								problem.Rho.field[i][j][k] = inpParams.physicalParams.soluteDensityValue;
							}
						}
					}
				}
			


			//Concentration (Begin)
				//Mask modification for grid (Begin)
				
				for (int i = 0; i<problem.MaskC.Nx; i++)
				{
					for (int j = 0; j<problem.MaskC.Ny; j++)
					{
						for (int k = 0; k<problem.MaskC.Nz; k++)
						{
							if (i == 0 && ConcentrationInletMask(i, j, k, problem.MaskC.Nx, problem.MaskC.Ny, problem.MaskC.Nz))
							{
								problem.MaskC.param[i][j][k].type = TPreDefinedBorderPoint;
								continue;
							}
							else if (!ConcentrationInletMask(i, j, k, problem.MaskC.Nx, problem.MaskC.Ny, problem.MaskC.Nz))
							{
								problem.MaskC.param[i][j][k].type = TPreNormalBorderPoint;
								continue;
							}

							int mijk = problem.MaskC.param[i][j][k].type ;
							if (mijk == TActualPoint) { problem.MaskC.param[i][j][k].type  = mijk; }
							else if (mijk == TPreDefinedBorderPoint)
							{
								if (i == 0 || j == 0 || k == 0) { problem.MaskC.param[i][j][k].type  = mijk; }
								else if (i == problem.MaskC.Nx || j == problem.MaskC.Ny || k == problem.MaskC.Nz)
								{ 
									problem.MaskC.param[i][j][k].type  = TPreNormalBorderPoint; 
								}
								else
								{
									//printf("N = %d, L = %d, M = %d\n", N, L, M);
									//printf("i = %d, j = %d, k = %d\n", i, j, k);
									throw "Error in TVPProblem::SolveBySplitting: Unknown concentration PreDefined point mask type.";
								}
							}
							else if (mijk == TFictivePoint)
							{
								if (((j == 0) || (j == problem.MaskC.Ny)) && (i != 0) && (i != problem.MaskC.Nx) && (k != 0) && (k != problem.MaskC.Nz))
								{
									problem.MaskC.param[i][j][k].type = TPreNormalBorderPoint;
								}
								else if (((k == 0) || (k == problem.MaskC.Nz)) && (i != 0) && (i != problem.MaskC.Nx) && (j != 0) && (j != problem.MaskC.Ny))
								{
									problem.MaskC.param[i][j][k].type  = TPreNormalBorderPoint;
								}
								else if (((i == 0) || (i == problem.MaskC.Nx)) && (k != 0) && (k != problem.MaskC.Nz) && (j != 0) && (j != problem.MaskC.Ny))
								{
									problem.MaskC.param[i][j][k].type  = TPreNormalBorderPoint;
								}
								else { problem.MaskC.param[i][j][k].type  = mijk; }
							}
							else
							{
								throw "Error in TVPProblem::SolveBySplitting: Unknown concentration point mask type.";
							}
						} //by k
					} //by j
				} //by i			  
				//Mask modification for grid (End)

				//Mask modification for space (Begin)
				for (int i = 0; i<problem.MaskC.Nx; i++)
				{
					for (int j = 0; j<problem.MaskC.Ny; j++)
					{
						for (int k = 0; k<problem.MaskC.Nz; k++)
						{
							int mijk = problem.MaskC.param[i][j][k].type ;
							if ((mijk == TActualPoint) || (mijk == TFictivePoint)) {problem.MaskC.param[i][j][k].type  = mijk; }
							else if (mijk == TPreDefinedBorderPoint)
							{
								problem.MaskC.param[i][j][k].type  = TFictivePoint;
							}
							else if (mijk == TPreNormalBorderPoint)
							{
								problem.MaskC.param[i][j][k].type  = TActualPoint;
							}
							else
							{
								throw "Error in TVPProblem::SolveBySplitting: Unknown concentration point mask type.";
							}
						}
					}
				}
				//Mask modification for space (End)
				//initial conditions (Begin)
				for (int i = 0; i<problem.MaskRho.Nx; i++)
				{
					for (int j = 0; j<problem.MaskRho.Ny; j++)
					{
						for (int k = 0; k<problem.MaskRho.Nz; k++)
						{
							problem.C1.field[i][j][k] = 0;
							problem.C.field[i][j][k] = 0;
						}
					}
				}

				for (int i = 0; i<problem.MaskRho.Nx; i++)
				{
					for (int j = 0; j < problem.MaskRho.Ny; j++)
					{
						for (int k = 0; k < problem.MaskRho.Nz; k++)
						{
							if (InitialConcentrationDistribution(i, j, k, problem.Rho.Nx, problem.Rho.Ny, problem.Rho.Nz))
							{
								problem.C.field[i][j][k] = 0.1;
							}
						}
					}
				}

				//initial conditions (End)		   		   
			}
			//Concentration (End)


			//Viscosity (Begin)	 
			
				//Mask modification for grid (Begin)
	
				for (int i = 0; i<problem.MaskV.Nx; i++)
				{
					for (int j = 0; j<problem.MaskV.Ny; j++)
					{
						for (int k = 0; k<problem.MaskV.Nz; k++)
						{
							int mijk = problem.MaskV.param[i][j][k].type;
							if (mijk == TActualPoint) { problem.MaskV.param[i][j][k].type = mijk; }
							else if (mijk == TPreDefinedBorderPoint) { problem.MaskV.param[i][j][k].type = TFictivePoint; }
							else if (mijk == TFictivePoint) { problem.MaskV.param[i][j][k].type = mijk; }
							else
							{
								throw "Error in TVPProblem::SolveBySplitting: Unknown viscosity point mask type.";
							}
						} //by k
					} //by j
				} //by i			  
				//Mask modification for grid (End)
				

				//Mask modification for space (Begin)
				for (int i = 0; i<problem.MaskV.Nx; i++)
				{
					for (int j = 0; j<problem.MaskV.Ny; j++)
					{
						for (int k = 0; k<problem.MaskV.Nz; k++)
						{
							int mijk = problem.MaskV.param[i][j][k].type;
							if ((mijk == TActualPoint) || (mijk == TFictivePoint)) { problem.MaskV.param[i][j][k].type = mijk; }
							else
							{
								throw "Error in TVPProblem::SolveBySplitting: Unknown viscosity point mask type.";
							}
						}
					}
				}
				//Mask modification for space (End)
			
				
				//initial conditions (Begin)
				for (int i = 0; i<problem.MaskRho.Nx; i++)
				{
					for (int j = 0; j<problem.MaskRho.Ny; j++)
					{
						for (int k = 0; k<problem.MaskRho.Nz; k++)
						{
							int mijk = problem.MaskRho.param[i][j][k].type;
							if (mijk == TActualPoint) { problem.MaskRho.param[i][j][k].type = 1 / inpParams.physicalParams.Re; }
							else if (mijk == TFictivePoint) { problem.MaskRho.param[i][j][k].type = 0; }
							else
							{
								throw "Error in TVPProblem::SolveBySplitting: Unknown viscosity point mask type.";
							}
						}
					}
				}

				for (int i = 0; i<problem.MaskRho.Nx; i++)
				{
					for (int j = 0; j < problem.MaskRho.Ny; j++)
					{
						for (int k = 0; k < problem.MaskRho.Nz; k++)
						{
							if (InitialViscosityDistribution(i, j, k, problem.MaskRho.Nx, problem.MaskRho.Ny, problem.MaskRho.Nz))
							{
								problem.MaskRho.param[i][j][k].type = inpParams.physicalParams.visValue;
							}
						}
					}
				}

				//initial conditions (End)		   		   
			
			//Viscosity (End)

		



		

		
				//Concentration (Begin)
		{
			
			
			//Transfer terms (Begin)
			cField uTr;
			cField vTr;
			cField wTr;
			cField rpTr;
			//Transfer terms (End)

			for (int eq = 1; eq<4;eq++)
			{
				printf("Velocity equation %d: Starting...\n", eq);

				//Current grid (Begin)

			/*	if (eq == 1) { gridPtr = fGridU; }
				else if (eq == 2) { gridPtr = fGridV; }
				else if (eq == 3) { gridPtr = fGridW; }
				else { throw "Error in TVPProblem::SolveBySplitting: Unknown equation number."; }*/
				

			/*	int N = grid.GetSeparator1().EndIndex;
				int L = grid.GetSeparator2().EndIndex;
				int M = grid.GetSeparator3().EndIndex;
				const T3dNumberMask& mask = grid.Mask;*/
				//Current grid (End)

				////Current vectors (Begin)
				//TRnRelease3dSpace* vecPrevPtr;
				//TRnRelease3dSpace* vecCurPtr;
				//if (eq == 1) { vecPrevPtr = &u_1; vecCurPtr = &u_; }
				//else if (eq == 2) { vecPrevPtr = &v_1; vecCurPtr = &v_; }
				//else if (eq == 3) { vecPrevPtr = &w_1; vecCurPtr = &w_; }
				//else { throw "Error in TVPProblem::SolveBySplitting: Unknown equation number."; }
				//TRnRelease3dSpace& f_1 = *vecPrevPtr;
				//TRnRelease3dSpace& f_ = *vecCurPtr;
				//Current vectors (End)

				//Transfer terms (Begin)
				cField uTr;
				cField vTr;
				cField wTr;
				cField rpTrRef;

				cField rpTr_X, rpTr_Y, rpTr_Z;///узнать чему должны быть равны значени€

				if (eq == 1) { rpTrRef = rpTr_X; }
				else if (eq == 2) { rpTrRef = rpTr_Y; }
				else if (eq == 3) { rpTrRef = rpTr_Z; }
				else { throw "Error in TVPProblem::SolveBySplitting: Unknown equation number."; }

				if (eq == 1 && timeStepNumber == 1)
				{
					/*char u_fileName[] = "u.txt";
					FILE* u_file = file_open(u_fileName, "r");
					printf("\nNonzero U:\n");*/ // про какой файл идет речь?


					for (int i = 0; i < problem.MaskRho.Nx; i++)
					{
						for (int j = 0; j < problem.MaskRho.Ny; j++)
						{
							for (int k = 0; k < problem.MaskRho.Nz; k++)
							{
								//fscanf_s(u_file, "%lf ", &u_1[i][j][k]);
								//if (u_1[i][j][k] != 0.00) printf("%lf ", u_1[i][j][k]);
							}
						}
					}
					//fclose(u_file);

				}

				cField rpTr = rpTrRef;

				//Transfer terms (End)


				//Viscosity (Begin)
				cField visTr;
				//Viscosity (End)
				cMask f_;

				//Current vector - to zero (Begin)
				//Current vector have to be calculated at all points (including border).   
				for (int i = 0;i<problem.MaskRho.Nx;i++)
				{
					for (int j = 0;j<problem.MaskRho.Ny;j++)
					{
						for (int k = 0;k<problem.MaskRho.Nz;k++)
						{
							f_.param[i][j][k] = 0;
						}
					}
				}
				//Current vector - to zero (End)


				//Convective items, right part, border conditions (Begin)
				for (int i = 0;i<N + 1;i++)
				{
					for (int j = 0;j<L + 1;j++)
					{
						for (int k = 0;k<M + 1;k++)
						{
							int mijk = mask[i][j][k];
							if (mijk == TFictivePoint)
							{
								uTr[i][j][k] = 0; vTr[i][j][k] = 0; wTr[i][j][k] = 0; rpTr[i][j][k] = 0;
								//Viscosity (Begin)
								visTr[i][j][k] = 0;
								//Viscosity (End)
							}
							else if (mijk == TActualPoint)
							{
								if (eq == 1)
								{
									uTr[i][j][k] = u_1[i][j][k];
									vTr[i][j][k] = (v_1[i][j][k] + v_1[i + 1][j][k] + v_1[i][j - 1][k] + v_1[i + 1][j - 1][k]) / 4;
									wTr[i][j][k] = (w_1[i][j][k] + w_1[i + 1][j][k] + w_1[i][j][k - 1] + w_1[i + 1][j][k - 1]) / 4;
									//Viscosity (Begin)
									visTr[i][j][k] = (vis[i][j][k] + vis[i + 1][j][k]) / 2;
									//Viscosity (End)
								}
								else if (eq == 2)
								{
									uTr[i][j][k] = (u_1[i][j][k] + u_1[i][j + 1][k] + u_1[i - 1][j][k] + u_1[i - 1][j + 1][k]) / 4;
									vTr[i][j][k] = v_1[i][j][k];
									wTr[i][j][k] = (w_1[i][j][k] + w_1[i][j + 1][k] + w_1[i][j][k - 1] + w_1[i][j + 1][k - 1]) / 4;
									//Viscosity (Begin)
									visTr[i][j][k] = (vis[i][j][k] + vis[i][j + 1][k]) / 2;
									//Viscosity (End)
								}
								else if (eq == 3)
								{
									uTr[i][j][k] = (u_1[i][j][k] + u_1[i][j][k + 1] + u_1[i - 1][j][k] + u_1[i - 1][j][k + 1]) / 4;
									vTr[i][j][k] = (v_1[i][j][k] + v_1[i][j][k + 1] + v_1[i][j - 1][k] + v_1[i][j - 1][k + 1]) / 4;
									wTr[i][j][k] = w_1[i][j][k];
									//Viscosity (Begin)
									visTr[i][j][k] = (vis[i][j][k] + vis[i][j][k + 1]) / 2;
									//Viscosity (End)
								}
								else { throw "Error in TVPProblem::SolveBySplitting: Unknown equation number."; }

								rpTr[i][j][k] = 0; //no external actions          				
							}
							else if (mijk == TDefinedBorderPoint)
							{
								uTr[i][j][k] = 0; vTr[i][j][k] = 0; wTr[i][j][k] = 0;
								//Viscosity (Begin)
								visTr[i][j][k] = 0;
								//Viscosity (End)
								long double defVal;
								if (eq == 1) { defVal = GetBorderConditionU(i, j, k, timeStepNumber); }
								else if (eq == 2) { defVal = GetBorderConditionV(i, j, k, timeStepNumber); }
								else if (eq == 3) { defVal = GetBorderConditionW(i, j, k, timeStepNumber); }
								else { throw "Error in TVPProblem::SolveBySplitting: Unknown equation number."; }
								rpTr[i][j][k] = defVal;
							}
							else if (mijk == TPreDefinedBorderPoint)
							{
								uTr[i][j][k] = 0; vTr[i][j][k] = 0; wTr[i][j][k] = 0;
								//Viscosity (Begin)
								visTr[i][j][k] = 0;
								//Viscosity (End)
								long double defVal;
								if (eq == 1) { defVal = GetBorderConditionU(i, j, k, timeStepNumber); }
								else if (eq == 2) { defVal = GetBorderConditionV(i, j, k, timeStepNumber); }
								else if (eq == 3) { defVal = GetBorderConditionW(i, j, k, timeStepNumber); }
								else { throw "Error in TVPProblem::SolveBySplitting: Unknown equation number."; }
								rpTr[i][j][k] = defVal;
							}
							else if (mijk == TNormalBorderPoint)
							{
								//throw "Error in TVPProblem::SolveBySplitting: Method is not implemented for this point border type of velocity.";

								long double nx = (grid.Normals[i][j][k])->X;
								long double ny = (grid.Normals[i][j][k])->Y;
								long double nz = (grid.Normals[i][j][k])->Z;

								uTr[i][j][k] = 0; vTr[i][j][k] = 0; wTr[i][j][k] = 0;
								//Viscosity (Begin)
								visTr[i][j][k] = 0;
								//Viscosity (End)


								long double defVal;
								if (eq == 1)
								{
									if ((nx == 0) || (ny != 0) || (nz != 0))
									{
										throw "Error in TVPProblem::SolveBySplitting: Method is not implemented for this point border type of velocity and this equation number and this normal vector.";
									}
									defVal = GetBorderConditionU(i, j, k, timeStepNumber);
								}
								else if (eq == 2)
								{
									//defVal = GetBorderConditionV(i,j,k,timeStepNumber);
									// NOTE: Fixed here for complex conditions
									throw "Error in TVPProblem::SolveBySplitting: Method is not implemented for this point border type of velocity and this equation number.";
								}
								else if (eq == 3)
								{
									//defVal = GetBorderConditionW(i,j,k,timeStepNumber);
									// NOTE: Fixed here for complex conditions
									throw "Error in TVPProblem::SolveBySplitting: Method is not implemented for this point border type of velocity and this equation number.";
								}
								else { throw "Error in TVPProblem::SolveBySplitting: Unknown equation number."; }
								rpTr[i][j][k] = defVal;
							}
							else if (mijk == TEquationBorderPoint)
							{
								throw "Error in TVPProblem::SolveBySplitting: Method is not implemented for this point border type of velocity.";

								//Viscosity (Begin)
								visTr[i][j][k] = 0;
								//Viscosity (End)

								long double nx = (grid.Normals[i][j][k])->X;
								long double ny = (grid.Normals[i][j][k])->Y;
								long double nz = (grid.Normals[i][j][k])->Z;

								const long double* hx = grid.GetSeparator1().Separation;
								const long double* hy = grid.GetSeparator2().Separation;
								const long double* hz = grid.GetSeparator2().Separation;

								T3dNormalGrid& gridP = dynamic_cast<T3dNormalGrid&>(*fGridP);
								const long double* hxP = gridP.GetSeparator1().Separation;
								const long double* hyP = gridP.GetSeparator2().Separation;
								const long double* hzP = gridP.GetSeparator3().Separation;


								long double h_1X = hx[i - 1];
								long double hX = hx[i];
								long double  hx2X = hx[i - 1] + hx[i];
								long double  hxhxhX = hx[i - 1] * hx[i] * (hx[i - 1] + hx[i]) / 2;

								long double h_1Y = hy[j - 1];
								long double hY = hy[j];
								long double hx2Y = hy[j - 1] + hy[j];
								long double hxhxhY = hy[j - 1] * hy[j] * (hy[j - 1] + hy[j]) / 2;

								long double h_1Z = hy[k - 1];
								long double hZ = hy[k];
								long double hx2Z = hy[k - 1] + hy[k];
								long double hxhxhZ = hy[k - 1] * hy[k] * (hy[k - 1] + hy[k]) / 2;


								long double eqVal = 0;

								if (eq == 1)
								{
									if ((nx == 0) || (ny != 0) || (nz != 0))
									{
										throw "Error in TVPProblem::SolveBySplitting: Method is not implemented for this point border type of velocity and this equation number and this normal vector.";
									}

									long double hxhxhUpX = hx[i] * hx[i + 1] * (hx[i] + hx[i + 1]) / 2;
									long double hxhxhDownX = hx[i - 2] * hx[i - 1] * (hx[i - 2] + hx[i - 1]) / 2;

									vTr[i][j][k] = (v_1[i][j][k] + v_1[i + 1][j][k] + v_1[i][j - 1][k] + v_1[i + 1][j - 1][k]) / 4;
									wTr[i][j][k] = (w_1[i][j][k] + w_1[i + 1][j][k] + w_1[i][j][k - 1] + w_1[i + 1][j][k - 1]) / 4;


									eqVal = eqVal - (vTr[i][j][k])*(u_1[i][j + 1][k] - u_1[i][j - 1][k]) / (hx2Y);
									eqVal = eqVal + (1 / fRe)*(hY*u_1[i][j - 1][k] - hx2Y*u_1[i][j][k] + h_1Y*u_1[i][j + 1][k]) / hxhxhY;

									eqVal = eqVal - (wTr[i][j][k])*(u_1[i][j][k + 1] - u_1[i][j][k - 1]) / (hx2Z);
									eqVal = eqVal + (1 / fRe)*(hZ*u_1[i][j][k - 1] - hx2Z*u_1[i][j][k] + h_1Z*u_1[i][j][k + 1]) / hxhxhZ;


									//eqVal = eqVal - (p[i+1][j][k] - p[i][j][k])/hxP[i];


									if (nx<0)
									{
										eqVal = eqVal - (u_1[i][j][k])*(u_1[i + 1][j][k] - u_1[i][j][k]) / (hx[i]);
										eqVal = eqVal + (1 / fRe)*(hx[i + 1] * u_1[i][j][k] - (hx[i] + hx[i + 1])*u_1[i + 1][j][k] + hx[i] * u_1[i + 2][j][k]) / hxhxhUpX;
									}
									else
									{
										eqVal = eqVal - (u_1[i][j][k])*(u_1[i][j][k] - u_1[i - 1][j][k]) / (hx[i - 1]);
										eqVal = eqVal + (1 / fRe)*(hx[i - 1] * u_1[i - 2][j][k] - (hx[i - 2] + hx[i - 1])*u_1[i - 1][j][k] + hx[i - 2] * u_1[i][j][k]) / hxhxhDownX;
									}

									eqVal = tau*eqVal + u_1[i][j][k];
								}
								else if (eq == 2)
								{
									throw "Error in TVPProblem::SolveBySplitting: Method is not implemented for this point border type of velocity and this equation number.";
								}
								else if (eq == 3)
								{
									throw "Error in TVPProblem::SolveBySplitting: Method is not implemented for this point border type of velocity and this equation number.";
								}
								else { throw "Error in TVPProblem::SolveBySplitting: Unknown equation number."; }
								uTr[i][j][k] = 0; vTr[i][j][k] = 0; wTr[i][j][k] = 0;
								rpTr[i][j][k] = eqVal;
							}
							else if ((mijk == TPreNormalBorderPoint) || (mijk == TPreEquationBorderPoint))
							{
								throw "Error in TVPProblem::SolveBySplitting: Method is not implemented for this point border type of velocity.";
							}
							else
							{
								throw "Error in TVPProblem::SolveBySplitting: Unknown point mask type.";
							}
						} //by k
					} //by j
				} //by i

				TRnRelease3dSpace& rpTrTmp = dynamic_cast<TRnRelease3dSpace &>(f_1.CreateInstance());

				const T3dNumberMask& maskC = (*gridCPtr).Mask;
				for (int i = 1;i<N;i++)
				{
					for (int j = 1;j<L;j++)
					{
						for (int k = 1;k<M;k++)
						{
							maskC[i][j][k] = TActualPoint;
						}
					}
				}

				// eq is equal to TImmersedBoundary::COORD_*
				if (eq == 1) {
					GetForce(rpTrTmp, *gridCPtr, grid, uTr, eq, timeStepNumber);
					//InterpolateUpdatePathes(grid, uTr, eq, timeStepNumber);
				}
				else if (eq == 2) {
					GetForce(rpTrTmp, *gridCPtr, grid, vTr, eq, timeStepNumber);
					//InterpolateUpdatePathes(grid, vTr, eq, timeStepNumber);
				}
				else if (eq == 3) {
					GetForce(rpTrTmp, *gridCPtr, grid, wTr, eq, timeStepNumber);
					//InterpolateUpdatePathes(grid, wTr, eq, timeStepNumber);
				}

				for (int i = 0;i<N + 1;i++)
				{
					for (int j = 0;j<L + 1;j++)
					{
						for (int k = 0;k<M + 1;k++)
						{
							int mijk = mask[i][j][k];
							if (mijk == TActualPoint)
							{
								rpTr[i][j][k] = rpTrTmp[i][j][k];
							}
						}
					}
				}

				//Convective items, right part, border conditions (End)

				//Transfer equation (Begin)
				//if (turbulenceMode==true)
				//{
				////Viscosity (Begin) 
				const long double* hx = grid.GetSeparator1().Separation;
				const long double* hy = grid.GetSeparator2().Separation;
				const long double* hz = grid.GetSeparator3().Separation;

				long double maxDVis = 0;

				for (int i = 0;i<N + 1;i++)
				{
					for (int j = 0;j<L + 1;j++)
					{
						for (int k = 0;k<M + 1;k++)
						{
							int mijk = mask[i][j][k];
							if (mijk == TActualPoint)
							{
								long double dVisX = (visTr[i + 1][j][k] - visTr[i - 1][j][k]) / (hx[i] + hx[i - 1]);
								long double dVisY = (visTr[i][j + 1][k] - visTr[i][j - 1][k]) / (hy[j] + hy[j - 1]);
								long double dVisZ = (visTr[i][j][k + 1] - visTr[i][j][k - 1]) / (hz[k] + hz[k - 1]);

								if (fabsl(dVisX)>maxDVis) { maxDVis = fabsl(dVisX); }
								if (fabsl(dVisY)>maxDVis) { maxDVis = fabsl(dVisY); }
								if (fabsl(dVisZ)>maxDVis) { maxDVis = fabsl(dVisZ); }

								long double dU = 0;
								long double dV = 0;
								long double dW = 0;

								TRnRelease3dSpace* velTrPtr = NULL;
								if (eq == 1)
								{
									velTrPtr = &uTr;
									dU = (uTr[i + 1][j][k] - uTr[i - 1][j][k]) / (hx[i] + hx[i - 1]);
									dV = (vTr[i + 1][j][k] - vTr[i - 1][j][k]) / (hx[i] + hx[i - 1]);
									dW = (wTr[i + 1][j][k] - wTr[i - 1][j][k]) / (hx[i] + hx[i - 1]);

									visTr[i][j][k] = (density[i][j][k] + density[i + 1][j][k]) / 2 * visTr[i][j][k];
								}
								else if (eq == 2)
								{
									velTrPtr = &vTr;
									dU = (uTr[i][j + 1][k] - uTr[i][j - 1][k]) / (hy[j] + hy[j - 1]);
									dV = (vTr[i][j + 1][k] - vTr[i][j - 1][k]) / (hy[j] + hy[j - 1]);
									dW = (wTr[i][j + 1][k] - wTr[i][j - 1][k]) / (hy[j] + hy[j - 1]);

									visTr[i][j][k] = (density[i][j][k] + density[i][j + 1][k]) / 2 * visTr[i][j][k];
								}
								else if (eq == 3)
								{
									velTrPtr = &wTr;
									dU = (uTr[i][j][k + 1] - uTr[i][j][k - 1]) / (hz[k] + hz[k - 1]);
									dV = (vTr[i][j][k + 1] - vTr[i][j][k - 1]) / (hz[k] + hz[k - 1]);
									dW = (wTr[i][j][k + 1] - wTr[i][j][k - 1]) / (hz[k] + hz[k - 1]);

									visTr[i][j][k] = (density[i][j][k] + density[i][j][k + 1]) / 2 * visTr[i][j][k];
								}
								else { throw "Error in TVPProblem::SolveBySplitting: Unknown equation number."; }
								TRnRelease3dSpace& velTr = *velTrPtr;

								long double dVelX = (velTr[i + 1][j][k] - velTr[i - 1][j][k]) / (hx[i] + hx[i - 1]);
								long double dVelY = (velTr[i][j + 1][k] - velTr[i][j - 1][k]) / (hy[j] + hy[j - 1]);
								long double dVelZ = (velTr[i][j][k + 1] - velTr[i][j][k - 1]) / (hz[k] + hz[k - 1]);


								rpTr[i][j][k] = rpTr[i][j][k] + density[i][j][k] * (dVisX*(dVelX + dU) + dVisY*(dVelY + dV) + dVisZ*(dVelZ + dW));
							}
						}
					}
				}

				//printf("Max vis. diff.: %lf\n",maxDVis);

			//Convective items, right part, border conditions (Begin)
		/*	for (int i = 0; i<problem.MaskRho.Nx; i++)
			{
				for (int j = 0; j<problem.MaskRho.Ny; j++)
				{
					for (int k = 0; k<problem.MaskRho.Nz; k++)
					{
						int mijk = mask[i][j][k];
						if (mijk == TFictivePoint)
						{
							uTr.field[i][j][k] = 0; vTr.field[i][j][k] = 0; wTr.field[i][j][k] = 0; rpTr.field[i][j][k] = 0;
						}
						else if (mijk == TActualPoint)
						{
							long double uu = (u[i - 1][j][k] + u[i][j][k]) / 2;
							uTr[i][j][k] = uu;

							long double vv = (v[i][j - 1][k] + v[i][j][k]) / 2;
							vTr[i][j][k] = vv;

							long double ww = (w[i][j][k - 1] + w[i][j][k]) / 2;
							wTr[i][j][k] = ww;

							rpTr[i][j][k] = 0; //no external actions
						}
						else if (mijk == TPreDefinedBorderPoint)
						{
							uTr[i][j][k] = 0; vTr[i][j][k] = 0; wTr[i][j][k] = 0;
							if (i == 0)
							{
								rpTr[i][j][k] = ConcentrationInletCondition(i, j, k, problem.MaskRho.Nx, problem.MaskRho.Ny, problem.MaskRho.Nz);
							}
							else
							{
								rpTr[i][j][k] = 0;
							}

							//rpTr[i][j][k] = 0;
						}
						else if (mijk == TPreNormalBorderPoint)
						{
							uTr[i][j][k] = 0; vTr[i][j][k] = 0; wTr[i][j][k] = 0;
							rpTr[i][j][k] = 0;
						}
						else if ((mijk == TDefinedBorderPoint) || (mijk == TNormalBorderPoint) || (mijk == TEquationBorderPoint) || (mijk == TPreEquationBorderPoint))
						{
							throw "Error in TVPProblem::SolveBySplitting: Method is not implemented for this point border type of concentration.";
						}
						else
						{
							throw "Error in TVPProblem::SolveBySplitting: Unknown concentration point mask type.";
						}
					} //by k
				} //by j
			} //by i
			//Convective items, right part, border conditions (End)
			SolveTransferEquation_(cnc, cnc_1, uTr, vTr, wTr, rpTr, gridC, 0.0);
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

			u_1.Assign(u);
			v_1.Assign(v);
			w_1.Assign(w);

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

	delete u_1Ptr;
	delete v_1Ptr;
	delete w_1Ptr;

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
	printf("Splitting method: OK.\n");
	}
	*/

	

	return 0;
}


bool initialDensityDistribution(int i, int j, int k, int N, int L, int M)
{
	return InputHole(i, j, k, N, L, M);
};


bool ConcentrationInletMask(int i, int j, int k, int L, int M, int K)
{
	return CylinderMask(i, j, k, L, M, K);
};

bool CylinderMask(int i, int j, int k, int L, int M, int K)
{
	return i == 0 && (j - L / 4 + 3)*(j - L / 4 + 3) + (k - M / 2)*(k - M / 2) < (L / 6)*(M / 6);
};

bool InitialConcentrationDistribution(int i, int j, int k, int L, int M, int K)
{
	return CylinderMask(i, j, k, L, M, K);
};

bool InitialViscosityDistribution(int i, int j, int k, int L, int M, int K)
{
	return CylinderMask(i, j, k, L, M, K);
};
