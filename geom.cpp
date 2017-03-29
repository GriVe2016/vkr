#include "stdafx.h"
#include "geom.h"
#include "param.h"
#include <math.h>
#include "immersedboundary.h"

long double CountDistribution(long double r)
{
	long double const_to = 0.3;

	if (r >= const_to)
		return 1.0;
	else
		return r / const_to;

}
 long double GetDistribution(int, int, int, long double, long double, long double)
{
	return 0.0;
}
bool CylinderMask(int i, int j, int k, int L, int M, int K)
{
	return i == 0 && (j - L / 4 + 3)*(j - L / 4 + 3) + (k - M / 2)*(k - M / 2) < (L / 6)*(M / 6);
};

bool InitialConcentrationDistribution(int i, int j, int k, int L, int M, int K)
{
	return CylinderMask(i, j, k, L, M, K);
};
bool initialDensityDistribution(int i, int j, int k, int L, int M, int K)
{
	return false;
};
bool InitialViscosityDistribution(int i, int j, int k, int L, int M, int K)
{
	return CylinderMask(i, j, k, L, M, K);
};

long double GetBorderConditionU(int lI, int lJ, int lK, int lN)
{
	return 0;
};
long double GetBorderConditionV(int lI, int lJ, int lK, int lN)
{
	return 0;
};
long double GetBorderConditionW(int lI, int lJ, int lK, int lN)
{
	return 0;
};
cGeometry::cGeometry(const sDomainParameters meshParam)
{

	btype iPart;
	int endIndex;
        btype hx,hy,hz;

	btype lengthX = meshParam.lengthX;
	btype lengthY = meshParam.lengthY;
	btype lengthZ = meshParam.lengthZ;
	iPart = 0;
	modf(lengthX / meshParam.h, &iPart);
	endIndex = (int)iPart;
    hx=  lengthX/iPart;
    int Nx = endIndex + 1;
		
	iPart = 0;
	modf(lengthY / meshParam.h, &iPart);
	endIndex = (int)iPart;
    hy=  lengthY/iPart;
	int Ny = endIndex + 1;

	iPart = 0;
	modf(lengthZ / meshParam.h, &iPart);
	endIndex = (int)iPart;
    hz=  lengthZ/iPart;
	int Nz = endIndex + 1;


		int NxU, NxV, NxW, NxP;
		int NyU, NyV, NyW, NyP;
		int NzU, NzV, NzW, NzP;
	// Размеры разнесенных сеток
		NxU = NxV = NxW = NxP = Nx + 1;
		NyU = NyV = NyW = NyP = Ny + 1;
		NzU = NzV = NzW = NzP = Nz + 1;
		--NxU;
		--NyV;
		--NzW;

		mesh = new cMesh();
		meshU = new cMesh();
		meshV = new cMesh();
		meshW = new cMesh();
		meshP = new cMesh();
    mesh->init(Nx, Ny, Nz,0.0,0.0,0.0,hx,hy,hz,lengthX,lengthY,lengthZ);
	meshU->init(NxU,NyU,NzU,0.0,-hy/2.0,-hz/2.0,hx,hy,hz,lengthX,lengthY+hy,lengthZ+hz);
	meshV->init(NxV,NyV,NzV,-hx/2.0,0.0,0.0,hx,hy,hz,lengthX+hx,lengthY,lengthZ+hz);
	meshW->init(NxW,NyW,NzW,-hx/2.0,-hy/2.0,0.0,hx,hy,hz,lengthX+hx,lengthY+hy,lengthZ);
	meshP->init(NxP,NyP,NzP,-hx/2.0,-hy/2.0,-hz/2.0,hx,hy,hz,lengthX+hx,lengthY+hy,lengthZ+hz);
	
	mask = new cMask();
	maskU = new cMask();
	maskV = new cMask();
	maskW = new cMask();
	maskP = new cMask();

	CreateGeometryGrid();
	CreateMultyGrids(NxU,NxV,NxW,NxP,NyU,NyV,NyW,NyP,NzU,NzV,NzW,NzP);

	U->init(maskU, meshU, NxU, NyU, NzU);
	U1->init(maskU, meshU, NxU, NyU, NzU);
	Us->init(maskU, meshU, NxU, NyU, NzU);
	Uprev->init(maskU, meshU, NxU, NyU, NzU);
	Uold->init(maskU, meshU, NxU, NyU, NzU);

	V->init(maskV, meshV, NxV, NyV, NzV);
	V1->init(maskV, meshV, NxV, NyV, NzV);
	Vs->init(maskV, meshV, NxV, NyV, NzV);
	Vprev->init(maskV, meshV, NxV, NyV, NzV);
	Vold->init(maskV, meshV, NxV, NyV, NzV);

	W->init(maskW, meshW, NxW, NyW, NzW);
	W1->init(maskW, meshW, NxW, NyW, NzW);
	Ws->init(maskW, meshW, NxW, NyW, NzW);
	Wprev->init(maskW, meshW, NxW, NyW, NzW);
	Wold->init(maskW, meshW, NxW, NyW, NzW);

	P->init(maskP, meshP, NxP, NyP, NzP);
	P1->init(maskP, meshP, NxP, NyP, NzP);
	Ps->init(maskP, meshP, NxP, NyP, NzP);
	Pprev->init(maskP, meshP, NxP, NyP, NzP);
	Pold->init(maskP, meshP, NxP, NyP, NzP);

	Fx->init(maskU, meshU, NxU, NyU, NzU);
	Fy->init(maskV, meshV, NxV, NyV, NzV);
	Fz->init(maskW, meshW, NxW, NyW, NzW);
	Fp->init(maskP, meshP, NxP, NyP, NzP);

	C->init(maskP, meshP, NxP, NyP, NzP);
	C1->init(maskP, meshP, NxP, NyP, NzP);
	Mu->init(maskP, meshP, NxP, NyP, NzP);
	Rho->init(maskP, meshP, NxP, NyP, NzP);


};

void cGeometry::CreateGeometryGrid()
{
	// Сначала все точки - расчетные
	for (int i = 0; i<mask->Nx; ++i)
	for (int j = 0; j<mask->Ny; ++j)
	for (int k = 0; k<mask->Nz; ++k)
	{
		mask->param[i][j][k].type = TActualPoint;
		mask->param[i][j][k].normal.set(0., 0., 0.);
	}

	// Левая и правая границы
	for (int j = 1; j < mask->Ny - 1; ++j)
	for (int k = 1; k < mask->Nz - 1; ++k)
	{
		mask->param[mask->Nx - 1][j][k].type = mask->param[0][j][k].type = TBorderPoint;
		mask->param[0][j][k].normal.set(-1., 0., 0.);
		mask->param[mask->Nx - 1][j][k].normal.set(1., 0., 0.);
	}

	// Ближняя и дальняя границы
	for (int i = 1; i<mask->Nx - 1; ++i)
	for (int k = 1; k<mask->Nz - 1; ++k)
	{
		mask->param[i][mask->Ny - 1][k].type = mask->param[i][0][k].type = TBorderPoint;
		mask->param[i][0][k].normal.set(0., -1., 0.);
		mask->param[i][mask->Ny - 1][k].normal.set(0., 1., 0.);
	}

	// Нижняя и верхняя границы
	for (int i = 1; i<mask->Nx - 1; ++i)
	for (int j = 1; j<mask->Ny - 1; ++j)
	{
		mask->param[i][j][mask->Nz - 1].type =mask->param[i][j][0].type = TBorderPoint;
		mask->param[i][j][0].normal.set(0, 0, -1);
		mask->param[i][j][mask->Nz - 1].normal.set(0, 0, 1);
	}

	// Фиктивные линии по ребрам
	for (int i = 0; i<mask->Nx; ++i)
	{
		mask->param[i][0][0].type = TFictivePoint;
		mask->param[i][0][mask->Nz - 1].type = TFictivePoint;
		mask->param[i][mask->Ny - 1][0].type = TFictivePoint;
		mask->param[i][mask->Ny - 1][mask->Nz - 1].type = TFictivePoint;
	}

	for (int j = 0; j<mask->Ny; ++j)
	{
		mask->param[0][j][0].type = TFictivePoint;
		mask->param[0][j][mask->Nz - 1].type = TFictivePoint;
		mask->param[mask->Nx - 1][j][0].type = TFictivePoint;
		mask->param[mask->Nx - 1][j][mask->Nz - 1].type = TFictivePoint;
	}

	for (int k = 0; k<mask->Nz; ++k)
	{
		mask->param[0][0][k].type = TFictivePoint;
		mask->param[0][mask->Ny - 1][k].type = TFictivePoint;
		mask->param[mask->Nx - 1][0][k].type = TFictivePoint;
		mask->param[mask->Nx - 1][mask->Ny - 1][k].type = TFictivePoint;
	}
}

void cGeometry::CreateMultyGrids(int NxU, int NxV, int NxW, int  NxP, int  NyU, int  NyV, int  NyW, int  NyP, int  NzU, int  NzV, int  NzW, int  NzP)
{
	// Заполнение маски для U //

	for (int i = 0; i<NxU; ++i)
	for (int j = 0; j<NyU; ++j)
	for (int k = 0; k<NzU; ++k)
	{
		maskU->param[i][j][k].type = TActualPoint;
		maskU->param[i][j][k].normal.set(0., 0., 0.);
	}

	// Границы области
	for (int j = 1; j<NyU - 1; ++j)
	for (int k = 1; k<NzU - 1; ++k)
	{
		maskU->param[0][j][k].type = maskU->param[NxU - 1][j][k].type = TNormalBorderPoint;
		//MaskU[0][j][k].mask = MaskU[NxU-1][j][k].mask = TEquationBorderPoint;
		maskU->param[0][j][k].normal.set(-1., 0., 0.);
		maskU->param[NxU - 1][j][k].normal.set(1., 0., 0.);
	}

	for (int i = 1; i<NxU - 1; ++i)
	for (int k = 1; k<NzU - 1; ++k)
	{
		maskU->param[i][0][k].type = maskU->param[i][NyU - 1][k].type = TPreDefinedBorderPoint;
		maskU->param[i][0][k].normal.set(0., -1., 0.);
		maskU->param[i][NyU - 1][k].normal.set(0., 1., 0.);
	}

	for (int i = 1; i<NxU - 1; ++i)
	for (int j = 1; j<NyU - 1; ++j)
	{
		maskU->param[i][j][0].type = maskU->param[i][j][NzU - 1].type = TPreDefinedBorderPoint;
		maskU->param[i][j][0].normal.set(0., 0., -1.);
		maskU->param[i][j][NzU - 1].normal.set(0., 0., 1.);
	}

	// Фиктивные линии по ребрам
	for (int i = 0; i<NxU; ++i)
	{
		maskU->param[i][0][0].type = TFictivePoint;
		maskU->param[i][0][NzU - 1].type = TFictivePoint;
		maskU->param[i][NyU - 1][0].type = TFictivePoint;
		maskU->param[i][NyU - 1][NzU - 1].type = TFictivePoint;
	}

	for (int j = 0; j<NyU; ++j)
	{
		maskU->param[0][j][0].type = TFictivePoint;
		maskU->param[0][j][NzU - 1].type = TFictivePoint;
		maskU->param[NxU - 1][j][0].type = TFictivePoint;
		maskU->param[NxU - 1][j][NzU - 1].type = TFictivePoint;
	}

	for (int k = 0; k<NzU; ++k)
	{
		maskU->param[0][0][k].type = TFictivePoint;
		maskU->param[0][NyU - 1][k].type = TFictivePoint;
		maskU->param[NxU - 1][0][k].type = TFictivePoint;
		maskU->param[NxU - 1][NyU - 1][k].type = TFictivePoint;
	}


	// Заполнение маски для V //

	for (int i = 0; i<NxV; ++i)
	for (int j = 0; j<NyV; ++j)
	for (int k = 0; k<NzV; ++k)
	{
		maskV->param[i][j][k].type = TActualPoint;
		maskV->param[i][j][k].normal.set(0., 0., 0.);
	}

	// Границы области
	for (int j = 1; j<NyV - 1; ++j)
	for (int k = 1; k<NzV - 1; ++k)
	{
		maskV->param[0][j][k].type = maskV->param[NxV - 1][j][k].type = TPreDefinedBorderPoint;
		maskV->param[0][j][k].normal.set(-1., 0., 0.);
		maskV->param[NxV - 1][j][k].normal.set(1., 0., 0.);
	}

	for (int i = 1; i<NxV - 1; ++i)
	for (int k = 1; k<NzV - 1; ++k)
	{
		maskV->param[i][0][k].type = maskV->param[i][NyV - 1][k].type = TDefinedBorderPoint;
		maskV->param[i][0][k].normal.set(0., -1., 0.);
		maskV->param[i][NyV - 1][k].normal.set(0., 1., 0.);
	}

	for (int i = 1; i<NxV - 1; ++i)
	for (int j = 1; j<NyV - 1; ++j)
	{
		maskV->param[i][j][0].type = maskV->param[i][j][NzV - 1].type = TPreDefinedBorderPoint;
		maskV->param[i][j][0].normal.set(0., 0., -1.);
		maskV->param[i][j][NzV - 1].normal.set(0., 0., 1.);
	}

	// Фиктивные линии по ребрам
	for (int i = 0; i<NxV; ++i)
	{
		maskV->param[i][0][0].type = TFictivePoint;
		maskV->param[i][0][NzV - 1].type = TFictivePoint;
		maskV->param[i][NyV - 1][0].type = TFictivePoint;
		maskV->param[i][NyV - 1][NzV - 1].type = TFictivePoint;
	}

	for (int j = 0; j<NyV; ++j)
	{
		maskV->param[0][j][0].type = TFictivePoint;
		maskV->param[0][j][NzV - 1].type = TFictivePoint;
		maskV->param[NxV - 1][j][0].type = TFictivePoint;
		maskV->param[NxV - 1][j][NzV - 1].type = TFictivePoint;
	}

	for (int k = 0; k<NzV; ++k)
	{
		maskV->param[0][0][k].type = TFictivePoint;
		maskV->param[0][NyV - 1][k].type = TFictivePoint;
		maskV->param[NxV - 1][0][k].type = TFictivePoint;
		maskV->param[NxV - 1][NyV - 1][k].type = TFictivePoint;
	}


	// Заполнение маски для W //

	for (int i = 0; i<NxW; ++i)
	for (int j = 0; j<NyW; ++j)
	for (int k = 0; k<NzW; ++k)
	{
		maskW->param[i][j][k].type = TActualPoint;
		maskW->param[i][j][k].normal.set(0., 0., 0.);
	}
	// Границы области
	for (int j = 1; j<NyW - 1; ++j)
	for (int k = 1; k<NzW - 1; ++k)
	{
		maskW->param[0][j][k].type = maskW->param[NxW - 1][j][k].type = TPreDefinedBorderPoint;
		maskW->param[0][j][k].normal.set(-1., 0., 0.);
		maskW->param[NxW - 1][j][k].normal.set(1., 0., 0.);
	}

	for (int i = 1; i<NxW - 1; ++i)
	for (int k = 1; k<NzW - 1; ++k)
	{
		maskW->param[i][0][k].type = maskW->param[i][NyW - 1][k].type = TPreDefinedBorderPoint;
		maskW->param[i][0][k].normal.set(0., -1., 0.);
		maskW->param[i][NyW - 1][k].normal.set(0., 1., 0.);
	}

	for (int i = 1; i<NxW - 1; ++i)
	for (int j = 1; j<NyW - 1; ++j)
	{
		maskW->param[i][j][0].type = maskW->param[i][j][NzW - 1].type = TDefinedBorderPoint;
		maskW->param[i][j][0].normal.set(0., 0., -1.);
		maskW->param[i][j][NzW - 1].normal.set(0., 0., 1.);
	}

	// Фиктивные линии по ребрам
	for (int i = 0; i<NxW; ++i)
	{
		maskW->param[i][0][0].type = TFictivePoint;
		maskW->param[i][0][NzW - 1].type = TFictivePoint;
		maskW->param[i][NyW - 1][0].type = TFictivePoint;
		maskW->param[i][NyW - 1][NzW - 1].type = TFictivePoint;
	}

	for (int j = 0; j<NyW; ++j)
	{
		maskW->param[0][j][0].type = TFictivePoint;
		maskW->param[0][j][NzW - 1].type = TFictivePoint;
		maskW->param[NxW - 1][j][0].type = TFictivePoint;
		maskW->param[NxW - 1][j][NzW - 1].type = TFictivePoint;
	}

	for (int k = 0; k<NzW; ++k)
	{
		maskW->param[0][0][k].type = TFictivePoint;
		maskW->param[0][NyW - 1][k].type = TFictivePoint;
		maskW->param[NxW - 1][0][k].type = TFictivePoint;
		maskW->param[NxW - 1][NyW - 1][k].type = TFictivePoint;
	}


	// Заполнение маски для P
	for (int i = 0; i<NxP; ++i)
	for (int j = 0; j<NyP; ++j)
	for (int k = 0; k<NzP; ++k)
	{
		maskP->param[i][j][k].type = TActualPoint;
		maskP->param[i][j][k].normal.set(0., 0., 0.);
	}

	// Границы области
	for (int j = 1; j<NyP - 1; ++j)
	for (int k = 1; k<NzP - 1; ++k)
	{
		maskP->param[0][j][k].type = maskP->param[NxP - 1][j][k].type = TPreDefinedBorderPoint;
		maskP->param[0][j][k].normal.set(-1., 0., 0.);
		maskP->param[NxP - 1][j][k].normal.set(1., 0., 0.);
	}

	for (int i = 1; i<NxP - 1; ++i)
	for (int k = 1; k<NzP - 1; ++k)
	{
		maskP->param[i][0][k].type = maskP->param[i][NyP - 1][k].type = TFictivePoint;
		maskP->param[i][0][k].normal.set(0., -1., 0.);
		maskP->param[i][NyP - 1][k].normal.set(0., 1., 0.);
	}

	for (int i = 1; i<NxP - 1; ++i)
	for (int j = 1; j<NyP - 1; ++j)
	{
		maskP->param[i][j][0].type = maskP->param[i][j][NzP - 1].type = TFictivePoint;
		maskP->param[i][j][0].normal.set(0., 0., -1.);
		maskP->param[i][j][NzP - 1].normal.set(0., 0., 1.);
	}

	// Фиктивные линии по ребрам
	for (int i = 0; i<NxP; ++i)
	{
		maskP->param[i][0][0].type = TFictivePoint;
		maskP->param[i][0][NzP - 1].type = TFictivePoint;
		maskP->param[i][NyP - 1][0].type = TFictivePoint;
		maskP->param[i][NyP - 1][NzP - 1].type = TFictivePoint;
	}

	for (int j = 0; j<NyP; ++j)
	{
		maskP->param[0][j][0].type = TFictivePoint;
		maskP->param[0][j][NzP - 1].type = TFictivePoint;
		maskP->param[NxP - 1][j][0].type = TFictivePoint;
		maskP->param[NxP - 1][j][NzP - 1].type = TFictivePoint;
	}

	for (int k = 0; k<NzP; ++k)
	{
		maskP->param[0][0][k].type = TFictivePoint;
		maskP->param[0][NyP - 1][k].type = TFictivePoint;
		maskP->param[NxP - 1][0][k].type = TFictivePoint;
		maskP->param[NxP - 1][NyP - 1][k].type = TFictivePoint;
	}
}


void cGeometry::setFToZero()
{
	for (int i = 0; i <Fx->Nx; i++)
		for (int j = 0; j < Fx->Ny + 1; j++)
			for (int k = 0; k < Fx->Nz + 1; k++)
				Fx->field[i][j][k] = 0.00;
	for (int i = 0; i < Fy->Nx + 1; i++)
		for (int j = 0; j < Fy->Ny; j++)
			for (int k = 0; k < Fy->Nz + 1; k++)
				Fy->field[i][j][k] = 0.00;
	for (int i = 0; i < Fz->Nx + 1; i++)
		for (int j = 0; j <  Fz->Ny + 1; j++)
			for (int k = 0; k <  Fz->Nz; k++)
				Fz->field[i][j][k] = 0.00;

	Uold = Fx;
	Vold = Fy;
	Wold = Fz;
	Pold = Fp;

	//Previous time step (Begin)
	U1 = Fx;
	V1 = Fy;
	W1 = Fz;
	//Previous time step (End)

	//Intermediate (Begin)
	U = Fx;
	V = Fy;
	W = Fz;
	//Intermediate (End)

	//For control stating (Begin)
	Us = Fx;
	Vs = Fy;
	Ws = Fz;
	//For control stating (End)
	
}
void cGeometry::SetInitCondDensity()
{
	cField densityPtr;
	if (densityPtr.field != NULL)
	{
		throw "Error in TVPProblem::SolveBySplitting: Density does not support restructing.";
	}


	

	//Mask modification for grid (Begin)
	int N = maskP->Nx;
	int L = maskP->Ny;
	int M = maskP->Nz;
	//printf("\nDENSITY:   N=%d, L=%d, M=%d\n",N,L,M);

	for (int i = 0; i < N + 1; i++)
	{
		for (int j = 0; j < L + 1; j++)
		{
			for (int k = 0; k < M + 1; k++)
			{
				int mijk = maskP->param[i][j][k].type;
				if (mijk == TActualPoint)
				{
					maskP->param[i][j][k].type = mijk;
				}
				else if (mijk == TPreDefinedBorderPoint)
				{
					if (i == 0 || j == 0 || k == 0)
					{
						maskP->param[i][j][k].type = mijk;
					}
					else if (i == N || j == L || k == M)
					{
						maskP->param[i][j][k].type = TPreNormalBorderPoint;
					}
					else
					{
						printf("N = %d, L = %d, M = %d\n", N, L, M);
						printf("i = %d, j = %d, k = %d\n", i, j, k);
						throw "Error in TVPProblem::SolveBySplitting: Unknown concentration PreDefined point mask type.";
					}
				}
				else if (mijk == TFictivePoint)
				{
					if (((j == 0) || (j == L)) && (i != 0) && (i != N) && (k != 0) && (k != M))
					{
						maskP->param[i][j][k].type = TPreNormalBorderPoint;
					}
					else if (((k == 0) || (k == M)) && (i != 0) && (i != N) && (j != 0) && (j != L))
					{
						maskP->param[i][j][k].type = TPreNormalBorderPoint;
					}
					else if (((i == 0) || (i == N)) && (k != 0) && (k != M) && (j != 0) && (j != L))
					{
						maskP->param[i][j][k].type = TPreNormalBorderPoint;
					}
					else
					{
						maskP->param[i][j][k].type = mijk;
					}
				}
				else
				{
					throw "Error in TVPProblem::SolveBySplitting: Unknown density point mask type.";
				}
			}
		}
	}

	//gridDensityPtr = new T3dNormalGrid(sep1, sep2, sep3, mask, normals, false);

	for (int i = 0; i < N + 1; i++)
	{
		for (int j = 0; j < L + 1; j++)
		{
			for (int k = 0; k < M + 1; k++)
			{
				int mijk = maskP->param[i][j][k].type;
				if ((mijk == TActualPoint) || (mijk == TFictivePoint))
				{
					maskP->param[i][j][k].type = mijk;
				}
				else if (mijk == TPreDefinedBorderPoint)
				{
					maskP->param[i][j][k].type = TFictivePoint;
				}
				else if (mijk == TPreNormalBorderPoint)
				{
					maskP->param[i][j][k].type = TActualPoint;
				}
				else
				{
					throw "Error in TVPProblem::SolveBySplitting: Unknown density point mask type.";
				}
			}
		}
	}

	/*TMask** masks = new TMask*[1];
	masks[0] = &mask;

	int spaceModelName = 88888;

	TRnRelease3dSpaceModel& spaceModel = *(new TRnRelease3dSpaceModel(1, masks));
	TRnSpace::AddModel(spaceModelName, &spaceModel, false);

	delete &mask;
	delete[] masks;
	delete &spaceModel;*/

	//densityPtr = new TRnRelease3dSpace(spaceModelName);

	//(*densityPtr) |= 0;
	cArray hx=meshP->X;
	cArray hy= meshP->Y;
	cArray hz= meshP->Z;


	long double hX_ = hx[1], hY_ = hy[1], hZ_ = hz[1];
	//printf("\nhX = %lf, hY = %lf, hZ  = %lf\n", (double)hX_, (double)hY_, (double) hZ_);
	long double r_, distr;

	long double envDensityValue = 1.0;
	//Initial conditions (Begin)
	for (int i = 0; i < N + 1; i++)
	{
		for (int j = 0; j < L + 1; j++)
		{
			for (int k = 0; k < M + 1; k++)
			{
				densityPtr.field[i][j][k] = envDensityValue;
			}
		}
	}
	long double soluteDensityValue = 1.5;
	for (int i = 0; i < N + 1; i++)
	{
		for (int j = 0; j < L; j++)
		{
			for (int k = 0; k < M; k++)
			{
				if (i == 0 /*&& InitialDensityDistribution(i, j, k, N, L, M)*/)
				{
					r_ = GetDistribution(i, j, k, hX_, hY_, hZ_);
					distr = CountDistribution(r_);
					densityPtr.field[i][j][k] = distr*soluteDensityValue + (1.00 - distr)*envDensityValue;
				}
				else
				{
					r_ = GetDistribution(i, j, k, hX_, hY_, hZ_);
					distr = CountDistribution(r_);
					densityPtr.field[i][j][k] = distr*soluteDensityValue + (1.00 - distr)*envDensityValue;
				}
			}
		}
	}
}

void cGeometry::SetInitCondVis()
{
	//Viscosity (Begin)	 
	cField visPtr;
		if ((visPtr.field != NULL))
		{
			throw "Error in TVPProblem::SolveBySplitting: Viscosity does not support restructing.";
		}

		//Mask modification for grid (Begin)
		int N = maskV->Nx;
		int L = maskV->Ny;
		int M = maskV->Nz;
		for (int i = 0;i<N + 1;i++)
		{
			for (int j = 0;j<L + 1;j++)
			{
				for (int k = 0;k<M + 1;k++)
				{
					int mijk = maskV->param[i][j][k].type;
					if (mijk == TActualPoint) { maskV->param[i][j][k].type = mijk; }
					else if (mijk == TPreDefinedBorderPoint) { maskV->param[i][j][k].type = TFictivePoint; }
					else if (mijk == TFictivePoint) { maskV->param[i][j][k].type = mijk; }
					else
					{
						throw "Error in TVPProblem::SolveBySplitting: Unknown viscosity point mask type.";
					}
				} //by k
			} //by j
		} //by i			  
		  //Mask modification for grid (End)
		//gridVisPtr = new T3dNormalGrid(sep1, sep2, sep3, mask, normals, false);

		//Mask modification for space (Begin)
		for (int i = 0;i<N + 1;i++)
		{
			for (int j = 0;j<L + 1;j++)
			{
				for (int k = 0;k<M + 1;k++)
				{
					int mijk = maskV->param[i][j][k].type;
					if ((mijk == TActualPoint) || (mijk == TFictivePoint)) { maskV->param[i][j][k].type = mijk; }
					else
					{
						throw "Error in TVPProblem::SolveBySplitting: Unknown viscosity point mask type.";
					}
				}
			}
		}
		//Mask modification for space (End)

		/*TMask** masks = new TMask*[1];
		masks[0] = &mask;

		int spaceModelName = 77777;

		TRnRelease3dSpaceModel& spaceModel = *(new TRnRelease3dSpaceModel(1, masks));
		TRnSpace::AddModel(spaceModelName, &spaceModel, false);

		delete &mask;
		delete[] masks;
		delete &spaceModel;

		visPtr = new TRnRelease3dSpace(spaceModelName);

		(*visPtr) |= 0;*/
		long double fRe=1;
		//Initial conditions (Begin)
		for (int i = 0;i<N + 1;i++)
		{
			for (int j = 0;j<L + 1;j++)
			{
				for (int k = 0;k<M + 1;k++)
				{
					int mijk = maskV->param[i][j][k].type;
					if (mijk == TActualPoint) { visPtr.field[i][j][k] = 1 / fRe; }
					else if (mijk == TFictivePoint) { visPtr.field[i][j][k] = 0; }
					else
					{
						throw "Error in TVPProblem::SolveBySplitting: Unknown viscosity point mask type.";
					}
				}
			}
		}

		cArray hx=meshV->X;
		cArray hy=meshV->Y;
		cArray hz= meshV->Z;

		long double hX_ = hx[1], hY_ = hy[1], hZ_ = hz[1];
		long double r_, distr = 0.;
		long double visValue = 1.5 / fRe;
		for (int i = 0;i<N + 1;i++)
		{
			for (int j = 0; j < L; j++)
			{
				for (int k = 0; k < M; k++)
				{
					if (i == 0 /*&& InitialDensityDistribution(i, j, k, N, L, M)*/)
					{
						r_ = GetDistribution(i, j, k, hX_, hY_, hZ_);
						distr = CountDistribution(r_);
						visPtr.field[i][j][k] = distr*visValue + (1.00 - distr) * 1 / fRe;
					}
					else
					{
						r_ = GetDistribution(i, j, k, hX_, hY_, hZ_);
						distr = CountDistribution(r_);
						visPtr.field[i][j][k] = distr*visValue + (1.00 - distr) * 1 / fRe;
					}
				}
			}
		}

		//Initial conditions (End)		   		   
}
	//Viscosity (End)




void cGeometry::SetInitCondRho(cInputParameters inpParams)
{
	
}


void cGeometry::SetInitCondC(cInputParameters inpParams)
{
	//Concentration (Begin)
	cField C_1Ptr, CPtr;
		if ((C_1Ptr.field != NULL) || (CPtr.field != NULL))
		{
			throw "Error in TVPProblem::SolveBySplitting: Concentration does not support restructing.";
		}

		//Mask modification for grid (Begin)
		int N = maskC->Nx;
		int L = maskC->Ny;
		int M = maskC->Nz;
		for (int i = 0;i<N + 1;i++)
		{
			for (int j = 0;j<L + 1;j++)
			{
				for (int k = 0;k<M + 1;k++)
				{
					if (i == 0 /*&& ConcentrationInletMask(i, j, k, N, L, M)*/)
					{
						maskC->param[i][j][k].type = TPreDefinedBorderPoint;
						continue;
					}
					else if (!false/*ConcentrationInletMask(i, j, k, N, L, M)*/)
					{
						maskC->param[i][j][k].type = TPreNormalBorderPoint;
						continue;
					}

					int mijk = maskC->param[i][j][k].type;
					if (mijk == TActualPoint) { maskC->param[i][j][k].type = mijk; }
					else if (mijk == TPreDefinedBorderPoint)
					{
						if (i == 0 || j == 0 || k == 0) { maskC->param[i][j][k].type = mijk; }
						else if (i == N || j == L || k == M) { maskC->param[i][j][k].type = TPreNormalBorderPoint; }
						else
						{
							printf("N = %d, L = %d, M = %d\n", N, L, M);
							printf("i = %d, j = %d, k = %d\n", i, j, k);
							throw "Error in TVPProblem::SolveBySplitting: Unknown concentration PreDefined point mask type.";
						}
					}
					else if (mijk == TFictivePoint)
					{
						if (((j == 0) || (j == L)) && (i != 0) && (i != N) && (k != 0) && (k != M))
						{
							maskC->param[i][j][k].type = TPreNormalBorderPoint;
						}
						else if (((k == 0) || (k == M)) && (i != 0) && (i != N) && (j != 0) && (j != L))
						{
							maskC->param[i][j][k].type = TPreNormalBorderPoint;
						}
						else if (((i == 0) || (i == N)) && (k != 0) && (k != M) && (j != 0) && (j != L))
						{
							maskC->param[i][j][k].type = TPreNormalBorderPoint;
						}
						else { maskC->param[i][j][k].type = mijk; }
					}
					else
					{
						throw "Error in TVPProblem::SolveBySplitting: Unknown concentration point mask type.";
					}
				} //by k
			} //by j
		} //by i			  
		  //Mask modification for grid (End)
		//gridCPtr = new T3dNormalGrid(sep1, sep2, sep3, mask, normals, false);

		//Mask modification for space (Begin)
		for (int i = 0;i<N + 1;i++)
		{
			for (int j = 0;j<L + 1;j++)
			{
				for (int k = 0;k<M + 1;k++)
				{
					int mijk = maskC->param[i][j][k].type;
					if ((mijk == TActualPoint) || (mijk == TFictivePoint)) { maskC->param[i][j][k].type = mijk; }
					else if (mijk == TPreDefinedBorderPoint)
					{
						maskC->param[i][j][k].type = TFictivePoint;
					}
					else if (mijk == TPreNormalBorderPoint)
					{
						maskC->param[i][j][k].type = TActualPoint;
					}
					else
					{
						throw "Error in TVPProblem::SolveBySplitting: Unknown concentration point mask type.";
					}
				}
			}
		}
		//Mask modification for space (End)

		/*TMask** masks = new TMask*[1];
		masks[0] = &mask;

		int spaceModelName = 55555;

		TRnRelease3dSpaceModel& spaceModel = *(new TRnRelease3dSpaceModel(1, masks));
		TRnSpace::AddModel(spaceModelName, &spaceModel, false);

		delete &mask;
		delete[] masks;
		delete &spaceModel;

		C_1Ptr = new TRnRelease3dSpace(spaceModelName);
		CPtr = new TRnRelease3dSpace(spaceModelName);

		(*C_1Ptr) |= 0;
		(*CPtr) |= 0;*/

		//Initial conditions (Begin)
		for (int i = 0;i<N + 1;i++)
		{
			for (int j = 0;j<L + 1;j++)
			{
				for (int k = 0;k<M + 1;k++)
				{
					C_1Ptr.field[i][j][k] = 0;
					CPtr.field[i][j][k] = 0;
				}
			}
		}
		
		cArray hx = mesh->X;
		cArray hy = mesh->Y;
		cArray hz = mesh->Z;

		long double hX_ = hx[1], hY_ = hy[1], hZ_ = hz[1];
		long double r_, distr = 0.;
		long double soluteConcentratonValue = 0.6;
		for (int i = 0;i<N + 1;i++)
		{
			for (int j = 0; j < L; j++)
			{
				for (int k = 0; k < M; k++)
				{
					if (i == 0 && InitialConcentrationDistribution(i, j, k, N, L, M))
					{
						r_ = GetDistribution(i, j, k, hX_, hY_, hZ_);
						distr = CountDistribution(r_);
						//printf("Distr: i,j,k=%d,%d,%d = %lf ", i, j, k, distr);
						C_1Ptr.field[i][j][k] = soluteConcentratonValue*distr;
					}
					else
					{
						r_ = GetDistribution(i, j, k, hX_, hY_, hZ_);
						distr = CountDistribution(r_);
						//printf("Distr: i,j,k=%d,%d,%d = %lf ", i, j, k, distr);
						C_1Ptr.field[i][j][k] = soluteConcentratonValue*distr;
					}
				}
			}
		}

		//Initial conditions (End)		   		   
}
	//Concentration (End)

void cGeometry::SolveByScalarRunning(int lM, long double* lU, const long double* lA, const long double* lB, const long double* lC, const long double* lF)
{
	int M = lM;
	long double* u = lU;
	const long double* a = lA;
	const long double* b = lB;
	const long double* c = lC;
	const long double* f = lF;

	long double *A = new long double[M + 1];
	long double *B = new long double[M + 1];

	A[0] = 0;
	B[0] = u[0];
	A[M] = B[M] = 0;

	for (int m = 1; m<M; m++)
	{
		A[m] = (-c[m]) / (a[m] * A[m - 1] + b[m]);
		B[m] = (f[m] - a[m] * B[m - 1]) / (a[m] * A[m - 1] + b[m]);
	}

	for (int m = M - 1; m>0; m--)
	{
		u[m] = A[m] * u[m + 1] + B[m];
	}

	delete[] A;
	delete[] B;
}





void cGeometry::SolveTransferEquation__(cField& lVecCur, cField& lVecPrev,cField&  lU, cField& lV, cField& lW, cField& lRP,const cMesh& lGrid, cField& lDiffusionCoefficient, int lSchemeType)
{
	//Preparing variables (Begin)
	cField& f = lVecCur;
	cField& f_1 = lVecPrev;
	cField& u = lU;
	cField& v = lV;
	cField& w = lW;
	cField& rp = lRP;

	cField& f_1_3 = lVecPrev;
	cField& f_2_3 = lVecPrev;

	const cMesh& grid = lGrid;

	

	int N=mask->Nx;
	int L=mask->Ny;
	int M=mask->Nz;

	cArray hx=mesh->X;
	cArray hy=mesh->Y;
	cArray hz=mesh->Z;

	long double tau=0.3;// = fTimeSeparator->SeparationValue;

	cField& diffCoefVar = lDiffusionCoefficient;

	//Preparing variables (End)

	//Direction order (Begin)
	int directionOrder = 1;
	int directionX;
	int directionY;
	int directionZ;
	if (directionOrder == 1) { directionX = 1;directionY = 2;directionZ = 3; }
	else if (directionOrder == 2) { directionX = 1;directionY = 3;directionZ = 2; }
	else if (directionOrder == 3) { directionX = 2;directionY = 1;directionZ = 3; }
	else if (directionOrder == 4) { directionX = 2;directionY = 3;directionZ = 1; }
	else if (directionOrder == 5) { directionX = 3;directionY = 1;directionZ = 2; }
	else if (directionOrder == 6) { directionX = 3;directionY = 2;directionZ = 1; }
	else { throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown direction order."; }
	//Direction order (End)

	//Terms using (Begin)
	int flagConvectiveX = 1; int flagDiffusionX = 1;
	int flagConvectiveY = 1; int flagDiffusionY = 1;
	int flagConvectiveZ = 1; int flagDiffusionZ = 1;
	//Terms using (Begin)

	//Scheme type (Begin)
	//0 - central;
	//1 - upstream;
	int schemeType = lSchemeType;
	//Scheme type (End)

	for (int directionGlobal = 1;directionGlobal < 4;directionGlobal++) //splitting scheame steps
	{
		for (int directionFake = 1;directionFake < 4;directionFake++) //for border condtions recalculating
		{
			int direction;
			if (directionFake == 1) { direction = directionGlobal; }
			else
			{
				if (directionGlobal == 1) { direction = directionFake; }
				else if (directionGlobal == 2)
				{
					if (directionFake == 2) { direction = 1; }
					else if (directionFake == 3) { direction = 3; }
					else { throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown fake direction."; }
				}
				else if (directionGlobal == 3) { direction = directionFake - 1; }
				else { throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown global direction."; }
			}


			int dir1EndIndex;
			int dir2EndIndex;
			int dirActEndIndex;

			if (direction == directionX)
			{
				dir1EndIndex = M;
				dir2EndIndex = L;
				dirActEndIndex = N;
			}
			else if (direction == directionY)
			{
				dir1EndIndex = M;
				dir2EndIndex = N;
				dirActEndIndex = L;
			}
			else if (direction == directionZ)
			{
				dir1EndIndex = L;
				dir2EndIndex = N;
				dirActEndIndex = M;
			}
			else { throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown direction."; }

			int flagConvective = 0;
			if (
				((direction == directionX) && (flagConvectiveX == 1)) ||
				((direction == directionY) && (flagConvectiveY == 1)) ||
				((direction == directionZ) && (flagConvectiveZ == 1))
				) {
				flagConvective = 1;
			}

			int flagDiffusion = 0;
			if (
				((direction == directionX) && (flagDiffusionX == 1)) ||
				((direction == directionY) && (flagDiffusionY == 1)) ||
				((direction == directionZ) && (flagDiffusionZ == 1))
				) {
				flagDiffusion = 1;
			}


			for (int dir1 = 0;dir1 < dir1EndIndex + 1;dir1++)
			{
				for (int dir2 = 0;dir2 < dir2EndIndex + 1;dir2++)
				{
					int dirActSecBeg = -1;
					int maskActSecBeg;

					int dirActSecEnd;
					int maskActSecEnd;

					for (int dirActSec = 0; dirActSec < dirActEndIndex + 1; dirActSec++)
					{
						int maskActSec;
						if (direction == directionX) { maskActSec = mask->param[dirActSec][dir2][dir1].type; }
						else if (direction == directionY) { maskActSec = mask->param[dir2][dirActSec][dir1].type; }
						else if (direction == directionZ) { maskActSec = mask->param[dir2][dir1][dirActSec].type; }
						else { throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown direction."; }

						if (dirActSecBeg < 0)
						{
							if (maskActSec == TActualPoint)
							{
								if (dirActSec == 0)
								{
									throw "Error in THydrodynamicsProblem::SolveTransferEquation: Initial border not found.";
								}
								else
								{
									dirActSecBeg = dirActSec - 1;
									if (direction == directionX) { maskActSecBeg = mask->param[dirActSecBeg][dir2][dir1].type; }
									else if (direction == directionY) { maskActSecBeg = mask->param[dir2][dirActSecBeg][dir1].type; }
									else if (direction == directionZ) { maskActSecBeg = mask->param[dir2][dir1][dirActSecBeg].type; }
									else { throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown direction."; }
								}
							}
						}
						else if (maskActSec != TActualPoint)
						{
							dirActSecEnd = dirActSec;
							maskActSecEnd = maskActSec;

							//Сhecking sector (Begin)
							if ((maskActSecBeg == TActualPoint) || (maskActSecBeg == TFictivePoint) || (maskActSecEnd == TActualPoint) || (maskActSecEnd == TFictivePoint))
							{
								throw "Error in THydrodynamicsProblem::SolveTransferEquation: Incorrect sector border.";
							}

							if ((dirActSecBeg < 0) || (dirActSecEnd - dirActSecBeg) < 2)
							{
								throw "Error in THydrodynamicsProblem::SolveTransferEquation: Incorrect sector size.";
							}
							//Сhecking sector (End)



							//Only borders recalculation and continue (Begin)
							if (direction != directionGlobal)
							{
								cField* f_AssignPtr;
								if (directionGlobal == 1) { f_AssignPtr = &f_1_3; }
								else if (directionGlobal == 2) { f_AssignPtr = &f_2_3; }
								else if (directionGlobal == 3) { f_AssignPtr = &f; }
								else { throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown global direction."; }
								cField& f_Assign = *f_AssignPtr;

								//Sector begin (Begin)
								if ((maskActSecBeg == TDefinedBorderPoint) || (maskActSecBeg == TEquationBorderPoint))
								{
									if (direction == directionX) { f_Assign.field[dirActSecBeg][dir2][dir1] = rp.field[dirActSecBeg][dir2][dir1]; }
									else if (direction == directionY) { f_Assign.field[dir2][dirActSecBeg][dir1] = rp.field[dir2][dirActSecBeg][dir1]; }
									else if (direction == directionZ) { f_Assign.field[dir2][dir1][dirActSecBeg] = rp.field[dir2][dir1][dirActSecBeg]; }
									else { throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown direction."; }
								}
								else if (maskActSecBeg == TPreDefinedBorderPoint)
								{
									if (direction == directionX) { f_Assign.field[dirActSecBeg][dir2][dir1] = 2 * rp.field[dirActSecBeg][dir2][dir1] - f_Assign.field[dirActSecBeg + 1][dir2][dir1]; }
									else if (direction == directionY) { f_Assign.field[dir2][dirActSecBeg][dir1] = 2 * rp.field[dir2][dirActSecBeg][dir1] - f_Assign.field[dir2][dirActSecBeg + 1][dir1]; }
									else if (direction == directionZ) { f_Assign.field[dir2][dir1][dirActSecBeg] = 2 * rp.field[dir2][dir1][dirActSecBeg] - f_Assign.field[dir2][dir1][dirActSecBeg + 1]; }
									else { throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown direction."; }
								}
								else if ((maskActSecBeg == TNormalBorderPoint) || (maskActSecBeg == TPreNormalBorderPoint))
								{
									//! need to check acceptable normal direction					   

									if (direction == directionX) { f_Assign.field[dirActSecBeg][dir2][dir1] = -(hx[dirActSecBeg])*rp.field[dirActSecBeg][dir2][dir1] + f_Assign.field[dirActSecBeg + 1][dir2][dir1]; }
									else if (direction == directionY) { f_Assign.field[dir2][dirActSecBeg][dir1] = -(hy[dirActSecBeg])*rp.field[dir2][dirActSecBeg][dir1] + f_Assign.field[dir2][dirActSecBeg + 1][dir1]; }
									else if (direction == directionZ) { f_Assign.field[dir2][dir1][dirActSecBeg] = -(hz[dirActSecBeg])*rp.field[dir2][dir1][dirActSecBeg] + f_Assign.field[dir2][dir1][dirActSecBeg + 1]; }
									else { throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown direction."; }
								}
								else { throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown border point type of begin sector."; }
								//Sector begin (End)


								//Sector end (Begin)
								if ((maskActSecEnd == TDefinedBorderPoint) || (maskActSecEnd == TEquationBorderPoint))
								{
									if (direction == directionX) { f_Assign.field[dirActSecEnd][dir2][dir1] = rp.field[dirActSecEnd][dir2][dir1]; }
									else if (direction == directionY) { f_Assign.field[dir2][dirActSecEnd][dir1] = rp.field[dir2][dirActSecEnd][dir1]; }
									else if (direction == directionZ) { f_Assign.field[dir2][dir1][dirActSecEnd] = rp.field[dir2][dir1][dirActSecEnd]; }
									else { throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown direction."; }
								}
								else if (maskActSecEnd == TPreDefinedBorderPoint)
								{
									if (direction == directionX) { f_Assign.field[dirActSecEnd][dir2][dir1] = 2 * rp.field[dirActSecEnd][dir2][dir1] - f_Assign.field[dirActSecEnd - 1][dir2][dir1]; }
									else if (direction == directionY) { f_Assign.field[dir2][dirActSecEnd][dir1] = 2 * rp.field[dir2][dirActSecEnd][dir1] - f_Assign.field[dir2][dirActSecEnd - 1][dir1]; }
									else if (direction == directionZ) { f_Assign.field[dir2][dir1][dirActSecEnd] = 2 * rp.field[dir2][dir1][dirActSecEnd] - f_Assign.field[dir2][dir1][dirActSecEnd - 1]; }
									else { throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown direction."; }
								}
								else if ((maskActSecEnd == TNormalBorderPoint) || (maskActSecEnd == TPreNormalBorderPoint))
								{
									//! need to check acceptable normal direction					   

									if (direction == directionX) { f_Assign.field[dirActSecEnd][dir2][dir1] = (hx[dirActSecEnd - 1])*rp.field[dirActSecEnd][dir2][dir1] + f_Assign.field[dirActSecEnd - 1][dir2][dir1]; }
									else if (direction == directionY) { f_Assign.field[dir2][dirActSecEnd][dir1] = (hy[dirActSecEnd - 1])*rp.field[dir2][dirActSecEnd][dir1] + f_Assign.field[dir2][dirActSecEnd - 1][dir1]; }
									else if (direction == directionZ) { f_Assign.field[dir2][dir1][dirActSecEnd] = (hz[dirActSecEnd - 1])*rp.field[dir2][dir1][dirActSecEnd] + f_Assign.field[dir2][dir1][dirActSecEnd - 1]; }
									else { throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown direction."; }
								}
								else { throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown border point type of end sector."; }
								//Sector end (End)

								dirActSecBeg = -1; //sector has been calculated
								continue;
							} //if it is fake direction
							  //Only borders recalculation and continue (End)





							  //Need for fictive points (Begin)
							int fictiveLeft = 0;
							if ((maskActSecBeg != TDefinedBorderPoint) && (maskActSecBeg != TEquationBorderPoint)) { fictiveLeft = 1; }

							int fictiveRight = 0;
							if ((maskActSecEnd != TDefinedBorderPoint) && (maskActSecEnd != TEquationBorderPoint)) { fictiveRight = 1; }

							int dimActSec = dirActSecEnd - dirActSecBeg + 1 + fictiveLeft + fictiveRight;
							//Need for fictive points (End)


							//Preparing scalar objects (Begin)
							int scalDim = dimActSec - 1;
							long double* scalVar = new long double[scalDim + 1];
							long double* scalA = new long double[scalDim + 1];
							long double* scalB = new long double[scalDim + 1];
							long double* scalC = new long double[scalDim + 1];
							long double* scalRP = new long double[scalDim + 1];

							scalA[0] = 0;  scalA[scalDim] = 0;
							scalB[0] = 0;  scalB[scalDim] = 0;
							scalC[0] = 0;  scalC[scalDim] = 0;
							scalRP[0] = 0;  scalRP[scalDim] = 0;


							//Borders points (Begin)
							if ((maskActSecBeg == TDefinedBorderPoint) || (maskActSecBeg == TEquationBorderPoint))
							{
								if (direction == directionX) { scalVar[0] = rp.field[dirActSecBeg][dir2][dir1]; }
								else if (direction == directionY) { scalVar[0] = rp.field[dir2][dirActSecBeg][dir1]; }
								else if (direction == directionZ) { scalVar[0] = rp.field[dir2][dir1][dirActSecBeg]; }
								else { throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown direction."; }
							}
							else if (maskActSecBeg == TPreDefinedBorderPoint)
							{
								scalVar[0] = 0;
								scalA[1] = 1;
								scalB[1] = 0.5;
								scalC[1] = 0.5;
								if (direction == directionX) { scalRP[1] = rp.field[dirActSecBeg][dir2][dir1]; }
								else if (direction == directionY) { scalRP[1] = rp.field[dir2][dirActSecBeg][dir1]; }
								else if (direction == directionZ) { scalRP[1] = rp.field[dir2][dir1][dirActSecBeg]; }
								else { throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown direction."; }
							}
							else if ((maskActSecBeg == TNormalBorderPoint) || (maskActSecBeg == TPreNormalBorderPoint))
							{
								//! need to check acceptable normal direction					   
								scalVar[0] = 0;
								scalA[1] = 1;
								if (direction == directionX)
								{
									scalB[1] = -1 / hx[dirActSecBeg];
									scalC[1] = 1 / hx[dirActSecBeg];
									scalRP[1] = rp.field[dirActSecBeg][dir2][dir1];
								}
								else if (direction == directionY)
								{
									scalB[1] = -1 / hy[dirActSecBeg];
									scalC[1] = 1 / hy[dirActSecBeg];
									scalRP[1] = rp.field[dir2][dirActSecBeg][dir1];
								}
								else if (direction == directionZ)
								{
									scalB[1] = -1 / hz[dirActSecBeg];
									scalC[1] = 1 / hz[dirActSecBeg];
									scalRP[1] = rp.field[dir2][dir1][dirActSecBeg];
								}
								else { throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown direction."; }
							}
							else { throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown border point type of begin sector."; }



							if ((maskActSecEnd == TDefinedBorderPoint) || (maskActSecEnd == TEquationBorderPoint))
							{
								if (direction == directionX) { scalVar[scalDim] = rp.field[dirActSecEnd][dir2][dir1]; }
								else if (direction == directionY) { scalVar[scalDim] = rp.field[dir2][dirActSecEnd][dir1]; }
								else if (direction == directionZ) { scalVar[scalDim] = rp.field[dir2][dir1][dirActSecEnd]; }
								else { throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown direction."; }
							}
							else if (maskActSecEnd == TPreDefinedBorderPoint)
							{
								scalVar[scalDim] = 0;
								scalC[scalDim - 1] = 1;
								scalA[scalDim - 1] = 0.5;
								scalB[scalDim - 1] = 0.5;
								if (direction == directionX) { scalRP[scalDim - 1] = rp.field[dirActSecEnd][dir2][dir1]; }
								else if (direction == directionY) { scalRP[scalDim - 1] = rp.field[dir2][dirActSecEnd][dir1]; }
								else if (direction == directionZ) { scalRP[scalDim - 1] = rp.field[dir2][dir1][dirActSecEnd]; }
								else { throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown direction."; }
							}
							else if ((maskActSecEnd == TNormalBorderPoint) || (maskActSecEnd == TPreNormalBorderPoint))
							{
								scalVar[scalDim] = 0;
								scalC[scalDim - 1] = 1;
								if (direction == directionX)
								{
									scalA[scalDim - 1] = -1 / hx[dirActSecEnd - 1];
									scalB[scalDim - 1] = 1 / hx[dirActSecEnd - 1];
									scalRP[scalDim - 1] = rp.field[dirActSecEnd][dir2][dir1];
								}
								else if (direction == directionY)
								{
									scalA[scalDim - 1] = -1 / hy[dirActSecEnd - 1];
									scalB[scalDim - 1] = 1 / hy[dirActSecEnd - 1];
									scalRP[scalDim - 1] = rp.field[dir2][dirActSecEnd][dir1];
								}
								else if (direction == directionZ)
								{
									scalA[scalDim - 1] = -1 / hz[dirActSecEnd - 1];
									scalB[scalDim - 1] = 1 / hz[dirActSecEnd - 1];
									scalRP[scalDim - 1] = rp.field[dir2][dir1][dirActSecEnd];
								}
								else { throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown direction."; }
							}
							else { throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown border point type of end sector."; }
							//Borders points (End)


							//Actual points (Begin)
							for (int scalLoop = 1 + fictiveLeft;scalLoop < scalDim - fictiveRight;scalLoop++)
							{
								int dirAct = dirActSecBeg + scalLoop - fictiveLeft;

								int i; int j; int k;
								cField* vecTrPtr;
								if (direction == directionX) { i = dirAct; j = dir2; k = dir1; vecTrPtr = &u; }
								else if (direction == directionY) { i = dir2; j = dirAct; k = dir1; vecTrPtr = &v; }
								else if (direction == directionZ) { i = dir2; j = dir1; k = dirAct; vecTrPtr = &w; }
								else { throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown direction."; }
								cField& vecTr = *vecTrPtr;

								if (mask->param[i][j][k].type != TActualPoint) { throw "Error in THydrodynamicsProblem::SolveTransferEquation: Incorrect sector innerr point."; }

								long double h_1X; long double hX; long double hx2X; long double hxhxhX;
								h_1X = hx[i - 1];
								hX = hx[i];
								hx2X = hx[i - 1] + hx[i];
								hxhxhX = hx[i - 1] * hx[i] * (hx[i - 1] + hx[i]) / 2;

								long double h_1Y; long double hY; long double hx2Y; long double hxhxhY;
								h_1Y = hy[j - 1];
								hY = hy[j];
								hx2Y = hy[j - 1] + hy[j];
								hxhxhY = hy[j - 1] * hy[j] * (hy[j - 1] + hy[j]) / 2;

								long double h_1Z; long double hZ; long double hx2Z; long double hxhxhZ;
								h_1Z = hz[k - 1];
								hZ = hz[k];
								hx2Z = hz[k - 1] + hz[k];
								hxhxhZ = hz[k - 1] * hz[k] * (hz[k - 1] + hz[k]) / 2;

								long double h_1; long double h; long double hx2; long double hxhxh;
								if (direction == directionX) { h_1 = h_1X; h = hX; hx2 = hx2X; hxhxh = hxhxhX; }
								else if (direction == directionY) { h_1 = h_1Y; h = hY; hx2 = hx2Y; hxhxh = hxhxhY; }
								else if (direction == directionZ) { h_1 = h_1Z; h = hZ; hx2 = hx2Z; hxhxh = hxhxhZ; }
								else { throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown direction."; }


								long double diffCoef = diffCoefVar.field[i][j][k];

								if (schemeType == 0)
								{
									scalA[scalLoop] = flagConvective*(0.5*(-vecTr.field[i][j][k] / hx2)) + flagDiffusion*(0.5*(-(diffCoef)*h / hxhxh));
									scalB[scalLoop] = 1 / tau + flagDiffusion*(0.5*(diffCoef)*hx2 / hxhxh);
									scalC[scalLoop] = flagConvective*(0.5*(vecTr.field[i][j][k] / hx2)) + flagDiffusion*(0.5*(-(diffCoef)*h_1 / hxhxh));
								}
								else if (schemeType == 1)
								{
									long double vecFlow = vecTr.field[i][j][k];
									if (vecFlow < 0)
									{
										scalA[scalLoop] = flagDiffusion*(0.5*(-(diffCoef)*h / hxhxh));
										scalB[scalLoop] = 1 / tau + flagConvective*(0.5*(-vecTr.field[i][j][k] / h)) + flagDiffusion*(0.5*(diffCoef)*hx2 / hxhxh);
										scalC[scalLoop] = flagConvective*(0.5*(vecTr.field[i][j][k] / h)) + flagDiffusion*(0.5*(-(diffCoef)*h_1 / hxhxh));
									}
									else
									{
										scalA[scalLoop] = flagConvective*(0.5*(-vecTr.field[i][j][k] / h)) + flagDiffusion*(0.5*(-(diffCoef)*h / hxhxh));
										scalB[scalLoop] = 1 / tau + flagConvective*(0.5*(vecTr.field[i][j][k] / h)) + flagDiffusion*(0.5*(diffCoef)*hx2 / hxhxh);
										scalC[scalLoop] = flagDiffusion*(0.5*(-(diffCoef)*h_1 / hxhxh));
									}
								}
								else { throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown scheme type."; }

								long double F = (1 / tau)*f_1.field[i][j][k] + rp.field[i][j][k];

								//Help variables for simple writing F (Begin)
								long double FC;
								long double FD;
								cField* fDir = NULL;
								//Help variables for simple writing F (End)

								if (direction == directionX)
								{
									if (schemeType == 0) { F = F + flagConvectiveX*(0.5*(-u.field[i][j][k])*(f_1.field[i + 1][j][k] - f_1.field[i - 1][j][k]) / hx2X); }
									else if (schemeType == 1)
									{
										long double vecFlow = u.field[i][j][k];
										if (vecFlow < 0) { F = F + flagConvectiveX*(0.5*(-u.field[i][j][k])*(f_1.field[i + 1][j][k] - f_1.field[i][j][k]) / hX); }
										else { F = F + flagConvectiveX*(0.5*(-u.field[i][j][k])*(f_1.field[i][j][k] - f_1.field[i - 1][j][k]) / hX); }
									}
									else { throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown scheme type."; }

									F = F + flagDiffusionX*(0.5*(diffCoef)*((hX*f_1.field[i - 1][j][k] - hx2X*f_1.field[i][j][k] + h_1X*f_1.field[i + 1][j][k]) / hxhxhX));


									//Y (Begin)						   
									if (directionY > direction) { fDir = &f_1; }
									else if (directionY == 1) { fDir = &f_1_3; }
									else if (directionY == 2) { fDir = &f_2_3; }

									FC = 0; FD = 0;
									for (int item = 1;item <= 2;item++)
									{
										cField* f_Ptr = fDir; if (item == 2) { f_Ptr = &f_1; }
										cField& f_ = *f_Ptr;
										if (schemeType == 0) { FC = FC + 0.5*(-v.field[i][j][k])*(f_.field[i][j + 1][k] - f_.field[i][j - 1][k]) / hx2Y; }
										else if (schemeType == 1)
										{
											long double vecFlow = v.field[i][j][k];
											if (vecFlow < 0) { FC = FC + 0.5*(-v.field[i][j][k])*(f_.field[i][j + 1][k] - f_.field[i][j][k]) / hY; }
											else { FC = FC + 0.5*(-v.field[i][j][k])*(f_.field[i][j][k] - f_.field[i][j - 1][k]) / hY; }
										}
										else { throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown scheme type."; }

										FD = FD + 0.5*(diffCoef)*((hY*f.field[i][j - 1][k] - hx2Y*f_.field[i][j][k] + h_1Y*f_.field[i][j + 1][k]) / hxhxhY);
									}
									F = F + flagConvectiveY*FC + flagDiffusionY*FD;
									//Y (End)

									//Z (Begin)						   
									if (directionZ > direction) { fDir = &f_1; }
									else if (directionZ == 1) { fDir = &f_1_3; }
									else if (directionZ == 2) { fDir = &f_2_3; }

									FC = 0; FD = 0;
									for (int item = 1;item <= 2;item++)
									{
										cField* f_Ptr = fDir; if (item == 2) { f_Ptr = &f_1; }
										cField& f_ = *f_Ptr;
										if (schemeType == 0) { FC = FC + 0.5*(-w.field[i][j][k])*(f_.field[i][j][k + 1] - f_.field[i][j][k - 1]) / hx2Z; }
										else if (schemeType == 1)
										{
											long double vecFlow = w.field[i][j][k];
											if (vecFlow < 0) { FC = FC + 0.5*(-w.field[i][j][k])*(f_.field[i][j][k + 1] - f_.field[i][j][k]) / hZ; }
											else { FC = FC + 0.5*(-w.field[i][j][k])*(f_.field[i][j][k] - f_.field[i][j][k - 1]) / hZ; }
										}
										else { throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown scheme type."; }

										FD = FD + 0.5*(diffCoef)*((hZ*f_.field[i][j][k - 1] - hx2Z*f_.field[i][j][k] + h_1Z*f_.field[i][j][k + 1]) / hxhxhZ);
									}
									F = F + flagConvectiveZ*FC + flagDiffusionZ*FD;
									//Z (End)						   
								}
								else if (direction == directionY)
								{
									if (schemeType == 0) { F = F + flagConvectiveY*(0.5*(-v.field[i][j][k])*(f_1.field[i][j + 1][k] - f_1.field[i][j - 1][k]) / hx2Y); }
									else if (schemeType == 1)
									{
										long double vecFlow = v.field[i][j][k];
										if (vecFlow < 0) { F = F + flagConvectiveY*(0.5*(-v.field[i][j][k])*(f_1.field[i][j + 1][k] - f_1.field[i][j][k]) / hY); }
										else { F = F + flagConvectiveY*(0.5*(-v.field[i][j][k])*(f_1.field[i][j][k] - f_1.field[i][j - 1][k]) / hY); }
									}
									else { throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown scheme type."; }


									F = F + flagDiffusionY*(0.5*(diffCoef)*((hY*f_1.field[i][j - 1][k] - hx2Y*f_1.field[i][j][k] + h_1Y*f_1.field[i][j + 1][k]) / hxhxhY));

									//X (Begin)						   
									if (directionX > direction) { fDir = &f_1; }
									else if (directionX == 1) { fDir = &f_1_3; }
									else if (directionX == 2) { fDir = &f_2_3; }

									FC = 0; FD = 0;
									for (int item = 1;item <= 2;item++)
									{
										cField* f_Ptr = fDir; if (item == 2) { f_Ptr = &f_1; }
										cField& f_ = *f_Ptr;
										if (schemeType == 0) { FC = FC + 0.5*(-u.field[i][j][k])*(f_.field[i + 1][j][k] - f_.field[i - 1][j][k]) / hx2X; }
										else if (schemeType == 1)
										{
											long double vecFlow = u.field[i][j][k];
											if (vecFlow < 0) { FC = FC + 0.5*(-u.field[i][j][k])*(f_.field[i + 1][j][k] - f_.field[i][j][k]) / hX; }
											else { FC = FC + 0.5*(-u.field[i][j][k])*(f_.field[i][j][k] - f_.field[i - 1][j][k]) / hX; }
										}
										else { throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown scheme type."; }

										FD = FD + 0.5*(diffCoef)*((hX*f_.field[i - 1][j][k] - hx2X*f_.field[i][j][k] + h_1X*f_.field[i + 1][j][k]) / hxhxhX);
									}
									F = F + flagConvectiveX*FC + flagDiffusionX*FD;
									//X (End)


									//Z (Begin)						   
									if (directionZ > direction) { fDir = &f_1; }
									else if (directionZ == 1) { fDir = &f_1_3; }
									else if (directionZ == 2) { fDir = &f_2_3; }

									FC = 0; FD = 0;
									for (int item = 1;item <= 2;item++)
									{
										cField* f_Ptr = fDir; if (item == 2) { f_Ptr = &f_1; }
										cField& f_ = *f_Ptr;
										if (schemeType == 0) { FC = FC + 0.5*(-w.field[i][j][k])*(f_.field[i][j][k + 1] - f_.field[i][j][k - 1]) / hx2Z; }
										else if (schemeType == 1)
										{
											long double vecFlow = w.field[i][j][k];
											if (vecFlow < 0) { FC = FC + 0.5*(-w.field[i][j][k])*(f_.field[i][j][k + 1] - f_.field[i][j][k]) / hZ; }
											else { FC = FC + 0.5*(-w.field[i][j][k])*(f_.field[i][j][k] - f_.field[i][j][k - 1]) / hZ; }
										}
										else { throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown scheme type."; }

										FD = FD + 0.5*(diffCoef)*((hZ*f_.field[i][j][k - 1] - hx2Z*f_.field[i][j][k] + h_1Z*f_.field[i][j][k + 1]) / hxhxhZ);
									}
									F = F + flagConvectiveZ*FC + flagDiffusionZ*FD;
									//Z (End)
								}
								else if (direction == directionZ)
								{
									if (schemeType == 0) { F = F + flagConvectiveZ*(0.5*(-w.field[i][j][k])*(f_1.field[i][j][k + 1] - f_1.field[i][j][k - 1]) / hx2Z); }
									else if (schemeType == 1)
									{
										long double vecFlow = w.field[i][j][k];
										if (vecFlow < 0) { F = F + flagConvectiveZ*(0.5*(-w.field[i][j][k])*(f_1.field[i][j][k + 1] - f_1.field[i][j][k]) / hZ); }
										else { F = F + flagConvectiveZ*(0.5*(-w.field[i][j][k])*(f_1.field[i][j][k] - f_1.field[i][j][k - 1]) / hZ); }
									}
									else { throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown scheme type."; }

									F = F + flagDiffusionZ*(0.5*(diffCoef)*((hZ*f_1.field[i][j][k - 1] - hx2Z*f_1.field[i][j][k] + h_1Z*f_1.field[i][j][k + 1]) / hxhxhZ));

									//X (Begin)						   
									if (directionX > direction) { fDir = &f_1; }
									else if (directionX == 1) { fDir = &f_1_3; }
									else if (directionX == 2) { fDir = &f_2_3; }

									FC = 0; FD = 0;
									for (int item = 1;item <= 2;item++)
									{
										cField* f_Ptr = fDir; if (item == 2) { f_Ptr = &f_1; }
										cField& f_ = *f_Ptr;
										if (schemeType == 0) { FC = FC + 0.5*(-u.field[i][j][k])*(f_.field[i + 1][j][k] - f_.field[i - 1][j][k]) / hx2X; }
										else if (schemeType == 1)
										{
											long double vecFlow = u.field[i][j][k];
											if (vecFlow < 0) { FC = FC + 0.5*(-u.field[i][j][k])*(f_.field[i + 1][j][k] - f_.field[i][j][k]) / hX; }
											else { FC = FC + 0.5*(-u.field[i][j][k])*(f_.field[i][j][k] - f_.field[i - 1][j][k]) / hX; }
										}
										else { throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown scheme type."; }

										FD = FD + 0.5*(diffCoef)*((hX*f_.field[i - 1][j][k] - hx2X*f_.field[i][j][k] + h_1X*f_.field[i + 1][j][k]) / hxhxhX);
									}
									F = F + flagConvectiveX*FC + flagDiffusionX*FD;
									//X (End)

									//Y (Begin)						   
									if (directionY > direction) { fDir = &f_1; }
									else if (directionY == 1) { fDir = &f_1_3; }
									else if (directionY == 2) { fDir = &f_2_3; }

									FC = 0; FD = 0;
									for (int item = 1;item <= 2;item++)
									{
										cField* f_Ptr = fDir; if (item == 2) { f_Ptr = &f_1; }
										cField& f_ = *f_Ptr;
										if (schemeType == 0) { FC = FC + 0.5*(-v.field[i][j][k])*(f_.field[i][j + 1][k] - f_.field[i][j - 1][k]) / hx2Y; }
										else if (schemeType == 1)
										{
											long double vecFlow = v.field[i][j][k];
											if (vecFlow < 0) { FC = FC + 0.5*(-v.field[i][j][k])*(f_.field[i][j + 1][k] - f_.field[i][j][k]) / hY; }
											else { FC = FC + 0.5*(-v.field[i][j][k])*(f_.field[i][j][k] - f_.field[i][j - 1][k]) / hY; }
										}
										else { throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown scheme type."; }

										FD = FD + 0.5*(diffCoef)*((hY*f_.field[i][j - 1][k] - hx2Y*f_.field[i][j][k] + h_1Y*f_.field[i][j + 1][k]) / hxhxhY);
									}
									F = F + flagConvectiveY*FC + flagDiffusionY*FD;
									//Y (End)
								}
								else { throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown direction."; }
								scalRP[scalLoop] = F;
							} // by act direction
							  //Actual points (End)

							  //Preparing scalar objects (End)


							  //Scalar running (Begin)
							SolveByScalarRunning(scalDim, scalVar, scalA, scalB, scalC, scalRP);
							//Scalar running (End)


							//Assigninig (Begin)
							cField* f_AssignPtr;
							if (direction == 1) { f_AssignPtr = &f_1_3; }
							else if (direction == 2) { f_AssignPtr = &f_2_3; }
							else if (direction == 3) { f_AssignPtr = &f; }
							else { throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown direction."; }
							cField& f_Assign = *f_AssignPtr;

							for (int scalLoop = 0 + fictiveLeft;scalLoop < scalDim + 1 - fictiveRight;scalLoop++)
							{
								int dirAct = dirActSecBeg + scalLoop - fictiveLeft;
								int i; int j; int k;
								if (direction == directionX) { i = dirAct; j = dir2; k = dir1; }
								else if (direction == directionY) { i = dir2; j = dirAct; k = dir1; }
								else if (direction == directionZ) { i = dir2; j = dir1; k = dirAct; }
								else { throw "Error in THydrodynamicsProblem::SolveTransferEquation: Unknown direction."; }

								f_Assign.field[i][j][k] = scalVar[scalLoop];
							}
							//Assigninig (End)

							//Deleting scalar objects (Begin)
							delete[] scalVar;
							delete[] scalA;
							delete[] scalB;
							delete[] scalC;
							delete[] scalRP;
							//Deleting scalar objects (End)

							dirActSecBeg = -1; //sector has been calculated
						} // sector calculating
					} // by sectors of act direction
				} // by direction 2
			} // by direction 1
		} // by fake direction (x,y,z) - for border conditions recalculating
	} // by global direction (x,y,z)
	delete &f_1_3;
	delete &f_2_3;
}

void cGeometry:: prepareMesh()
{
	
		
	int N = mesh->Nx;
	int L = mesh->Ny;
	int M = mesh->Nz;
		for (int i = 0;i<N + 1;i++)
		{
			for (int j = 0;j<L + 1;j++)
			{
				for (int k = 0;k<M + 1;k++)
				{
					int mijk = mask->param[i][j][k].type;
					if ((mijk == TActualPoint) || (mijk == TFictivePoint)) { continue; }
					else if (mijk == TDefinedBorderPoint)
					{
						mask->param[i][j][k].type = TFictivePoint;
					}
					else if (mijk == TEquationBorderPoint)
					{
						mask->param[i][j][k].type = TFictivePoint;
					}
					else if ((mijk == TPreDefinedBorderPoint) || (mijk == TNormalBorderPoint) || (mijk == TPreNormalBorderPoint) || (mijk == TPreEquationBorderPoint))
					{
						mask->param[i][j][k].type = TActualPoint;
					}
					else
					{
						throw "Error in TVPProblem::PrepareVariables: Unknown point mask type.";
					}
				}
			}
		}

		
}