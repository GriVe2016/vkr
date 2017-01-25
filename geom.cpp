
#include "stdafx.h"
#include "nodetype.h"
#include "geom.h"
#include "service.h"
#include <math.h>

cGeometry::cGeometry(sDomainParameters meshParam)
{

	btype iPart;
	int endIndex;

	lengthX = meshParam.lengthX;
	lengthY = meshParam.lengthY;
	lengthZ = meshParam.lengthZ;
	iPart = 0;
	modf(lengthX / meshParam.h, &iPart);
	endIndex = (int)iPart;
	Nx = endIndex + 1;
		
	iPart = 0;
	modf(lengthY / meshParam.h, &iPart);
	endIndex = (int)iPart;
	Ny = endIndex + 1;

	iPart = 0;
	modf(lengthZ / meshParam.h, &iPart);
	endIndex = (int)iPart;
	Nz = endIndex + 1;

	// Размеры разнесенных сеток
	NxU = NxV = NxW = NxP = Nx + 1;
	NyU = NyV = NyW = NyP = Ny + 1;
	NzU = NzV = NzW = NzP = Nz + 1;
	--NxU;
	--NyV;
	--NzW;

	initMemoryXYZ();
	CreateGeometryNodes();
	CreateGeometryGrid();
	CreateMultyNodes();
	CreateMultyGrids();



	U.init (NxU, NyU, NzU, _x, &MaskU);
	U1.init(NxU, NyU, NzU, _x, &MaskU);
	Us.init(NxU, NyU, NzU, _x, &MaskU);
	Fx.init(NxU, NyU, NzU, _x, &MaskU);

	V.init (NxV, NyV, NzV, _y, &MaskV);
	V1.init(NxV, NyV, NzV, _y, &MaskV);
	Vs.init(NxV, NyV, NzV, _y, &MaskV);
	Fy.init(NxV, NyV, NzV, _y, &MaskV);

	W.init (NxW, NyW, NzW, _z, &MaskW);
	W1.init(NxW, NyW, NzW, _z, &MaskW);
	Ws.init(NxW, NyW, NzW, _z, &MaskW);
	Fz.init(NxW, NyW, NzW, _z, &MaskW);

	P.init(NxP, NyP, NzP, _p, &MaskP);
	P1.init(NxP, NyP, NzP, _p, &MaskP);
	//Ps.init(NxP, NyP, NzP, _p, &MaskP);
	C.init(NxP, NyP, NzP, _p, &MaskP);
	C1.init(NxP, NyP, NzP, _p, &MaskP);
	Mu.init(NxP, NyP, NzP, _p, &MaskP);
	Rho.init(NxP, NyP, NzP, _p, &MaskP);
	//initMemoryUVWP();

}

cGeometry::~cGeometry()
{

}
/*
void cGeometry::setFToZero()
{
	for (int i = 0; i < Nx; i++)
	for (int j = 0; j < Ny + 1; j++)
	for (int k = 0; k < Nz + 1; k++)
		Fx[i][j][k] = 0.00;

	for (int i = 0; i < Nx + 1; i++)
	for (int j = 0; j < Ny; j++)
	for (int k = 0; k < Nz + 1; k++)
		Fy[i][j][k] = 0.00;

	for (int i = 0; i < Nx + 1; i++)
	for (int j = 0; j < Ny + 1; j++)
	for (int k = 0; k < Nz; k++)
		Fz[i][j][k] = 0.00;
}
*/

void cGeometry::initMemoryXYZ()
{
	// Узлы сетки и шаги
	allocateArray(X, Nx);
	allocateArray(Y, Ny);
	allocateArray(Z, Nz);
	/*
	allocateArray(Hx, Nx);
	allocateArray(Hy, Ny);
	allocateArray(Hz, Nz);
	*/
	// Маска
	Mask.init( Nx, Ny, Nz);

	// смещенные сетки
	allocateArray(XU, NxU);
	allocateArray(YU, NyU);
	allocateArray(ZU, NzU);
	/*
	allocateArray(HxU, NxU);
	allocateArray(HyU, NyU);
	allocateArray(HzU, NzU);
	*/
	allocateArray(XV, NxV);
	allocateArray(YV, NyV);
	allocateArray(ZV, NzV);
	/*
	allocateArray(HxV, NxV);
	allocateArray(HyV, NyV);
	allocateArray(HzV, NzV);
	*/

	allocateArray(XW, NxW);
	allocateArray(YW, NyW);
	allocateArray(ZW, NzW);
	/*
	allocateArray(HxW, NxW);
	allocateArray(HyW, NyW);
	allocateArray(HzW, NzW);
	*/
	allocateArray(XP, NxP);
	allocateArray(YP, NyP);
	allocateArray(ZP, NzP);
	/*
	allocateArray(HxP, NxP);
	allocateArray(HyP, NyP);
	glAllocVector(HzP, NzP);
    */
	MaskU.init(NxU, NyU, NzU);
	MaskV.init(NxV, NyV, NzV);
	MaskW.init(NxW, NyW, NzW);
	MaskP.init(NxP, NyP, NzP);
	MaskRho.init(NxP, NyP, NzP);
	MaskC.init(NxP, NyP, NzP);
	
}

void cGeometry::CreateGeometryNodes()
{
	// Узлы, шаги сетки
	for (int i = 0; i<Nx; ++i)
		X[i] = lengthX*i / (Nx - 1);
	/*
	for (int i = 1; i<Nx; ++i)
		Hx[i] = X[i] - X[i - 1];
	Hx[0] = Hx[1];
    */

	for (int j = 0; j<Ny; ++j)
		Y[j] = lengthY*j / (Ny - 1);
	/*
	for (int j = 1; j<Ny; ++j)
		Hy[j] = Y[j] - Y[j - 1];
	Hy[0] = Hy[1];
    */
	for (int k = 0; k<Nz; ++k)
		Z[k] = lengthZ*k / (Nz - 1);
	/*
	for (int k = 1; k<Nz; ++k)
		Hz[k] = Z[k] - Z[k - 1];
	Hy[0] = Hy[1];
	*/
}

void cGeometry::CreateGeometryGrid()
{
	// Сначала все точки - расчетные
	for (int i = 0; i<Nx; ++i)
	for (int j = 0; j<Ny; ++j)
	for (int k = 0; k<Nz; ++k)
	{
		Mask.param[i][j][k].type = TActualPoint;
		Mask.param[i][j][k].normal.set(0., 0., 0.);
	}

	// Левая и правая границы
	for (int j = 1; j<Ny - 1; ++j)
	for (int k = 1; k<Nz - 1; ++k)
	{
		Mask.param[Nx - 1][j][k].type = Mask.param[0][j][k].type = TBorderPoint;
		Mask.param[0][j][k].normal.set(-1., 0., 0.);
		Mask.param[Nx - 1][j][k].normal.set(1., 0., 0.);
	}

	// Ближняя и дальняя границы
	for (int i = 1; i<Nx - 1; ++i)
	for (int k = 1; k<Nz - 1; ++k)
	{
		Mask.param[i][Ny - 1][k].type = Mask.param[i][0][k].type = TBorderPoint;
		Mask.param[i][0][k].normal.set(0., -1., 0.);
		Mask.param[i][Ny - 1][k].normal.set(0., 1., 0.);
	}

	// Нижняя и верхняя границы
	for (int i = 1; i<Nx - 1; ++i)
	for (int j = 1; j<Ny - 1; ++j)
	{
		Mask.param[i][j][Nz - 1].type = Mask.param[i][j][0].type = TBorderPoint;
		Mask.param[i][j][0].normal.set(0, 0, -1);
		Mask.param[i][j][Nz - 1].normal.set(0, 0, 1);
	}

	// Фиктивные линии по ребрам
	for (int i = 0; i<Nx; ++i)
	{
		Mask.param[i][0][0].type = TFictivePoint;
		Mask.param[i][0][Nz - 1].type = TFictivePoint;
		Mask.param[i][Ny - 1][0].type = TFictivePoint;
		Mask.param[i][Ny - 1][Nz - 1].type = TFictivePoint;
	}

	for (int j = 0; j<Ny; ++j)
	{
		Mask.param[0][j][0].type = TFictivePoint;
		Mask.param[0][j][Nz - 1].type = TFictivePoint;
		Mask.param[Nx - 1][j][0].type = TFictivePoint;
		Mask.param[Nx - 1][j][Nz - 1].type = TFictivePoint;
	}

	for (int k = 0; k<Nz; ++k)
	{
		Mask.param[0][0][k].type = TFictivePoint;
		Mask.param[0][Ny - 1][k].type = TFictivePoint;
		Mask.param[Nx - 1][0][k].type = TFictivePoint;
		Mask.param[Nx - 1][Ny - 1][k].type = TFictivePoint;
	}
}

void cGeometry::CreateMultyNodes()
{
	// Узлы для компоненты U скорости
	for (int i = 0; i<NxU; ++i)
		XU[i] = X[i];

	for (int j = 1; j<NyU - 1; ++j)
		YU[j] = (Y[j - 1] + Y[j]) / 2.0;
	YU[0] = Y[0] - (Y[1] - Y[0]) / 2.0;
	YU[NyU - 1] = Y[Ny - 1] + (Y[Ny - 1] - Y[Ny - 2]) / 2.0;

	for (int k = 1; k<NzU - 1; ++k)
		ZU[k] = (Z[k - 1] + Z[k]) / 2.0;
	ZU[0] = Z[0] - (Z[1] - Z[0]) / 2.0;
	ZU[NzU - 1] = Z[Nz - 1] + (Z[Nz - 1] - Z[Nz - 2]) / 2.0;

	/*
	for (int i = 1; i<NxU; ++i)
		HxU[i] = XU[i] - XU[i - 1];
	HxU[0] = HxU[1];

	for (int j = 1; j<NyU; ++j)
		HyU[j] = YU[j] - YU[j - 1];
	HyU[0] = HyU[1];

	for (int k = 1; k<NzU; ++k)
		HzU[k] = ZU[k] - ZU[k - 1];
	HzU[0] = HzU[1];
	*/

	// Узлы для компоненты V скорости
	for (int i = 1; i<NxV - 1; ++i)
		XV[i] = (X[i - 1] + X[i]) / 2.0;
	XV[0] = X[0] - (X[1] - X[0]) / 2.0;
	XV[NxV - 1] = X[Nx - 1] + (X[Nx - 1] - X[Nx - 2]) / 2.0;

	for (int j = 0; j<NyV; ++j)
		YV[j] = Y[j];

	for (int k = 1; k<NzV - 1; ++k)
		ZV[k] = (Z[k - 1] + Z[k]) / 2.0;
	ZV[0] = Z[0] - (Z[1] - Z[0]) / 2.0;
	ZV[NzV - 1] = Z[Nz - 1] + (Z[Nz - 1] - Z[Nz - 2]) / 2.0;

	/*
	for (int i = 1; i<NxV; ++i)
		HxV[i] = XV[i] - XV[i - 1];
	HxV[0] = HxV[1];

	for (int j = 1; j<NyV; ++j)
		HyV[j] = YV[j] - YV[j - 1];
	HyV[0] = HyV[1];

	for (int k = 1; k<NzV; ++k)
		HzV[k] = ZV[k] - ZV[k - 1];
	HzV[0] = HzV[1];
	*/

	// Узлы для компоненты W скорости
	for (int i = 1; i<NxW - 1; ++i)
		XW[i] = (X[i - 1] + X[i]) / 2.0;
	XW[0] = X[0] - (X[1] - X[0]) / 2.0;
	XW[NxW - 1] = X[Nx - 1] + (X[Nx - 1] - X[Nx - 2]) / 2.0;

	for (int j = 1; j<NyW - 1; ++j)
		YW[j] = (Y[j - 1] + Y[j]) / 2.0;
	YW[0] = Y[0] - (Y[1] - Y[0]) / 2.0;
	YW[NyW - 1] = Y[Ny - 1] + (Y[Ny - 1] - Y[Ny - 2]) / 2.0;

	for (int k = 0; k<NzW; ++k)
		ZW[k] = Z[k];

	/*
	for (int i = 1; i<NxW; ++i)
		HxW[i] = XW[i] - XW[i - 1];
	HxW[0] = HxW[1];

	for (int j = 1; j<NyW; ++j)
		HyW[j] = YW[j] - YW[j - 1];
	HyW[0] = HyW[1];

	for (int k = 1; k<NzW; ++k)
		HzW[k] = ZW[k] - ZW[k - 1];
	HzW[0] = HzW[1];
	*/
	// Узлы для давления
	for (int i = 1; i<NxP - 1; ++i)
		XP[i] = (X[i - 1] + X[i]) / 2.0;
	XP[0] = X[0] - (X[1] - X[0]) / 2.0;
	XP[NxP - 1] = X[Nx - 1] + (X[Nx - 1] - X[Nx - 2]) / 2.0;

	for (int j = 1; j<NyP - 1; ++j)
		YP[j] = (Y[j - 1] + Y[j]) / 2.0;
	YP[0] = Y[0] - (Y[1] - Y[0]) / 2.0;
	YP[NyP - 1] = Y[Ny - 1] + (Y[Ny - 1] - Y[Ny - 2]) / 2.0;

	for (int k = 1; k<NzP - 1; ++k)
		ZP[k] = (Z[k - 1] + Z[k]) / 2.0;
	ZP[0] = Z[0] - (Z[1] - Z[0]) / 2.0;
	ZP[NzP - 1] = Z[Nz - 1] + (Z[Nz - 1] - Z[Nz - 2]) / 2.0;

	/*
	for (int i = 1; i<NxP; ++i)
		HxP[i] = XP[i] - XP[i - 1];
	HxP[0] = HxP[1];

	for (int j = 1; j<NyP; ++j)
		HyP[j] = YP[j] - YP[j - 1];
	HyP[0] = HyP[1];

	for (int k = 1; k<NzP; ++k)
		HzP[k] = ZP[k] - ZP[k - 1];
	HzP[0] = HzP[1];
	*/
}
//*************************************
void cGeometry::CreateMultyGrids()
{
	// Заполнение маски для U //

	for (int i = 0; i<NxU; ++i)
	for (int j = 0; j<NyU; ++j)
	for (int k = 0; k<NzU; ++k)
	{
		MaskU.param[i][j][k].type = TActualPoint;
		MaskU.param[i][j][k].normal.set(0., 0., 0.);
	}

	// Границы области
	for (int j = 1; j<NyU - 1; ++j)
	for (int k = 1; k<NzU - 1; ++k)
	{
		MaskU.param[0][j][k].type = MaskU.param[NxU - 1][j][k].type = TNormalBorderPoint;
		//MaskU[0][j][k].mask = MaskU[NxU-1][j][k].mask = TEquationBorderPoint;
		MaskU.param[0][j][k].normal.set(-1., 0., 0.);
		MaskU.param[NxU - 1][j][k].normal.set(1., 0., 0.);
	}

	for (int i = 1; i<NxU - 1; ++i)
	for (int k = 1; k<NzU - 1; ++k)
	{
		MaskU.param[i][0][k].type = MaskU.param[i][NyU - 1][k].type = TPreDefinedBorderPoint;
		MaskU.param[i][0][k].normal.set(0., -1., 0.);
		MaskU.param[i][NyU - 1][k].normal.set(0., 1., 0.);
	}

	for (int i = 1; i<NxU - 1; ++i)
	for (int j = 1; j<NyU - 1; ++j)
	{
		MaskU.param[i][j][0].type = MaskU.param[i][j][NzU - 1].type = TPreDefinedBorderPoint;
		MaskU.param[i][j][0].normal.set(0., 0., -1.);
		MaskU.param[i][j][NzU - 1].normal.set(0., 0., 1.);
	}

	// Фиктивные линии по ребрам
	for (int i = 0; i<NxU; ++i)
	{
		MaskU.param[i][0][0].type = TFictivePoint;
		MaskU.param[i][0][NzU - 1].type = TFictivePoint;
		MaskU.param[i][NyU - 1][0].type = TFictivePoint;
		MaskU.param[i][NyU - 1][NzU - 1].type = TFictivePoint;
	}

	for (int j = 0; j<NyU; ++j)
	{
		MaskU.param[0][j][0].type = TFictivePoint;
		MaskU.param[0][j][NzU - 1].type = TFictivePoint;
		MaskU.param[NxU - 1][j][0].type = TFictivePoint;
		MaskU.param[NxU - 1][j][NzU - 1].type = TFictivePoint;
	}

	for (int k = 0; k<NzU; ++k)
	{
		MaskU.param[0][0][k].type = TFictivePoint;
		MaskU.param[0][NyU - 1][k].type = TFictivePoint;
		MaskU.param[NxU - 1][0][k].type = TFictivePoint;
		MaskU.param[NxU - 1][NyU - 1][k].type = TFictivePoint;
	}


	// Заполнение маски для V //

	for (int i = 0; i<NxV; ++i)
	for (int j = 0; j<NyV; ++j)
	for (int k = 0; k<NzV; ++k)
	{
		MaskV.param[i][j][k].type = TActualPoint;
		MaskV.param[i][j][k].normal.set(0., 0., 0.);
	}

	// Границы области
	for (int j = 1; j<NyV - 1; ++j)
	for (int k = 1; k<NzV - 1; ++k)
	{
		MaskV.param[0][j][k].type = MaskV.param[NxV - 1][j][k].type = TPreDefinedBorderPoint;
		MaskV.param[0][j][k].normal.set(-1., 0., 0.);
		MaskV.param[NxV - 1][j][k].normal.set(1., 0., 0.);
	}

	for (int i = 1; i<NxV - 1; ++i)
	for (int k = 1; k<NzV - 1; ++k)
	{
		MaskV.param[i][0][k].type = MaskV.param[i][NyV - 1][k].type = TDefinedBorderPoint;
		MaskV.param[i][0][k].normal.set(0., -1., 0.);
		MaskV.param[i][NyV - 1][k].normal.set(0., 1., 0.);
	}

	for (int i = 1; i<NxV - 1; ++i)
	for (int j = 1; j<NyV - 1; ++j)
	{
		MaskV.param[i][j][0].type = MaskV.param[i][j][NzV - 1].type = TPreDefinedBorderPoint;
		MaskV.param[i][j][0].normal.set(0., 0., -1.);
		MaskV.param[i][j][NzV - 1].normal.set(0., 0., 1.);
	}

	// Фиктивные линии по ребрам
	for (int i = 0; i<NxV; ++i)
	{
		MaskV.param[i][0][0].type = TFictivePoint;
		MaskV.param[i][0][NzV - 1].type = TFictivePoint;
		MaskV.param[i][NyV - 1][0].type = TFictivePoint;
		MaskV.param[i][NyV - 1][NzV - 1].type = TFictivePoint;
	}

	for (int j = 0; j<NyV; ++j)
	{
		MaskV.param[0][j][0].type = TFictivePoint;
		MaskV.param[0][j][NzV - 1].type = TFictivePoint;
		MaskV.param[NxV - 1][j][0].type = TFictivePoint;
		MaskV.param[NxV - 1][j][NzV - 1].type = TFictivePoint;
	}

	for (int k = 0; k<NzV; ++k)
	{
		MaskV.param[0][0][k].type = TFictivePoint;
		MaskV.param[0][NyV - 1][k].type = TFictivePoint;
		MaskV.param[NxV - 1][0][k].type = TFictivePoint;
		MaskV.param[NxV - 1][NyV - 1][k].type = TFictivePoint;
	}


	// Заполнение маски для W //

	for (int i = 0; i<NxW; ++i)
	for (int j = 0; j<NyW; ++j)
	for (int k = 0; k<NzW; ++k)
	{
		MaskW.param[i][j][k].type = TActualPoint;
		MaskW.param[i][j][k].normal.set(0., 0., 0.);
	}
	// Границы области
	for (int j = 1; j<NyW - 1; ++j)
	for (int k = 1; k<NzW - 1; ++k)
	{
		MaskW.param[0][j][k].type = MaskW.param[NxW - 1][j][k].type = TPreDefinedBorderPoint;
		MaskW.param[0][j][k].normal.set(-1., 0., 0.);
		MaskW.param[NxW - 1][j][k].normal.set(1., 0., 0.);
	}

	for (int i = 1; i<NxW - 1; ++i)
	for (int k = 1; k<NzW - 1; ++k)
	{
		MaskW.param[i][0][k].type = MaskW.param[i][NyW - 1][k].type = TPreDefinedBorderPoint;
		MaskW.param[i][0][k].normal.set(0., -1., 0.);
		MaskW.param[i][NyW - 1][k].normal.set(0., 1., 0.);
	}

	for (int i = 1; i<NxW - 1; ++i)
	for (int j = 1; j<NyW - 1; ++j)
	{
		MaskW.param[i][j][0].type = MaskW.param[i][j][NzW - 1].type = TDefinedBorderPoint;
		MaskW.param[i][j][0].normal.set(0., 0., -1.);
		MaskW.param[i][j][NzW - 1].normal.set(0., 0., 1.);
	}

	// Фиктивные линии по ребрам
	for (int i = 0; i<NxW; ++i)
	{
		MaskW.param[i][0][0].type = TFictivePoint;
		MaskW.param[i][0][NzW - 1].type = TFictivePoint;
		MaskW.param[i][NyW - 1][0].type = TFictivePoint;
		MaskW.param[i][NyW - 1][NzW - 1].type = TFictivePoint;
	}

	for (int j = 0; j<NyW; ++j)
	{
		MaskW.param[0][j][0].type = TFictivePoint;
		MaskW.param[0][j][NzW - 1].type = TFictivePoint;
		MaskW.param[NxW - 1][j][0].type = TFictivePoint;
		MaskW.param[NxW - 1][j][NzW - 1].type = TFictivePoint;
	}

	for (int k = 0; k<NzW; ++k)
	{
		MaskW.param[0][0][k].type = TFictivePoint;
		MaskW.param[0][NyW - 1][k].type = TFictivePoint;
		MaskW.param[NxW - 1][0][k].type = TFictivePoint;
		MaskW.param[NxW - 1][NyW - 1][k].type = TFictivePoint;
	}


	// Заполнение маски для P
	for (int i = 0; i<NxP; ++i)
	for (int j = 0; j<NyP; ++j)
	for (int k = 0; k<NzP; ++k)
	{
		MaskP.param[i][j][k].type = TActualPoint;
		MaskP.param[i][j][k].normal.set(0., 0., 0.);
	}

	// Границы области
	for (int j = 1; j<NyP - 1; ++j)
	for (int k = 1; k<NzP - 1; ++k)
	{
		MaskP.param[0][j][k].type = MaskP.param[NxP - 1][j][k].type = TPreDefinedBorderPoint;
		MaskP.param[0][j][k].normal.set(-1., 0., 0.);
		MaskP.param[NxP - 1][j][k].normal.set(1., 0., 0.);
	}

	for (int i = 1; i<NxP - 1; ++i)
	for (int k = 1; k<NzP - 1; ++k)
	{
		MaskP.param[i][0][k].type = MaskP.param[i][NyP - 1][k].type = TFictivePoint;
		MaskP.param[i][0][k].normal.set(0., -1., 0.);
		MaskP.param[i][NyP - 1][k].normal.set(0., 1., 0.);
	}

	for (int i = 1; i<NxP - 1; ++i)
	for (int j = 1; j<NyP - 1; ++j)
	{
		MaskP.param[i][j][0].type = MaskP.param[i][j][NzP - 1].type = TFictivePoint;
		MaskP.param[i][j][0].normal.set(0., 0., -1.);
		MaskP.param[i][j][NzP - 1].normal.set(0., 0., 1.);
	}

	// Фиктивные линии по ребрам
	for (int i = 0; i<NxP; ++i)
	{
		MaskP.param[i][0][0].type = TFictivePoint;
		MaskP.param[i][0][NzP - 1].type = TFictivePoint;
		MaskP.param[i][NyP - 1][0].type = TFictivePoint;
		MaskP.param[i][NyP - 1][NzP - 1].type = TFictivePoint;
	}

	for (int j = 0; j<NyP; ++j)
	{
		MaskP.param[0][j][0].type = TFictivePoint;
		MaskP.param[0][j][NzP - 1].type = TFictivePoint;
		MaskP.param[NxP - 1][j][0].type = TFictivePoint;
		MaskP.param[NxP - 1][j][NzP - 1].type = TFictivePoint;
	}

	for (int k = 0; k<NzP; ++k)
	{
		MaskP.param[0][0][k].type = TFictivePoint;
		MaskP.param[0][NyP - 1][k].type = TFictivePoint;
		MaskP.param[NxP - 1][0][k].type = TFictivePoint;
		MaskP.param[NxP - 1][NyP - 1][k].type = TFictivePoint;
	}
}


void cGeometry::setFToZero()
{
	Fx = 0.0;
	Fy = 0.0;
	Fz = 0.0;
	return;
}