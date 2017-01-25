#include "stdafx.h"
#include "nodetype.h"
#include "service.h"



cArray::cArray(int _n, int _c)
{
	N = _n;
	coord = _c;
	A = new btype [N];
}

cArray::~cArray()
{
	delete [] A;
}

void cMask::init(int _n, int _m, int _l)
{
	Nx= _n;
	Ny = _m;
	Nz = _l;
	param = new nodeMask**[_n];
	for (int i = 0; i < _n; ++i)
	{
		param[i] = new nodeMask*[_m];
		for (int j = 0; j < _m; ++j)
		{
			param[i][j] = new nodeMask[_l];
		}
	}
}

cMask::cMask(const cMask& mask)
{
	Nx = mask.Nx;
	Ny = mask.Ny;
	Nz = mask.Nz;
	param = new nodeMask**[Nx];
	for (int i = 0; i < Nx; ++i)
	{
		param[i] = new nodeMask*[Ny];
		for (int j = 0; j < Ny; ++j)
		{
			param[i][j] = new nodeMask[Nz];
		}
	}

	for (int i = 0; i < Nx; i++)
	   for (int j = 0; j < Ny; j++)
	   for (int k = 0; k < Nz; k++)
		   param[i][j][k] = mask.param[i][j][k];

}

cField::cField(int _n, int _m, int _l)
{
	Nx = _n;
	Ny = _m;
	Nz = _l;
	//coord = _c;
	mask = NULL;
	field = new btype**[Nx];
	for (int i = 0; i<Nx; ++i)
	{
		field[i] = new btype*[Ny];
		for (int j = 0; j<Ny; ++j)
		{
			field[i][j] = new btype[Nz];
		}
	}

	for (int i = 0; i<Nx; ++i)
		for (int j = 0; j<Ny; ++j)
		for (int k = 0; k < Nz; ++k)
			field[i][j][k] = 0.0;
}

void cField::init(int _n, int _m, int _l, int _c,cMask *_mask)
{
	Nx = _n;
	Ny = _m;
	Nz = _l;
	coord = _c;
	mask = _mask;
	field = new btype**[Nx];
	for (int i = 0; i<Nx; ++i)
	{
		field[i] = new btype*[Ny];
		for (int j = 0; j<Ny; ++j)
		{
			field[i][j] = new btype [Nz];
		}
	}
}
//cField(const cField&);
cField::~cField()
{
	
	for (int i = 0; i<Nx; ++i)
	   for (int j = 0; j<Ny; ++j)
		   delete[] field[i][j];

	   for (int i = 0; i<Nx; ++i)
		   delete[] field[i];
	   delete[] field;
}

cField& cField::operator = (const btype& _c)
{
	for (int i = 0; i < Nx;i++)
	for (int j = 0; j < Ny; j++)
	for (int k = 0; k < Nz; k++)
		field[i][j][k] = _c;
	return *this;
}

cField& cField::operator = (const cField& _c)
{
	for (int i = 0; i < Nx; i++)
	for (int j = 0; j < Ny; j++)
	for (int k = 0; k < Nz; k++)
		field[i][j][k] = _c.field[i][j][k];
	return *this;
}

cMask::~cMask()
{

	for (int i = 0; i<Nx; ++i)
	for (int j = 0; j<Ny; ++j)
		delete[] param[i][j];

	for (int i = 0; i<Nx; ++i)
		delete[] param[i];
	delete[] param;
}
/*
template <typename T> cArray<T>::cArray(const cArray<T>& B)
{

}

template <typename T> cArray<T> cArray<T>::operator + (cArray<T>& B)
{
}

template <typename T> cArray<T> cArray<T>::operator - (cArray<T>& B)
{
}

template <typename T> cArray<T>& cArray<T>::operator = (const cArray<T>& B)
{
}

*/
