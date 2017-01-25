#include "stdafx.h"
#include "nodetype.h"

#ifndef NORMAL_VECTOR
#define NORMAL_VECTOR

#define _x 1
#define _y 2
#define _z 3
#define _p 4

struct NormalVector
{
	btype x, y, z;
	NormalVector() 
	{
		x = 0.;
		y = 0.;
		z = 0.;
	}
	NormalVector(btype a, btype b, btype c) : x(a), y(b), z(c) {}
	void set(btype a, btype b, btype c)
	{
		x = a;
		y = b;
		z = c;
	}
};

const int TFictivePoint = 0;
const int TActualPoint = 1;
const int TBorderPoint = 2;
const int TNormalBorderPoint = 3;
const int TPreDefinedBorderPoint = 4;
const int TDefinedBorderPoint = 5;
const int TPreNormalBorderPoint = 6;

struct nodeMask
{
	int type;
	NormalVector normal;
};


template <typename T> void allocateArray(T *&arr, int n)
{
	arr = new T[n];
}

template <typename T> void allocateCube(T ***&arr, int n, int m, int l)
{
	arr = new T**[n];
	for (int i = 0; i<n; ++i)
	{
		arr[i] = new T*[m];
		for (int j = 0; j<m; ++j)
		{
			arr[i][j] = new T[l];
		}
	}
}



class cArray
{
	btype *A;
	int N;
	int coord;
public:
	cArray(int,int);
	//cArray(const cArray&);
	~cArray();

	/*
	btype operator [] (int i) const
	{
		return A[i];
	}*/
	//cArray<T> operator + (cArray<T>&);
	//cArray<T> operator - (cArray<T>&);
	//cArray& operator = (const cArray&);


};
class cMask
{
public:	//btype ***F;
	nodeMask ***param;
	int Nx, Ny, Nz;
	//int coord;


	cMask()
	{
		param = NULL;
		Nx = Ny = Nz = 0;
	}
	void init(int,int,int);
	cMask(const cMask&);
	~cMask();

	/*
	btype operator [] (int i) const
	{
	return A[i];
	}
	*/
	//cArray<T> operator + (cArray<T>&);
	//cArray<T> operator - (cArray<T>&);
	//cArray& operator = (const cArray&);


};

class cField
{
   public:
	btype ***field;
	//nodeMask ***mask;
	int Nx,Ny,Nz;
	int coord;
	cMask *mask;

	cField()
	{
		field = NULL;
		Nx = Ny = Nz = coord = 0;
	};

	cField(int _n, int _m, int _l);

	void init(int _n, int _m, int _l,int,cMask *_mask);
	//cField(const cField&);
	~cField();

	cField& operator = (const btype&);
	cField& operator = (const cField&);

	/*
	btype operator [] (int i) const
	{
		return A[i];
	}
	*/
	//cArray<T> operator + (cArray<T>&);
	//cArray<T> operator - (cArray<T>&);
	//cArray& operator = (const cArray&);


};



/*
template <typename T> class cArray
{
	T *A;
	int n;
public:
	cArray(int );
	cArray(const cArray<T>&);
	~cArray();


	T operator [] (int i) const
	{
		return A[i];
	}
	cArray<T> operator + (cArray<T>&);
	cArray<T> operator - (cArray<T>&);
	cArray<T>& operator = (const cArray<T>&);


};
*/

#endif