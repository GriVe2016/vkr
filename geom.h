#include "param.h"

#ifndef NODETYPE
#define NODETYPE
   typedef double btype;



class NormalVector
{
public:
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
const int TEquationBorderPoint=7;//=?????
const int TPreEquationBorderPoint = 8;

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
	btype bgn, length, h;
	int N;
public:
	cArray() {};
	void init(int _N,btype _bgn,btype _length,btype _h)
        {
            N = _N;
            bgn = _bgn; 
            length = _length;
            h = _h;
        }
        
	

	btype operator [] (int i) const
	{
		return bgn + i*h;
	}
	//cArray<T> operator + (cArray<T>&);
	//cArray<T> operator - (cArray<T>&);
	//cArray& operator = (const cArray&);

	                                                                            
};

class cMesh 
{
public:
    int Nx,Ny,Nz;
    cArray X,Y,Z;
	
 
	cMesh()
	{
		Nx = Ny = Nz = 0;
	};
	void init(int _Nx, int _Ny, int _Nz, btype x_bgn, btype y_bgn, btype z_bgn, btype x_h, 
		btype y_h, btype z_h, btype x_length, btype y_length, btype z_length) 
	{
		Nx = _Nx;
		Ny = _Ny;
		Nz = _Nz;
		X.init(_Nx, x_bgn, x_length, x_h);
		Y.init(_Ny, y_bgn, y_length, y_h);
		Z.init(_Nz, z_bgn, z_length, z_h);
    }
	
};


class cMask
{
public:	//btype ***F;
	nodeMask ***param;
	int Nx, Ny, Nz;
	

	
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
	int Nx, Ny, Nz;
	cMask *mask;
	cMesh *mesh;
	cField() {};
	void init(cMask *_mask, cMesh *_mesh, int _Nx, int _Ny, int _Nz)
	{
		Nx = _Nx;
		Ny = _Ny;
		Nz = _Nz;
		mask = _mask;
		field = new btype**[Nx];
		for (int i = 0; i<Nx; ++i)
		{
			field[i] = new btype*[Ny];
			for (int j = 0; j<Ny; ++j)
			{
				field[i][j] = new btype[Nz];
			}
		}
		
	}

	//~cField();

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



class cGeometry
{
public:	
	cMask *mask, *maskU, *maskV, *maskW, *maskP, *maskC, *maskRho;
	cMesh *mesh, *meshU, *meshV, *meshW, *meshP;
	

	cField *U, *U1, *Us, *Uprev, *Uold;
	cField *V, *V1, *Vs, *Vprev, *Vold;
	cField *W, *W1, *Ws, *Wprev, *Wold;
	cField *P, *P1, *Ps, *Pprev, *Pold;
	cField *Fx, *Fy, *Fz, *Fp;
	cField *C, *C1;
	cField *Rho;
	cField *Mu;



	cGeometry(const sDomainParameters meshParam);
	void CreateGeometryGrid();
	void CreateMultyGrids(int NxU,int NxV, int NxW,int  NxP,int  NyU,int  NyV,int  NyW,int  NyP,int  NzU,int  NzV,int  NzW,int  NzP);
	void setFToZero();
	void SetInitCondVis();
	void SetInitCondDensity();
	void SetInitCondRho(cInputParameters inpParams);
	void SetInitCondC( cInputParameters inpParams);
	void SetInitCondV(cInputParameters inpParams);
	void prepareMesh();
	void SolveTransferEquation__(cField& lVecCur, cField& lVecPrev,
		cField&  lU, cField& lV, cField& lW, cField& lRP,
		const cMesh& lGrid, cField& lDiffusionCoefficient, int lSchemeType);
	void SolveByScalarRunning(int lM, long double* lU, const long double* lA, const long double* lB, const long double* lC, const long double* lF);
};




#endif /* NODETYPE */
