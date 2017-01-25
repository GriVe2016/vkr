
#include "nodetype.h"
#include "param.h"
#include "service.h"


struct sCoord
{
	cField *currnt, *pred, *interm, *f;
	int Nx, Ny, Nz;
};
class cGeometry
{ 

	void initMemoryXYZ();
	void CreateGeometryNodes();
	void CreateGeometryGrid();
	void CreateMultyNodes();
	void CreateMultyGrids();
	
    public:

		sCoord sX, sY, sZ;

		int Nx, Ny, Nz;
		int NxU, NyU, NzU;
		int NxV, NyV, NzV;
		int NxW, NyW, NzW;
		int NxP, NyP, NzP;

		btype h, lengthX, lengthY, lengthZ;

		btype *X, *Y, *Z;
		btype *XU, *YU, *ZU;
		btype *XV, *YV, *ZV;
		btype *XW, *YW, *ZW;
		btype *XP, *YP, *ZP;

		cField U, U1, Us, Uprev, Uold;
		cField V, V1, Vs, Vprev, Vold;
		cField W, W1, Ws, Wprev, Wold;
		cField P, P1, Ps, Pprev, Pold;
		cField Fx, Fy, Fz, Fp;
		cField C, C1;
		cField Rho;
		cField Mu;

		cMask Mask;
		cMask MaskU, MaskV, MaskW, MaskP, MaskRho,MaskC;


        void setFToZero();
		//void initBy(int, int, int, btype ***, btype ***);
		cGeometry(sDomainParameters meshParam);
	    ~cGeometry();



};

