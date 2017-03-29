
#include "param.h"
#include "geom.h"
#include <math.h>

class TNode
{
public:
	long double x, y, z;
	long double xPrev, yPrev, zPrev;
	long double xRef, yRef, zRef;
	long double xForce, yForce, zForce;
	long double xVel, yVel, zVel;
	long double p;
	int type;
	long double stretchingStiffness;

	//TImmersedBoundary *boundary;
	struct Neighbor {
		TNode *next;
		long double initialDistanceNext;
		TNode *prev;
		long double initialDistancePrev;
		//unordered_map<AXIS, long double> initialCurvatures;
	} neighbors;

	// debug
	//map<tuple<int, int, int>, TNode*> surround;

	/*
	long double GetStiffnessX(int timestep) { return 1; };
	long double GetStiffnessY(int timestep) { return 1; };
	long double GetStiffnessZ(int timestep) { return 1; };
    */

	/*
	struct Getter {
		std::function<long double(TNode*)> coord;
		std::function<long double(TNode*)> ref;
		int type;
	};

	static long double getX(TNode* node) { return node->x; };
	static long double getY(TNode* node) { return node->y; };
	static long double getZ(TNode* node) { return node->z; };

	static long double getXRef(TNode* node) { return node->xRef; };
	static long double getYRef(TNode* node) { return node->yRef; };
	static long double getZRef(TNode* node) { return node->zRef; };

	static Getter Ox;
	static Getter Oy;
	static Getter Oz;
    */
	
	friend long double distance(TNode* P1, TNode* P2)
	{
		long double DX = P1->x - P2->x;
		long double DY = P1->y - P2->y;
		long double DZ = P1->z - P2->z;

		return sqrt(DX * DX + DY * DY + DZ * DZ);
	}
	

	/*
	static long double curvature(TNode* master, TNode* next, TNode* prev, TNode::Getter getters)
	{
		return getters.coord(next) - 2 * getters.coord(master) + getters.coord(prev);
	}

	virtual long double GetStretchingForce(TNode *neighbor, TNode::Getter getter) { return 0; };
	virtual long double GetBendingForce(TNode *next, TNode *prev, TNode::Getter getter) { return 0; };
	virtual long double GetTargetForce(TNode::Getter getter) { return 0; };
	*/
	TNode() {};
	~TNode() {};
};


class cImmersedBoundary
{
public:
	friend class cInputParameters;
	long double height;
	long double radius;
	long double step;
	long double y_begin_center;
	long double y_end_center;
	long double z_begin_center;
	long double z_end_center;
	int type;
	long double stiffness;
	long double stretchingStiffness;// = 10000;//10000;
	long double bendingStiffness;//= 0;//10000;//20; //1800;//15000;
	//long double step = 0.01;
	int nodesCount;
	int radius_nodes, height_nodes;
	TNode **nodes;
	int max_path_step, out_path_step;
	int n_pathes;
	btype **X_pathes;
	btype **Y_pathes;
	btype **Z_pathes;

	friend bool InputHole(int i, int j, int k, int N, int L, int M)
		// сделать дружественной погруженной границе
	{
		btype yc = L / 2;
		btype zc = M / 2;
		btype R = L / 2;
		if ((j - yc)*(j - yc) + (k - zc)*(k - zc) < R*R) return 1;
		else return 0;
	}

	cImmersedBoundary(cInputParameters* inpParams);
	~cImmersedBoundary(){};

private:
	void copyInputParams(sImmersedBoundaryParameters p);

};