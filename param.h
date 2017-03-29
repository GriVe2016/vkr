//#include "geom.h"

#ifndef PARAMS
#define PARAMS

const int aneurysm = 1001;
const int thrombus = 1002;

struct sProblemParameters
{
	int type;
};

struct sDomainParameters
{
	long double h;
	long double lengthX, lengthY, lengthZ;

};

struct sAneurysmParameters
{
	long double x_center, y_center, z_center;
	long double radius;
	long double radiusAdd;
	long double stifnessMin;
	long double type;
};


struct sTimeParameters
{
	long double endValue, startValue, separationValue;

	void UpdateSeparationValue();

};

struct sImmersedBoundaryParameters
{
	long double height;
	long double radius;
	long double step;
	long double y_begin_center;
	long double y_end_center;
	long double z_begin_center;
	long double z_end_center;
	long double baseStifness;
	int type;
};

struct sPhysicalParameters
{
	long double Re;
	long double PressureInput;
	long double PressureOutput;
	long double visValue;
	long double soluteDensityValue;
	long double envDensityValue;

};

class cInputParameters
{
public:
	sDomainParameters domainParams;
	sTimeParameters timeParams;
	sImmersedBoundaryParameters immersedBoundaryParams;
	sPhysicalParameters physicalParams;
	sProblemParameters problemParams;
	sAneurysmParameters aneurysmParams;
	cInputParameters(char* inpfile);
};


#endif