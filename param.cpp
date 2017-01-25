#include "stdafx.h"
#include "param.h"
#include "math.h"
#include "immersedboundary.h"
#include <iostream>
#include <fstream>


using namespace std;

cInputParameters::cInputParameters(char* inpfile)
{
	FILE *mf;
	int  forComment = 20;
	// Строки для считывания о обработки ошибок
	char str[80];
	char *estr;

	try
	{

		// Открытие файла с режимом доступа «только чтение» и привязка к нему потока данных
		printf(" Open data file: ");

		 fopen_s(&mf,inpfile, "r");

		// Проверка открытия файла
		if (mf == NULL) { throw inpfile; }
		

		//Чтение (построчно) данных из файла в бесконечном цикле 
		//while (1)
		//{
		// Чтение строк  из файла
		// ProblemParameters
		//---------------------------------
		estr = fgets(str, sizeof(str), mf); //1
		puts(str);
		//---------------------------------
		estr = fgets(str, sizeof(str), mf); //2
		if (estr == NULL) throw 2;
		puts(str);
		problemParams.type = atoi(&str[forComment]);
		printf(" %d\n", problemParams.type);
		// DomainParameters
		//---------------------------------
		estr = fgets(str, sizeof(str), mf); //1
		puts(str);
		//---------------------------------
		estr = fgets(str, sizeof(str), mf); //2
		if (estr==NULL) throw 2;
		puts(str);
		domainParams.h = atof(&str[forComment]);
		printf(" %lf\n", domainParams.h);
		//--------------------------------------
		estr = fgets(str, sizeof(str), mf); //3
		if (estr == NULL) throw 3;
		puts(str);
		domainParams.lengthX = atof(&str[forComment]);
		printf(" %lf\n", domainParams.lengthX);
		//--------------------------------------
		estr = fgets(str, sizeof(str), mf); //4
		if (estr == NULL) throw 4;
		puts(str);
		domainParams.lengthY = atof(&str[forComment]);
		printf(" %lf\n", domainParams.lengthY);
		//--------------------------------------
		estr = fgets(str, sizeof(str), mf);//5
		puts(str);
		domainParams.lengthZ = atof(&str[forComment]);
		printf(" %lf\n", domainParams.lengthZ);

		// sTimeParameters
		//---------------------------------
		estr = fgets(str, sizeof(str), mf);
		puts(str);
		//--------------------------------------
		estr = fgets(str, sizeof(str), mf);
		puts(str);
		timeParams.startValue = atof(&str[forComment]);
		printf(" %lf\n", timeParams.startValue);
		//--------------------------------------
		estr = fgets(str, sizeof(str), mf);
		puts(str);
		timeParams.endValue = atof(&str[forComment]);
		printf(" %lf\n", timeParams.endValue);
		//--------------------------------------
		estr = fgets(str, sizeof(str), mf);
		puts(str);
		timeParams.separationValue = atof(&str[forComment]);
		printf(" %lf\n", timeParams.separationValue);

		// sImmersedBoundaryParameters
		//--------------------------------------------------
		estr = fgets(str, sizeof(str), mf);
		puts(str);
		//--------------------------------------
		estr = fgets(str, sizeof(str), mf);
		puts(str);
		immersedBoundaryParams.height = atof(&str[forComment]);
		printf(" %lf\n", immersedBoundaryParams.height);
		//--------------------------------------
		estr = fgets(str, sizeof(str), mf);
		puts(str);
		immersedBoundaryParams.radius = atof(&str[forComment]);
		printf(" %lf\n", immersedBoundaryParams.radius);
		//--------------------------------------
		estr = fgets(str, sizeof(str), mf);
		puts(str);
		immersedBoundaryParams.step = atof(&str[forComment]);
		printf(" %lf\n", immersedBoundaryParams.step);
		//--------------------------------------
		estr = fgets(str, sizeof(str), mf);
		puts(str);
		immersedBoundaryParams.y_begin_center = atof(&str[forComment]);
		printf(" %lf\n", immersedBoundaryParams.y_begin_center);
		//--------------------------------------
		estr = fgets(str, sizeof(str), mf);
		puts(str);
		immersedBoundaryParams.z_begin_center = atof(&str[forComment]);
		printf(" %lf\n", immersedBoundaryParams.z_begin_center);
		//--------------------------------------
		estr = fgets(str, sizeof(str), mf);
		puts(str);
		immersedBoundaryParams.y_end_center = atof(&str[forComment]);
		printf(" %lf\n", immersedBoundaryParams.y_end_center);
		//--------------------------------------
		estr = fgets(str, sizeof(str), mf);
		puts(str);
		immersedBoundaryParams.z_end_center = atof(&str[forComment]);
		printf(" %lf\n", immersedBoundaryParams.z_end_center);
		//--------------------------------------
		estr = fgets(str, sizeof(str), mf);
		puts(str);
		immersedBoundaryParams.type = atoi(&str[forComment]);
		printf(" %d\n", immersedBoundaryParams.type);

		// sPhysicalParameters
		//---------------------------------------------------
		estr = fgets(str, sizeof(str), mf);
		puts(str);
		//---------------------------------------------------
		estr = fgets(str, sizeof(str), mf);
		puts(str);
		physicalParams.Re = atof(&str[forComment]);
		printf(" %lf\n", physicalParams.Re);
		//---------------------------------------------------
		estr = fgets(str, sizeof(str), mf);
		puts(str);
		physicalParams.PressureInput = atof(&str[forComment]);
		printf(" %lf\n", physicalParams.PressureInput);
		//---------------------------------------------------
		estr = fgets(str, sizeof(str), mf);
		puts(str);
		physicalParams.PressureOutput = atof(&str[forComment]);
		printf(" %lf\n", physicalParams.PressureOutput);
		//---------------------------------------------------
		estr = fgets(str, sizeof(str), mf);
		puts(str);
		long double viscisityCoef;
		viscisityCoef = atof(&str[forComment]);
		printf(" %lf\n", viscisityCoef);
		physicalParams.visValue = viscisityCoef / physicalParams.Re;
		//---------------------------------------------------
		estr = fgets(str, sizeof(str), mf);
		puts(str);
		physicalParams.soluteDensityValue = atof(&str[forComment]);
		printf(" %lf\n", physicalParams.soluteDensityValue);
		//---------------------------------------------------
		estr = fgets(str, sizeof(str), mf);
		puts(str);
		physicalParams.envDensityValue = atof(&str[forComment]);
		printf(" %lf\n", physicalParams.envDensityValue);

		// sAneurythmParameters
		//---------------------------------------------------
		estr = fgets(str, sizeof(str), mf);
		puts(str);
		//---------------------------------------------------
		estr = fgets(str, sizeof(str), mf);
		puts(str);
		aneurysmParams.x_center = atof(&str[forComment]);
		printf(" %lf\n", aneurysmParams.x_center);
		//---------------------------------------------------
		estr = fgets(str, sizeof(str), mf);
		puts(str);
		aneurysmParams.y_center = atof(&str[forComment]);
		printf(" %lf\n", aneurysmParams.y_center);
		//---------------------------------------------------
		estr = fgets(str, sizeof(str), mf);
		puts(str);
		aneurysmParams.z_center = atof(&str[forComment]);
		printf(" %lf\n", aneurysmParams.z_center);
		//---------------------------------------------------
		estr = fgets(str, sizeof(str), mf);
		puts(str);
		aneurysmParams.radius = atof(&str[forComment]);
		printf(" %lf\n", aneurysmParams.radius);
		//---------------------------------------------------
		estr = fgets(str, sizeof(str), mf);
		puts(str);
		aneurysmParams.radiusAdd = atof(&str[forComment]);
		printf(" %lf\n", aneurysmParams.radiusAdd);
		//---------------------------------------------------
		estr = fgets(str, sizeof(str), mf);
		puts(str);
		aneurysmParams.stifnessMin = atof(&str[forComment]);
		printf(" %lf\n", aneurysmParams.stifnessMin);
		//---------------------------------------------------
		estr = fgets(str, sizeof(str), mf);
		puts(str);
		aneurysmParams.type = atoi(&str[forComment]);
		printf(" %lf\n", aneurysmParams.type);

		if (fclose(mf) == EOF) printf("ошибка\n");
 		else printf("выполнено\n");
	}
	catch (char* inpfile)
	{
		cerr << "can't open input file "<<inpfile<<"\n";
	}
	catch (int line)
	{
		cerr << "reading input file error line"<< line<<"\n";
	}
	//timeParams.UpdateSeparationValue();
}

void sTimeParameters::UpdateSeparationValue()
{
	long double resDev = (endValue - startValue) / separationValue;
	long double iPart = 0;

	modfl(resDev, &iPart);

	int fEndIndex = (int)iPart;
	long double eps = 0.00000000001;
	if (fabsl(resDev - iPart) > eps) { fEndIndex++; }

	separationValue = (endValue - startValue) / ((long double)fEndIndex);

}