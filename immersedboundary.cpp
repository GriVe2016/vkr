#include "stdafx.h"
#include "immersedboundary.h"
#include <math.h>


void cImmersedBoundary::copyInputParams(sImmersedBoundaryParameters p)
{
	height=p.height;
	radius=p.radius;
	step=p.step;
	y_begin_center=p.y_begin_center;
	y_end_center=p.y_end_center;
	z_begin_center=p.z_begin_center;
	z_end_center=p.z_end_center;
	type=p.type;
}

cImmersedBoundary::cImmersedBoundary(cInputParameters* inpParams)
{
	copyInputParams(inpParams->immersedBoundaryParams);
	//bool make_shift = true;
	//double shift_angle = 0.;
	//double shift_step = shift_angle*M_PI / 180.;
	//double shift = 0.0;
	long double iPart = 0.;
	long double ld_Pi = 4.*atanl(1.);

	modfl(this->height / this->step, &iPart);
	int endIndex = (int)iPart;
	this->height_nodes = endIndex + 1;

	modfl(this->radius*2.*ld_Pi / this->step, &iPart);
	endIndex = (int)iPart;
	this->radius_nodes = endIndex;
	this->nodesCount = this->height_nodes*this->radius_nodes;
	this->nodes = new TNode*[this->nodesCount];

		for (int n = 0; n < this->height_nodes; ++n) {
		for (int i = 0; i < this->radius_nodes; i++) {
			this->nodes[i + n*this->radius_nodes] = new TNode();
			//this->nodes[i + n * 120]->;
		}
	}
		/*
		if (inpParams->problemParams.type == aneurysm)
		{

			if (type == 1)
			{
				double rc;
				long double iRadius = (long double) this->radius_nodes;
				long double iHeight = (long double) this->height_nodes - 1.0;
				int cx = this->height_nodes / 2., cy = this->radius_nodes / 4., Rpip = this->radius_nodes / 10., ix, iy;

				for (int n = 0; n < this->height_nodes; ++n) {
					for (int i = 0; i < this->radius_nodes; i++) {
						long double x = this->height*double(n) / iHeight;
						long double y = this->y_begin_center + this->radius * cos(-ld_Pi * 2 * (iRadius - i) / iRadius); //+ shift_step
						long double z = this->z_begin_center + this->radius * sin(-ld_Pi * 2 * (iRadius - i) / iRadius); //+ shift_step

						this->nodes[i + n*this->radius_nodes]->x = x;
						this->nodes[i + n*this->radius_nodes]->xRef = x;
						this->nodes[i + n*this->radius_nodes]->xVel = 0;
						this->nodes[i + n*this->radius_nodes]->xForce = 0;

						this->nodes[i + n*this->radius_nodes]->y = y;
						this->nodes[i + n*this->radius_nodes]->yRef = y;
						this->nodes[i + n*this->radius_nodes]->yVel = 0;
						this->nodes[i + n*this->radius_nodes]->yForce = 0;

						this->nodes[i + n*this->radius_nodes]->z = z;
						this->nodes[i + n*this->radius_nodes]->zRef = z;
						this->nodes[i + n*this->radius_nodes]->zVel = 0;
						this->nodes[i + n*this->radius_nodes]->zForce = 0;


						ix = abs(n - cx);
						iy = abs(i - cy);
						rc = sqrt((double)(ix*ix + iy*iy));
						if (rc <= Rpip)
							this->nodes[i + n * this->radius_nodes]->stretchingStiffness = inpParams->aneurysmParams.stifnessMin; //+(this->stiffness - stiffness_plus)*sqrt(rc);    //sin(rc*0.5*M_PI / Rpip);
						else
							this->nodes[i + n *this->radius_nodes]->stretchingStiffness = this->stiffness;


					}
				}
			}
			else if (type == 2) // spline
				//--------------------------------  
			{
				char filename[] = "line.msh";
				FILE* SplFile;

				fopen_s(&SplFile, filename, "r"); // do exception !!!
				if (SplFile == NULL)
				{
					printf("No line.msh file.");
				}
				else
				{
					int nSpl;
					int tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8;

					fscanf_s(SplFile, "%d %d %d %d %d", &tmp1, &tmp2, &tmp3, &nSpl, &tmp4);
					double *xSpl;
					xSpl = new double[nSpl];
					double *zSpl;
					zSpl = new double[nSpl];
					double dtmp1;

					for (int i = 0; i < nSpl; i++)
					{
						fscanf_s(SplFile, "%d %d %d %d %d %d %d %d", &tmp1, &tmp2, &tmp3, &tmp4, &tmp5, &tmp6, &tmp7, &tmp8);
						fscanf_s(SplFile, "%lf %lf %lf", &xSpl[i], &zSpl[i], &dtmp1);
					}
					fclose(SplFile);

					//this->step = sqrtl(pow(xSpl[1] - xSpl[0], 2) + pow(zSpl[1] - zSpl[0], 2));
					//printf("\n==========================\n");
					//printf("this->step", this->step);

					double rc;
					long double iRadius = (long double) this->radius_nodes;
					long double iHeight = (long double) this->height_nodes - 1.;
					int cx = this->height_nodes / 2., cy = 3.*this->radius_nodes / 4., Rpip = this->radius_nodes / 15., ix, iy;
					int iSpl = 0;
					long double Alpha, z_plus;

					for (int n = 0; n < this->height_nodes; ++n) {
						long double x;
						if (n == 0)
						{
							x = xSpl[0];
							z_plus = zSpl[0];
							printf("-->z_b=%lf   ", z_plus);
						}
						else if (n == this->height_nodes - 1)
						{
							x = xSpl[nSpl - 1];
							z_plus = zSpl[nSpl - 1];
							printf("-->z_e=%lf   ", z_plus);
						}
						else // Finding intersection point
							//--------------------------------
						{
							x = this->height*float(n) / iHeight;
							while (!(xSpl[iSpl] <= x && x < xSpl[iSpl + 1]))
								iSpl++;
							Alpha = (x - xSpl[iSpl + 1]) / (xSpl[iSpl] - xSpl[iSpl + 1]);
							z_plus = zSpl[iSpl] * Alpha + (1. - Alpha)*zSpl[iSpl + 1];
							if (iSpl != 0) iSpl--;
						}
						for (int i = 0; i < this->radius_nodes; i++) {


							long double y = this->y_begin_center + this->radius * cos(-ld_Pi * 2 * (iRadius - i) / iRadius); // +shift_step);
							long double z = z_plus + this->radius * sin(-ld_Pi * 2 * (iRadius - i) / iRadius); // +shift_step);

							this->nodes[i + n*this->radius_nodes]->x = x;
							this->nodes[i + n*this->radius_nodes]->xRef = x;
							this->nodes[i + n*this->radius_nodes]->xVel = 0;
							this->nodes[i + n*this->radius_nodes]->xForce = 0;

							this->nodes[i + n*this->radius_nodes]->y = y;
							this->nodes[i + n*this->radius_nodes]->yRef = y;
							this->nodes[i + n*this->radius_nodes]->yVel = 0;
							this->nodes[i + n*this->radius_nodes]->yForce = 0;

							this->nodes[i + n*this->radius_nodes]->z = z;
							this->nodes[i + n*this->radius_nodes]->zRef = z;
							this->nodes[i + n*this->radius_nodes]->zVel = 0;
							this->nodes[i + n*this->radius_nodes]->zForce = 0;

							/*
							ix = abs(n - cx);
							iy = abs(i - cy);
							rc = sqrt((double)(ix*ix + iy*iy));
							if (rc <= Rpip)
							//this->nodes[i + n * 120]->stretchingStiffness = this->stiffness*rc*rc + sin(rc/Rpip)*stiffness_plus;
							this->nodes[i + n * this->radius_nodes]->stretchingStiffness = stiffness_plus; //+(this->stiffness - stiffness_plus)*sqrt(rc);    //sin(rc*0.5*M_PI / Rpip);
							else*/
							//this->nodes[i + n *this->radius_nodes]->stretchingStiffness = this->stiffness;

							// Verison 1 (not good)
							/*
							if ((n>=49 && n<=69) && (i>=20 && i<=39))
							this->nodes[i + n * 120]->stretchingStiffness = this->stiffness/100;
							else
							this->nodes[i + n * 120]->stretchingStiffness=this->stiffness;
							*/
						//}
						//this->zCenter += 0.005;
					//}
					//this->zCenter -= 0.005;
					//printf("\nzCenter=%lf\n", this->zCenter);
				//}
/*

			}
			else if (type == 3)
			{
				//printf("\nBOUNDARY TYPE = %1d step=%lf\n", tube_type, shift_step);


				double an_begin = 0.25;
				double an1 = 0.35, an2 = 0.65;
				double an_end = 0.75;
				double step_div = (an1 - an_begin) / this->step;
				double Alpha_step = 1. / step_div;
				double Alpha = 1.;
				bool Alpha_flag = true;

				double rc, stiffness_plus = this->stiffness / 30;
				double iRadius = (double) this->radius_nodes;
				double iHeight = (double) this->height_nodes - 1.;
				//int cx = this->height_nodes / 2., cy = this->radius_nodes / 4., Rpip = this->radius_nodes / 10., ix, iy;

				for (int n = 0; n < this->height_nodes; ++n)
				{
					//if (n != 0) shift += shift_step;
					//printf("\nn=%d shift=%lf ", n, shift);

					for (int i = 0; i < this->radius_nodes; i++)
					{
						long double x = this->height*float(n) / iHeight;
						long double y = this->y_begin_center + this->radius * cos(-ld_Pi * 2 * (iRadius - i) / iRadius); // +shift);
						long double z = this->z_begin_center + this->radius * sin(-ld_Pi * 2 * (iRadius - i) / iRadius); // +shift);
						if (i == 0)
						{
							printf("  --> %lf %lf  s=%lf ", y, z);
							if ((n + 1) % 2 == 0) printf("\n");
						}

						this->nodes[i + n*this->radius_nodes]->x = x;
						this->nodes[i + n*this->radius_nodes]->xRef = x;
						this->nodes[i + n*this->radius_nodes]->xVel = 0;
						this->nodes[i + n*this->radius_nodes]->xForce = 0;

						this->nodes[i + n*this->radius_nodes]->y = y;
						this->nodes[i + n*this->radius_nodes]->yRef = y;
						this->nodes[i + n*this->radius_nodes]->yVel = 0;
						this->nodes[i + n*this->radius_nodes]->yForce = 0;

						this->nodes[i + n*this->radius_nodes]->z = z;
						this->nodes[i + n*this->radius_nodes]->zRef = z;
						this->nodes[i + n*this->radius_nodes]->zVel = 0;
						this->nodes[i + n*this->radius_nodes]->zForce = 0;


						//ix = abs(n - cx);
						//iy = abs(i - cy);
						//rc = sqrt((double)(ix*ix + iy*iy));
						if (x >= an_begin && x <= an1)
						{
							this->nodes[i + n * this->radius_nodes]->stretchingStiffness = Alpha*this->stiffness + (1. - Alpha)*stiffness_plus; //+(this->stiffness - stiffness_plus)*sqrt(rc);    //sin(rc*0.5*M_PI / Rpip);
							if (i == this->radius_nodes - 1)
							{
								Alpha -= Alpha_step;
								printf("-->%lf  ", this->nodes[i + n * this->radius_nodes]->stretchingStiffness);
							}
						}
						else if (x > an1 && x < an2)
						{
							this->nodes[i + n * this->radius_nodes]->stretchingStiffness = stiffness_plus;
							if (Alpha_flag && i == this->radius_nodes - 1)
							{
								Alpha_flag = false;
								Alpha = 0.;
								printf("-->%lf  ", this->nodes[i + n * this->radius_nodes]->stretchingStiffness);
							}
						}
						else if (x >= an2 && x <= an_end)
						{
							this->nodes[i + n * this->radius_nodes]->stretchingStiffness = Alpha*this->stiffness + (1. - Alpha)*stiffness_plus; //+(this->stiffness - stiffness_plus)*sqrt(rc);    //sin(rc*0.5*M_PI / Rpip);
							if (i == this->radius_nodes - 1)
							{
								Alpha += Alpha_step;
								printf("-->%lf  ", this->nodes[i + n * this->radius_nodes]->stretchingStiffness);
							}
						}
						else
							this->nodes[i + n *this->radius_nodes]->stretchingStiffness = this->stiffness;


					}
				}
			}

			else if (type == 4) // spline
				//--------------------------------  
			{
				char filename[] = "line.msh";
				FILE* SplFile;

				fopen_s(&SplFile, filename, "r"); // do exception !!!
				if (SplFile == NULL)
				{
					printf("No line.msh file.");
				}
				else
				{
					int nSpl;
					int tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8;

					fscanf_s(SplFile, "%d %d %d %d %d", &tmp1, &tmp2, &tmp3, &nSpl, &tmp4);
					double *xSpl;
					xSpl = new double[nSpl];
					double *ySpl;
					ySpl = new double[nSpl];
					double *zSpl;
					zSpl = new double[nSpl];
					double dtmp1;

					for (int i = 0; i <= nSpl - 1; i++)
					{
						fscanf_s(SplFile, "%d %d %d %d %d %d %d %d", &tmp1, &tmp2, &tmp3, &tmp4, &tmp5, &tmp6, &tmp7, &tmp8);
						fscanf_s(SplFile, "%lf %lf %lf", &xSpl[i], &zSpl[i], &dtmp1);
						//zSpl[i] += 1.0;
						ySpl[i] = 0.4;
						printf("\n %2d  %lf %lf %lf", i, xSpl[i], ySpl[i], zSpl[i]);
					}
					fclose(SplFile);

					//this->step = sqrtl(pow(xSpl[1] - xSpl[0], 2) + pow(zSpl[1] - zSpl[0], 2));
					//printf("\n==========================\n");
					//printf("this->step", this->step);
					long double iPart = 0.;
					this->step = 0.01;
					this->height = 1.0;
					modfl(this->height / this->step, &iPart);
					int endIndex = (int)iPart;
					this->height_nodes = endIndex + 1;


					//this->height_nodes = nSpl;

					double rc, stiffness_plus = this->stiffness / 28;
					double iRadius = (double) this->radius_nodes;
					double iHeight = (double) this->height_nodes - 1.;
					double dcx = 0.5;  //this->height_nodes / 2.;
					double dcy = 0.4;  //3.*this->radius_nodes / 4.;
					double dcz = 0.65;
					double dRmax = 0.25, dRmin = 0.12;
					double dRpip = dRmax;  //this->radius_nodes / 15.;

					int ix, iy;
					int iSpl = 0;
					long double Alpha, y_plus, z_plus;

					for (int n = 0; n < this->height_nodes; ++n) {
						long double x;
						if (n == 0)
						{
							x = xSpl[0];
							y_plus = ySpl[0];
							z_plus = zSpl[0];
							printf("-->z_b=%lf   ", z_plus);
						}
						else if (n == this->height_nodes - 1)
						{
							x = xSpl[nSpl - 1];
							y_plus = ySpl[nSpl - 1];
							z_plus = zSpl[nSpl - 1];
							printf("-->z_e=%lf   ", z_plus);
						}
						else // Finding intersection point
							//--------------------------------
						{
							x = this->height*float(n) / iHeight;
							while (!(xSpl[iSpl] <= x && x < xSpl[iSpl + 1]))
								iSpl++;
							Alpha = (x - xSpl[iSpl + 1]) / (xSpl[iSpl] - xSpl[iSpl + 1]);
							y_plus = ySpl[iSpl] * Alpha + (1. - Alpha)*ySpl[iSpl + 1];
							z_plus = zSpl[iSpl] * Alpha + (1. - Alpha)*zSpl[iSpl + 1];
							if (iSpl != 0) iSpl--;
						}
						for (int i = 0; i < this->radius_nodes; i++) {


							long double y = y_plus + this->radius * cos(-ld_Pi * 2 * (iRadius - i) / iRadius); // +shift_step);
							long double z = z_plus + this->radius * sin(-ld_Pi * 2 * (iRadius - i) / iRadius); // +shift_step);

							this->nodes[i + n*this->radius_nodes]->x = x;
							this->nodes[i + n*this->radius_nodes]->xRef = x;
							this->nodes[i + n*this->radius_nodes]->xVel = 0;
							this->nodes[i + n*this->radius_nodes]->xForce = 0;

							this->nodes[i + n*this->radius_nodes]->y = y;
							this->nodes[i + n*this->radius_nodes]->yRef = y;
							this->nodes[i + n*this->radius_nodes]->yVel = 0;
							this->nodes[i + n*this->radius_nodes]->yForce = 0;

							this->nodes[i + n*this->radius_nodes]->z = z;
							this->nodes[i + n*this->radius_nodes]->zRef = z;
							this->nodes[i + n*this->radius_nodes]->zVel = 0;
							this->nodes[i + n*this->radius_nodes]->zForce = 0;


							//ix = x - dcx;
							//iy = abs(i - cy);
							rc = sqrt((x - dcx)*(x - dcx) + (y - dcy)*(y - dcy) + (z - dcz)*(z - dcz));

							if (rc <= dRmax)
							{
								//Count Alpha
								double Alpha_rad = (rc - dRmax) / (dRmin - dRmax);
								if (0.0 <= Alpha_rad && Alpha_rad <= 1.0)
								{
									this->nodes[i + n * this->radius_nodes]->stretchingStiffness = (1.0 - Alpha_rad)*this->stiffness + Alpha_rad*stiffness_plus;
								}
								else
									this->nodes[i + n * this->radius_nodes]->stretchingStiffness = stiffness_plus;

							}
							else
								this->nodes[i + n *this->radius_nodes]->stretchingStiffness = this->stiffness;
							if (i == 0) printf("s(%2d)=%lf  ", n, this->nodes[i + n *this->radius_nodes]->stretchingStiffness);

							// Verison 1 (not good)
							/*
							if ((n>=49 && n<=69) && (i>=20 && i<=39))
							this->nodes[i + n * 120]->stretchingStiffness = this->stiffness/100;
							else
							this->nodes[i + n * 120]->stretchingStiffness=this->stiffness;
							*/
						//}
						//this->zCenter += 0.005;
					//}
					//this->zCenter -= 0.005;
					//printf("\nzCenter=%lf\n", this->zCenter);
				//}

			//}
/*
		} 

	for (int n = 0; n < this->height_nodes; ++n) {
		for (int i = 0; i < this->radius_nodes; i++) {
			//this->nodes[i + n*this->radius_nodes]->boundary = this;

			if (n <  this->height_nodes - 1)
			{
				this->nodes[i + n*this->radius_nodes]->neighbors.next = this->nodes[i + (n + 1) * this->radius_nodes];
				this->nodes[i + n*this->radius_nodes]->neighbors.initialDistanceNext =
					distance(this->nodes[i + n * this->radius_nodes], this->nodes[i + n * this->radius_nodes]->neighbors.next);
			}
			if (n > 0)
			{
				this->nodes[i + n*this->radius_nodes]->neighbors.prev = this->nodes[i + (n - 1) * this->radius_nodes];
				this->nodes[i + n*this->radius_nodes]->neighbors.initialDistancePrev = distance(this->nodes[i + n * this->radius_nodes], this->nodes[i + n * this->radius_nodes]->neighbors.prev);
			}
			/*
			if (this->nodes[i + n*this->radius_nodes]->neighbors.prev != nullptr && this->nodes[i + n * this->radius_nodes]->neighbors.next != nullptr)
			{
				this->nodes[i + n*this->radius_nodes]->neighbors.initialCurvatures[OX] = TNode::curvature(
					this->nodes[i + n*this->radius_nodes],
					this->nodes[i + n*this->radius_nodes]->neighbors.next,
					this->nodes[i + n*this->radius_nodes]->neighbors.prev,
					TNode::Ox
					);
				this->nodes[i + n*this->radius_nodes]->neighbors.initialCurvatures[OY] = TNode::curvature(
					this->nodes[i + n*this->radius_nodes],
					this->nodes[i + n*this->radius_nodes]->neighbors.next,
					this->nodes[i + n*this->radius_nodes]->neighbors.prev,
					TNode::Oy
					);
				this->nodes[i + n*this->radius_nodes]->neighbors.initialCurvatures[OZ] = TNode::curvature(
					this->nodes[i + n*this->radius_nodes],
					this->nodes[i + n*this->radius_nodes]->neighbors.next,
					this->nodes[i + n*this->radius_nodes]->neighbors.prev,
					TNode::Oz
					);
			}
			*/
		//}
	//}
	/*
	this->n_pathes = 6;
	this->max_path_step = 1000;
	this->out_path_step = 30;

	this->X_pathes = new double*[this->n_pathes];
	this->Y_pathes = new double*[this->n_pathes];
	this->Z_pathes = new double*[this->n_pathes];

	for (int i = 0; i < this->n_pathes; i++)
	{
		this->X_pathes[i] = new double[this->max_path_step];
		this->Y_pathes[i] = new double[this->max_path_step];
		this->Z_pathes[i] = new double[this->max_path_step];
	}

	for (int i = 0; i < this->n_pathes; i++)
	for (int j = 0; j < this->max_path_step; j++)
	{
		this->X_pathes[i][j] = 0.0;
		this->Y_pathes[i][j] = 0.0;
		this->Z_pathes[i][j] = 0.0;
	}

	this->Y_pathes[0][0] = 0.4;
	this->Z_pathes[0][0] = 0.495;

	this->Y_pathes[1][0] = 0.448;
	this->Z_pathes[1][0] = 0.482;

	this->Y_pathes[2][0] = 0.482;
	this->Z_pathes[2][0] = 0.448;

	this->Y_pathes[3][0] = 0.495;
	this->Z_pathes[3][0] = 0.4;

	this->Y_pathes[4][0] = 0.467;
	this->Z_pathes[4][0] = 0.333;

	this->Y_pathes[5][0] = 0.4;
	this->Z_pathes[5][0] = 0.305;
	*/
}
