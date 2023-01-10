#include "ProfileTable.h"
#include<stdio.h>
#include<math.h>

double Interp_Fast(double *Tab, double xmin, double xmax, int nx, double x)
{
	// Interpolate, x axis must be linear in [xmin xmax] with size nx
	// Speed: 1E-8 seconds per call for nx=100

	double y1, y2, y, dx, x1, x2;
	int id1, id2;

	dx = (xmax - xmin) / ((double)(nx - 1));
	id1 = (int)floor((x - xmin) / dx);
	id2 = id1 + 1;
	x1 = xmin + ((double)id1) * dx;
	x2 = x1 + dx;

	if ((id1 < 0) || (id1 > nx - 2))
	{
		// Do not extrapolate
		y = 0;
	}
	else
	{
		y1 = Tab[id1];
		y2 = Tab[id2];
		y = (y2 - y1) * (x - x1) / (x2 - x1) + y1;
	}

	return y;
}

double InterpProfile(double z, double m)
{
	int zid1, zid2, mid1, mid2;
	double z1, z2, m1, m2, lnm1, lnm2, r, dz, dlnm, lnm_min, lnm_max, lnm, f1, f2, F1, F2;

	lnm = log(m);

	lnm_min = 2.3026 * LgM_axis_min;
	lnm_max = 2.3026 * LgM_axis_max;
	
	dz = (Z_axis_max - Z_axis_min) / ((double)Z_axis_Size - 1);
	dlnm = (lnm_max - lnm_min) / ((double)M_axis_Size - 1);

	zid1 = (int)floor((z - Z_axis_min) / dz);
	mid1 = (int)floor((lnm - lnm_min) / dlnm);
	zid2 = zid1 + 1;
	mid2 = mid1 + 1;

	if ((zid1 < 0) || (zid1 > Z_axis_Size - 2) || (mid1 < 0) || (mid1 > M_axis_Size - 2))
	{
		// Do not extrapolate
		r = 0.0;
	}
	else
	{

		z1 = Z_axis_min + ((double)zid1) * dz;
		z2 = Z_axis_min + ((double)zid2) * dz;

		lnm1 = lnm_min + ((double)mid1) * dlnm;
		lnm2 = lnm_min + ((double)mid2) * dlnm;
		printf("%E\n",lnm1);
		printf("%E\n",lnm);
		printf("%E\n",lnm2);
		printf("%E\n",z1);
		printf("%E\n",z);
		printf("%E\n",z2);
		
		// Interpolate for lnm1
		f1 = ProfileTable[zid1][mid1];
		f2 = ProfileTable[zid2][mid1];

		F1 = (f2 - f1) * (z - z1) / (z2 - z1) + f1;

		// Interpolate for lnm2
		f1 = ProfileTable[zid1][mid2];
		f2 = ProfileTable[zid2][mid2];

		F2 = (f2 - f1) * (z - z1) / (z2 - z1) + f1;

		r = (F2 - F1) * (lnm - lnm1) / (lnm2 - lnm1) + F1;
	}

	// Convert erg to Joule unit
	return r;
}

int main()
{
	double m,z,p;
	// m=100;
	// z=190;
	scanf("%lf",&m);
	scanf("%lf",&z);
	p=InterpProfile(z,m);
	printf("%E\n",p);
	return 0;
}