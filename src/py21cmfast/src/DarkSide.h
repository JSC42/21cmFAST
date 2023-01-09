#include "ProfileTable.h5"
#define Integration_TimeStep 10000
#define Use_DarkSide 1
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

double Conditional_HMF(double growthf, double M, double Delta, double sigma2, struct CosmoParams *cosmo_params)
{
	// Return conditional HMF dn/dm, unit: 1/(Msun Mpc^3)
	// ---- inputs ----
	// growthf: growth factor
	// M: halo mass in Msun
	// Delta: Overdensity
	// sigma2: sigma2

	double OmegaMh2 = cosmo_params->OMm * pow(cosmo_params->hlittle, 2);
	double r;
	double RhoM = OmegaMh2 * 2.775E11; // Average matter density in Msun/Mpc^3
	float dNdM_old = dNdM_conditional((float)growthf, (float)log(M), 30.0, Deltac, (float)Delta, (float)sigma2);
	double dNdM_double = (double)dNdM_old;
	double RhoM_gird = RhoM * (1 + Delta);

	// This should be positive anyway, check where the - sign comes from
	return fabs(-RhoM_gird * dNdM_double / M / sqrt(2 * PI));
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

double FindBoost(double M1, double M2, double growthf, double z, double Delta, double sigma2, struct CosmoParams *cosmo_params, int hmf_model)
{
	// Return comoving emissivity at 21cm frequency for fiducial settings, in W/Hz/M^3
	// Using a brute force rectangular integrator for now, switch to GSL in next version
	// Set hmf_model = -1 to use conditional hmf
	// All masses are in Msun, integration is done in log space

	double lnm1, lnm2, dlnm, lnm, m, fun, hmf, result, ConversionFactor;
	int idx, TimeStep;
	ConversionFactor = 1.2111e-15;

	if (hmf_model >= 0)
	{
		TimeStep = 5 * Integration_TimeStep;
	}
	else
	{
		TimeStep = Integration_TimeStep;
	}

	// Don't waste time on out-range masses
	lnm1 = log(fmax(M1, pow(10,LgM_axis_min)));
	lnm2 = log(fmin(M2, pow(10,LgM_axis_max)));

	dlnm = (lnm2 - lnm1) / ((double)(TimeStep - 1));
	lnm = lnm1;
	result = 0.0;

	for (idx = 0; idx < TimeStep; idx++)
	{

		// The value of your function in requested mass
		m = exp(lnm);
		fun = InterpProfile(,m);

		if (hmf_model == -1)
		{
			hmf = Conditional_HMF(growthf, m, Delta, sigma2, cosmo_params);
		}
		else if (hmf_model == 0)
		{
			// PS HMF
			hmf = dNdM(growthf, m);
		}
		else if (hmf_model == 1)
		{
			// ST HMF
			hmf = dNdM_st(growthf, m);
		}
		else if (hmf_model == 2)
		{
			// Watson HMF
			hmf = dNdM_WatsonFOF(growthf, m);
		}
		else if (hmf_model == 3)
		{
			// Watson FoF HMF
			hmf = dNdM_WatsonFOF_z(z, growthf, m);
		}
		else
		{
			hmf = 0.0;
			LOG_ERROR("Wrong choice of HMF!");
			Throw(ValueError);
		}

		result = result + fun * hmf;

		lnm = lnm + dlnm;
	}

	result = result * dlnm;

	// Converting from Mpc to m
	return 1 + ConversionFactor * result/((pow(cosmo_params->OMm-cosmo_params->OMb),2) * (pow(cosmo_params->hlittle, 4)) pow(1+z,3));
}
