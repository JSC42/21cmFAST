#define Integration_TimeStep 2000
#define PBH_Table_Size 100

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
