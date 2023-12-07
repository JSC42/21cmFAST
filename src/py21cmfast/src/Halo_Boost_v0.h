// Halo Boost Factor module, can be run and testes outside 21cmFAST

// #include <stdio.h>
// #include <stdlib.h>
// #include <math.h>

#define HaloProfile_Length 2000
#define Boost_nmh 1000
#define Boost_Interp_Table_Size 100

double HaloProfile_Kernel(double z, double mh, double r, int ProfileType)
{
    /* Halo density profile
    ---- inputs ----
    z : redshift
    mh : halo mass in msun, dm+baryon
    r : distance to center in pc
    ProfileType : result type
                  0 - viral radius in pc
                  1 - Rho_DM in msun/pc^3
    */

    // printf("Check that the integrated mass converges to m\n");
    double OmM, OmC, m, OmR, OmL, h, pi, rho_cr0, zp, OmMz, d, Delta_C, log10_c, c, delta_c, rv1, rv2, rv3, r_vir, x, cx, rho_cr, RhoDM;

    // Some settings
    OmM = 0.30964168161;
    OmR = 9.1e-5;
    OmC = 0.260667;
    OmL = 0.69026731839;
    h = 0.6766;
    pi = 3.141592653589793;
    rho_cr0 = 2.775e-7 * pow(h, 2.); // critical density in msun/pc^3

    // Pre-requisites
    m = mh * OmC / OmM; // DM mass
    zp = 1. + z;
    OmMz = OmM * pow(zp, 3.) / (OmM * pow(zp, 3.) + OmL);
    d = OmMz - 1.;
    Delta_C = 18. * pow(pi, 2.) + 82. * d - 39. * pow(d, 2.);
    log10_c = 1.071 - 0.098 * (log10(m) - 12.);
    c = pow(10., log10_c) / zp; // concentration, see appdx.A of Zip.et for the additional (1+z) factor
    delta_c = Delta_C * pow(c, 3.) / (3. * (log(1. + c) - c / (1. + c)));

    rv1 = 0.784 * pow(m * h / 1.0e8, 1. / 3.);
    rv2 = pow(OmM * Delta_C / (OmMz * 18. * pow(pi, 2.)), -1. / 3.);
    rv3 = (10. / zp) / h * 1000.;
    r_vir = rv1 * rv2 * rv3;

    x = r / r_vir;
    if (x > 1.)
    {
        RhoDM = 0.;
        // printf("%E  %E\n", r, r_vir);
    }
    else
    {
        cx = c * x;
        // rho_cr = rho_cr0 * (OmL + OmM * zp**3 + OmR * zp**4);
        rho_cr = rho_cr0 * (OmL + OmM * pow(zp, 3.) + OmR * pow(zp, 4.));
        RhoDM = rho_cr * delta_c / (cx * pow(1. + cx, 2.));
    }
    if (ProfileType == 0)
    {
        return r_vir;
    }
    else
    {
        return RhoDM;
    }
}

double HaloProfile_Integrator_analytic(double z, double mh)
{
    /* Do the following HaloProfile integration analytically:
    \int dr r^2 \rho^2_{dm}
    ---- inputs ----
    z : redshift
    mh : halo mass in msun, dm+baryon
    ---- output unit ----
    msun^2/Mpc^3
    */

    // printf("Check that the integrated mass converges to m\n");
    double OmM, OmC, m, OmR, OmL, h, pi, rho_cr0, zp, OmMz, d, Delta_C, log10_c, c, delta_c, rv1, rv2, rv3, r_vir, rho_cr, res;

    // Some settings
    OmM = 0.30964168161;
    OmR = 9.1e-5;
    OmC = 0.260667;
    OmL = 0.69026731839;
    h = 0.6766;
    pi = 3.141592653589793;
    rho_cr0 = 2.775e-7 * pow(h, 2.); // critical density in msun/pc^3

    // Pre-requisites
    m = mh * OmC / OmM; // DM mass
    zp = 1. + z;
    OmMz = OmM * pow(zp, 3.) / (OmM * pow(zp, 3.) + OmL);
    d = OmMz - 1.;
    Delta_C = 18. * pow(pi, 2.) + 82. * d - 39. * pow(d, 2.);
    log10_c = 1.071 - 0.098 * (log10(m) - 12.);
    c = pow(10., log10_c) / zp; // concentration, see appdx.A of Zip.et for the additional (1+z) factor
    delta_c = Delta_C * pow(c, 3.) / (3. * (log(1. + c) - c / (1. + c)));

    rv1 = 0.784 * pow(m * h / 1.0e8, 1. / 3.);
    rv2 = pow(OmM * Delta_C / (OmMz * 18. * pow(pi, 2.)), -1. / 3.);
    rv3 = (10. / zp) / h * 1000.;
    r_vir = rv1 * rv2 * rv3;                                          // pc
    rho_cr = rho_cr0 * (OmL + OmM * pow(zp, 3.) + OmR * pow(zp, 4.)); // msun/pc^3

    res = pow(rho_cr * delta_c, 2.0);
    res = res * pow(r_vir / c / (1. + c), 3.0);
    res = res / 3.0 * (pow(1. + c, 3.) - 1.0);
    res = res / 1.0e18; // convert to msun^2/Mpc^3

    return res;
}

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

double diff(double x1, double x2)
{
    // Relative difference between 2 non-zero numbers
    double d1, d2, r;
    d1 = (x2 - x1) / x1;
    d2 = (x2 - x1) / x2;
    return fabs(fmax(d1, d2));
}

double Integrate(double *x, double *dydx, int nx)
{
    int idx;
    double dx, r;

    r = 0.0;
    for (idx = 1; idx < nx; idx++)
    {
        dx = x[idx] - x[idx - 1];
        r = r + dx * dydx[idx];
    }

    return r;
}

void linspace_vec(double x1, double x2, int nx, double *x)
{
    // Fill array x with numbers linearly distributed between [x1, x2], nx is length
    double dx, x_;
    int idx;
    dx = (x2 - x1) / ((double)nx - 1.0);
    x_ = x1;
    for (idx = 0; idx < nx; idx++)
    {
        x[idx] = x_;
        x_ = x_ + dx;
    }
}

void logspace_vec(double x1, double x2, int nx, double *x)
{
    double lx1, lx2, lx, dlx;
    int idx;
    lx1 = log(x1);
    lx2 = log(x2);

    dlx = (lx2 - lx1) / ((double)nx - 1.0);
    lx = lx1;

    for (idx = 0; idx < nx; idx++)
    {
        x[idx] = exp(lx);
        lx = lx + dlx;
    }
}

double HaloProfile_vec(double z, double mh, double r_min_ratio, double *r, double *rho, int nr)
{
    double r_max, r_, pi;
    int idx;
    double fun[nr], lnr[nr], converge;

    r_max = HaloProfile_Kernel(z, mh, 10.0, 0);
    logspace_vec(r_min_ratio * r_max, r_max, nr, r);

    for (idx = 0; idx < nr; idx++)
    {
        r_ = r[idx];
        rho[idx] = HaloProfile_Kernel(z, mh, r_, 1);

        lnr[idx] = log(r_);
        // The thing that we actually need
        fun[idx] = pow(r_, 3.0) * pow(rho[idx], 2.0);
    }

    converge = Integrate(lnr, fun, nr);

    return converge;
}

double HaloProfile(double z, double mh, double *r, double *rho, int nr)
{
    // Get a converged Halo Profile, unit : Msun/pc^3
    // This is not very efficient because a huge part of the calls might be wasted

    int count, count_max, proceed, idx;
    double dif, tol, r_min_ratio, c1, c2, lnr[nr], mfun[nr], pi, r_, mass;

    if (nr > 100000)
    {
        // Foolproof measure
        printf("nr too large, exitting");
        exit(1);
    }

    // Stop when each iteration changes integration by less than tol

    tol = 0.05;
    proceed = 1;
    count = 0;
    count_max = 20;
    r_min_ratio = 0.1;
    pi = 3.141592653589793;

    c1 = 0.0;

    while (proceed)
    {
        c2 = HaloProfile_vec(z, mh, r_min_ratio, r, rho, nr);
        dif = fabs(c2 - c1) / c2;

        // printf("count = %d, dif = %E, r1 = %E, r2 = %E\n", count, dif, r[0], r[nr-1]);
        if (dif < tol)
        {
            proceed = 0;
        }
        r_min_ratio = r_min_ratio / 10.0;
        count = count + 1;
        if (count > count_max)
        {
            printf("count_max exceeded!");
            proceed = 0;
        }

        // numerical things?
        if (c2 < c1)
        {
            printf("c2 should be increasing, perhaps due to resolution issue, try increasing nr, e.g. in HaloProfile_Integrator. Exiting.\n");
            printf("c1 = %E, c2 = %E, count = %d\n", c1, c2, count);
            exit(1);
        }

        c1 = c2;
    }

    // While we are at it, compute mass and see whether it agrees with input
    for (idx = 0; idx < nr; idx++)
    {
        r_ = r[idx];
        lnr[idx] = log(r_);
        mfun[idx] = 4 * pi * rho[idx] * pow(r_, 3.0);
    }

    mass = Integrate(lnr, mfun, nr);

    return mass;
}

double HaloProfile_Integrator(double z, double mh)
{
    // Do integration \int dr r^2 rho^2_{dm}
    // result unit: msun^2/Mpc^3
    int nr = HaloProfile_Length, idx;
    double r[nr], RhoM[nr], fun[nr], lnr[nr], result, OmM, OmC, r_, RhoC, swap;

    OmM = 0.30964168161;
    OmC = 0.260667;

    swap = HaloProfile(z, mh, r, RhoM, nr);

    for (idx = 0; idx < nr; idx++)
    {
        r_ = r[idx];
        lnr[idx] = log(r_);
        RhoC = RhoM[idx] * OmC / OmM;
        fun[idx] = pow(r_, 3.0) * pow(RhoC, 2.0);
    }

    result = Integrate(lnr, fun, nr) * 1.0e18;

    return result;
}

double Boost_Factor_v2(double z, double delta, double sigma2, struct CosmoParams *cosmo_params, int hmf_model)
{
    // Get the boost factor all in one go
    int idx;
    double m, dndm, growthf, lm1, lm2, dlm, lm, RhoM, RhoCr, dfcoll_dm, Halo, dHalo_dm, RhoM_grid, Ln10, HI, RhoDM, res;

    RhoCr = 2.775e11 * pow(cosmo_params->hlittle, 2); // Critical density, msun/Mpc^3
    RhoM = cosmo_params->OMm * RhoCr;
    RhoDM = ( cosmo_params->OMm - cosmo_params->OMb) * RhoCr;
    RhoM_grid = hmf_model==-1 ? RhoM * (1. + delta) : RhoM;
    Ln10 = log(10.);

    // mass range
    lm1 = -5.99;
    lm2 = 19.99;
    dlm = (lm2 - lm1)/((double) Boost_nmh - 1.0);

    growthf = dicke(z);
    fcoll = 0.0;
    Halo = 0.0;

    for (idx = 0; idx < Boost_nmh; idx++)
    {
        lm = lm1 + ((double)idx)*dlm;
        m = pow(10.0, lm);
        
        // 1 - hmf
        if (hmf_model == -1)
        {
            dndm = Conditional_HMF(growthf, m, delta, sigma2, cosmo_params);
        }
        else if (hmf_model == 0)
        {
            dndm = dNdM(growthf, m);
        }
        else if (hmf_model == 1)
        {
            dndm = dNdM_st(growthf, m);
        }
        else if (hmf_model == 2)
        {
            dndm = dNdM_WatsonFOF(growthf, m);
        }
        else if (hmf_model == 3)
        {
            dndm = dNdM_WatsonFOF_z(z, growthf, m);
        }
        else
        {
            LOG_ERROR("Wrong choice of hmf_model!");
            Throw(ValueError);
        }

        // 2 - fcoll
        dfcoll_dm = m * dndm / RhoM_grid;
        fcoll += m * dfcoll_dm * dlm * Ln10;

        // 3 - Halo
        HI = HaloProfile_Integrator_analytic(z, m); // \int dr r^2 \rho^2_{dm}
        dHalo_dm = 4. * PI * dndm * HI / (pow(RhoDM, 2.0) * pow(1. + z, 3.));
        Halo += m * dHalo_dm * dlm * Ln10;

    }
    
    res = pow(1. - fcoll, 2.) + Halo;
    res = fmax(1., res);
    
    return res;
}

/*
double Get_fcoll(double *m, double *dndm, double delta, int nm, int show_status)
{

    double OmM, h, RhoM, RhoM_Halo, dlnm, fun, result;
    int idx;

    OmM = 0.30964168161;
    h = 0.6766;
    RhoM = 2.775e11 * pow(h, 2.0) * OmM * (1 + delta);
    RhoM_Halo = 0.0;

    for (idx = 1; idx < nm; idx++)
    {
        dlnm = log(m[idx]) - log(m[idx - 1]);
        fun = dlnm * pow(m[idx], 2.0) * dndm[idx];
        RhoM_Halo = RhoM_Halo + fun;
    }

    result = RhoM_Halo / RhoM;
    if (show_status)
    {
        printf("RhoM_Halo = %E, RhoM = %E, fcoll = %E ", RhoM_Halo, RhoM, result);
    }

    return result;
}

double Boost_Factor_kernel(double z, double *m, double *dndm, double delta, int nm, int show_status)
{
    double fcoll, IGM, Halo, lnm[nm], fun[nm], mh, HI, h, OmC, RhoC, RhoC2, pi, Boost;
    int idx;
    // Takes about 1.5s for nr = 2000

    pi = 3.141592653589793;
    h = 0.6766;
    OmC = 0.260667;
    RhoC = 2.775e11 * pow(h, 2.0) * OmC; // comoving, Msun/Mpc^3
    RhoC2 = pow(RhoC, 2.0);

    // First get IGM part
    fcoll = Get_fcoll(m, dndm, delta, nm, show_status);
    IGM = pow(1 - fcoll, 2.0);

    // Now do halo
    for (idx = 0; idx < nm; idx++)
    {
        mh = m[idx];
        lnm[idx] = log(mh);
        HI = HaloProfile_Integrator(z, mh);
        fun[idx] = mh * dndm[idx] * HI;
    }

    Halo = Integrate(lnm, fun, nm);
    Halo = 4 * pi * Halo / (RhoC2 * pow(1.0 + z, 3.0));

    Boost = IGM + Halo;

    // might be some numerical stuff
    if (show_status)
    {
        printf("IGM = %E, Halo = %E, B = %E", IGM, Halo, Boost);
    }

    // if (Boost < 1.0)
    // {
    //     Boost = 1.0;
    // }
    return Boost;
}

int main()
{
    int nr = 1000, idx;
    double z, mh, m, r[nr], rho[nr], mass, dif;

    z = 3;
    mh = 1.0e12;

    m = HaloProfile(z, mh, r, rho, nr);
    dif = diff(m, mh);

    //for (idx = 0; idx < nr; idx++)
    //{
    //      printf("%E  %E\n", r[idx], rho[idx]);
    //}


    double B;
    int N = 100;
    for (idx = 0; idx < N; idx++)
    {
        B = Boost_Factor_Integrator(z, m_vec, dndm_vec, 0.0, 1000);
        printf("%d, B = %E\n", idx, B);
    }

    return 0;
}

// Now comes the part where we need 21cmFAST

double BoostFactor(double z, double growthf, double Delta, double sigma2, struct CosmoParams *cosmo_params, int hmf_model)
{
    double m, m1, m2, lnm, lnm1, lnm2, dlnm, m_vec[Boost_nmh], dndm_vec[Boost_nmh], result, hmf;
    int idx, show_status;

    // Min and Max for hmf, max is default interp table range, do not change

    m1 = 1.01e-6;
    // m1 = 1.01e0;

    m2 = 0.99e20;
    lnm1 = log(m1);
    lnm2 = log(m2);
    dlnm = (lnm2 - lnm1) / ((double)Boost_nmh - 1.0);
    lnm = lnm1;

    if (hmf_model != -1)
    {
        show_status = 1;
    }
    else
    {
        show_status = 0;
    }

    // Get a set of hmf
    for (idx = 0; idx < Boost_nmh; idx++)
    {
        m = exp(lnm);
        m_vec[idx] = m;
        if (hmf_model == -1)
        {
            hmf = Conditional_HMF(growthf, m, Delta, sigma2, cosmo_params);
        }
        else if (hmf_model == 0)
        {
            hmf = dNdM(growthf, m);
        }
        else if (hmf_model == 1)
        {
            hmf = dNdM_st(growthf, m);
        }
        else if (hmf_model == 2)
        {
            hmf = dNdM_WatsonFOF(growthf, m);
        }
        else if (hmf_model == 3)
        {
            hmf = dNdM_WatsonFOF_z(z, growthf, m);
        }
        else
        {
            hmf = 0.0;
            LOG_ERROR("Wrong choice of HMF!");
            Throw(ValueError);
        }
        dndm_vec[idx] = hmf;
        lnm = lnm + dlnm;
    }

    result = Boost_Factor_kernel(z, m_vec, dndm_vec, Delta, Boost_nmh, show_status);

    if (hmf_model != -1)
    {
        printf(", z = %E, B = %E\n", z, result);
    }
    return result;
}
*/

double DM_Term(double z, double Boost, double Delta, struct AstroParams *astro_params, struct FlagOptions *flag_options, double xe, double dt_dzp, int Type)
{
    // Ok I am doing everything in SI units
    double RhoC_c2_2, Pann, Q, P27, dEdVdt_Inj, nH, dxe_dt, dT_dt, fHe, dxe_dz, dT_dz, B, ntot, f1, f4, kB;
    if (flag_options->USE_HALO_BOOST)
    {
        B = Boost;
    }
    else
    {
        B = 1.0;
    }

    // RhoC_c2_2 = [OmC * Rho_cr *c^2]^2
    RhoC_c2_2 = 4.0610249075234377e-20;
    Q = 1.602176634E-19;
    kB = 1.38064852E-23;
    fHe = 0.08112582781456953;

    P27 = (double)astro_params->Pann27;
    Pann = P27 * 1.0E-27 * 1.0E-6 / (1.0E9 * Q);

    dEdVdt_Inj = B * Pann * RhoC_c2_2 * pow(1.0 + z, 6.0);

    nH = 0.19015670534605955 * pow(1.0 + z, 3.0) * (1.0 + Delta);

    // First do ionisation, I am not gonna do LyA yet
    f1 = (1.0 - xe) / 3.0;
    // ok this is somewhat inaccurate, we are not doing Helium here but we are assuming that xe is shared by H and He. should mention this in our paper
    dxe_dt = f1 * dEdVdt_Inj / (13.6 * Q * nH);

    // Do heating
    f4 = (1.0 + 2 * xe) / 3.0;
    ntot = nH * (1. + fHe + xe + 2.0 * xe * fHe); // proton, neutral He, e from H, e from He
    dT_dt = f4 * dEdVdt_Inj * 2.0 / (3.0 * kB * ntot);

    dxe_dz = dxe_dt * dt_dzp;
    dT_dz = dT_dt * dt_dzp;

    if (Type == 1)
    {
        return dxe_dz;
    }
    else
    {
        return dT_dz;
    }
}
