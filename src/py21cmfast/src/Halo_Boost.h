// Halo Boost Factor module, can be run and testes outside 21cmFAST

// #include <stdio.h>
// #include <stdlib.h>
// #include <math.h>
#define HaloProfile_Length 2000
#define Boost_nmh 1000

double HaloProfile_Kernel(double z, double mh, double r, int type)
{
    double h, OmM, OmL, OmB, Grav, pc, zp, pi, c, Omega_m_z, d, Delta_c, F, r_vir, p1, p2, p3, RhoM, H, rho_c, delta_c, x;
    h = 0.6766;
    OmM = 0.30964168161;
    OmB = 0.04897468161;
    OmL = 1.0 - OmM;
    Grav = 6.67259e-8;                     //  cm^3/g/s^2
    pc = 3.26 * 3e10 * 365 * 24 * 60 * 60; // cm
    zp = 1.0 + z;
    pi = 3.141592653589793;

    c = pow(10, (1.071 - 0.098 * (log10(mh) - 12))) / zp; // Concentration

    Omega_m_z = OmM * pow(zp, 3.0) / (OmM * pow(zp, 3.0) + OmL);

    d = Omega_m_z - 1;
    Delta_c = 18.0 * pow(pi, 2.0) + 82.0 * d - 39.0 * pow(d, 2.0);

    F = log(1 + c) - c / (1 + c);

    // Get r_vir
    p1 = 0.784 * pow((mh / 1e8 * h), 1.0 / 3.0);
    p2 = pow((OmM / Omega_m_z * Delta_c / 18 / (pi * pi)), -1.0 / 3.0);
    p3 = (10.0 / zp / h) * 1.e3;
    r_vir = p1 * p2 * p3;

    if (r > r_vir)
    {
        RhoM = 0;
    }
    else
    {
        H = h * 100.0 * sqrt(OmL + OmM * pow(zp, 3.0));
        rho_c = 3.0 * pow(H * 1.0e5 / 1.0e6 / pc, 2.0) / (8 * pi * Grav);
        rho_c = rho_c / 1.989e33 * pow(pc * 1.0e6, 3.0);
        delta_c = Delta_c / 3 * pow(c, 3.0) / F;
        x = r / r_vir;
        RhoM = rho_c * delta_c / (c * x * pow(1 + c * x, 2.0));
        RhoM = RhoM / 1.0e18;
    }
    if (type == 0)
    {
        return RhoM;
    }
    else
    {
        return r_vir;
    }
}

// void HaloProfile_vec(double z, double mh, double *r, double *rho, int LenR)

double diff(double x1, double x2)
{
    double d1, d2, r;
    d1 = (x2 - x1) / x1;
    d2 = (x2 - x1) / x2;
    if (d2 > d1)
    {
        r = d2;
    }
    else
    {
        r = d1;
    }

    return fabs(r);
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

    r_max = HaloProfile_Kernel(z, mh, 10.0, 1);
    logspace_vec(r_min_ratio * r_max, r_max, nr, r);

    for (idx = 0; idx < nr; idx++)
    {
        r_ = r[idx];
        rho[idx] = HaloProfile_Kernel(z, mh, r_, 0);

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

/*
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
*/

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
    ntot = nH * (1. + fHe + xe + 2.0*xe*fHe); // proton, neutral He, e from H, e from He
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
