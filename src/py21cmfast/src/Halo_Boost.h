// Halo Boost Factor module
#define print_debug_info 0
#define Boost_nmh 1000
#define Boost_Interp_Table_Size 200

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

double HaloProfile_Integrator(double z, double mh)
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
    res = res * 1.0e18; // convert to msun^2/Mpc^3

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
    // Allow some redundencies
    if (id1 < 0)
    {
        if (id1 == -1)
        {
            id1 = 0;
        }
        else
        {
            printf("crash imminent @ Interp_Fast, overflow detected, xmin = %E, xmax = %E, x = %E, id1 = %d\n", xmin, xmax, x, id1);
            LOG_ERROR("Interpolation overflow");
            Throw(ValueError);
        }
    }
    if (id1 > nx - 2)
    {
        if (id1 <= nx)
        {
            id1 = nx - 2;
        }
        else
        {
            printf("crash imminent @ Interp_Fast, overflow detected, xmin = %E, xmax = %E, x = %E, id1 = %d\n", xmin, xmax, x, id1);
            LOG_ERROR("Interpolation overflow");
            Throw(ValueError);
        }
    }
    id2 = id1 + 1;
    x1 = xmin + ((double)id1) * dx;
    x2 = x1 + dx;
    y1 = Tab[id1];
    y2 = Tab[id2];
    y = (y2 - y1) * (x - x1) / (x2 - x1) + y1;

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
    // printf("Using debug version for Conditional_HMF\n");
    // return fabs(-RhoM * dNdM_double / M / sqrt(2 * PI));

}

double BoostFactor(double z, double growthf, double delta, double sigma2, struct CosmoParams *cosmo_params, int hmf_model)
{
    // Get the boost factor all in one go, set growthf externaly to avoid repeated computations
    int idx;
    double m, dndm, lm1, lm2, dlm, lm, RhoM, RhoCr, dfcoll_dm, Halo, dHalo_dm, RhoM_grid, Ln10, HI, RhoDM, res, fcoll;

    RhoCr = 2.775e11 * pow(cosmo_params->hlittle, 2); // Critical density, msun/Mpc^3
    RhoM = cosmo_params->OMm * RhoCr;
    RhoDM = (cosmo_params->OMm - cosmo_params->OMb) * RhoCr;
    RhoM_grid = hmf_model == -1 ? RhoM * (1. + delta) : RhoM;
    Ln10 = log(10.);

    // mass range
    lm1 = -5.99;
    lm2 = 19.99;
    dlm = (lm2 - lm1) / ((double)Boost_nmh - 1.0);

    fcoll = 0.0;
    Halo = 0.0;

    for (idx = 0; idx < Boost_nmh; idx++)
    {
        lm = lm1 + ((double)idx) * dlm;
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
        HI = HaloProfile_Integrator(z, m); // \int dr r^2 \rho^2_{dm}
        dHalo_dm = 4. * PI * dndm * HI / (pow(RhoDM, 2.0) * pow(1. + z, 3.));
        Halo += m * dHalo_dm * dlm * Ln10;
    }

    res = pow(1. - fcoll, 2.) + Halo;
    
    res = fmax(1., res);
    /*
    if (print_debug_info && hmf_model != -1)
    {
        FILE *OutputFile;
        OutputFile = fopen("Fcoll_tmp.txt", "a");
        fprintf(OutputFile, "%f   %E\n", z, fcoll);
        fclose(OutputFile);
    }
    */

    return res;
}

double EoR_Drive_DM(double z, double Boost, double Delta, struct AstroParams *astro_params, struct FlagOptions *flag_options, double xe, double dt_dzp, int Type)
{
    /*
    Get dxe/dz and dT/dz from DM injection
    I am doing everything in SI units
    */
    double RhoC_c2_2, Pann, Q, P27, dEdVdt_Inj, nH, dxe_dt, dT_dt, fHe, dxe_dz, dT_dz, B, ntot, kB, f_HIon, f_Heat, f_LyA;
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
    f_HIon = (1.0 - xe) / 3.0;
    // ok this is somewhat inaccurate, we are not doing Helium here but we are assuming that xe is shared by H and He. should mention this in our paper
    dxe_dt = f_HIon * dEdVdt_Inj / (13.6 * Q * nH);

    // Do heating
    f_Heat = (1.0 + 2 * xe) / 3.0;
    ntot = nH * (1. + fHe + xe + 2.0 * xe * fHe); // proton, neutral He, e from H, e from He
    dT_dt = f_Heat * dEdVdt_Inj * 2.0 / (3.0 * kB * ntot);

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
