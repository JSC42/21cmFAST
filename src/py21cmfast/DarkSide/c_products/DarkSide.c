// A down-sized test version of DarkSide.h, used here only for testing interpolation
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define print_debug_info 0
#define Boost_nmh 1000
#define Boost_Interp_Table_Size 200
#include "Tables.h"

double Read_EFF_Tab(int idxz, int idxE, double *Tab, int UseLog)
{    /* 
    A way to read 2D EFF table, I can't get pointer to work
    ---- inputs ----
    idxz : z index
    idxE : Ek index
    Tab : EFF Table
    UseLog : return log of fc, useful if u want to do subsequent interpolation in log space
    */
    int id;
    double r;
    if ((idxE > Ek_axis_Size - 1) || (idxz > Redshift_Size - 1))
    {
        printf("Error: array index overflow.\n");
        exit(1);
    }
    id = idxE * Redshift_Size + idxz;
    r = Tab[id];
    if (r < 0)
    {
        printf("Error: fc is negative.\n");
        exit(1);
    }
    if (r < 1E-100)
    {
        r=1E-100;
    }
    if (UseLog)
    {
        r=log(r);
    }
    return r;
}

int Find_Sign(double x)
{
    // Find sign of a number
    int result;
    if (x > 0.0)
    {
        result = 1;
    }
    else if (x < 0.0)
    {
        result = -1;
    }
    else
    {
        // maybe x is too close to 0?
        result = 0;
    }
    return result;
}

int Find_Index_1D(double *Tab, double x, int nx)
{
    /* 
    From table Tab with length nx, find closest left element index for x
    return -1 if not in range
    */
    double x1, x2, x3;
    int id1, id2, id3, Stop, s1, s2, s3, idx, count;
    Stop = 0;
    id1 = 0;
    id3 = nx - 1;

    count = 0;
    while (Stop == 0)
    {
        count = count + 1;
        id2 = (int)round((((double)id1) + ((double)id3)) / 2.0);
        if (id3 == id1 + 1)
        {
            idx = id1;
            Stop = 1;
        }

        x1 = Tab[id1];
        x2 = Tab[id2];
        x3 = Tab[id3];

        s1 = Find_Sign(x - x1);
        s2 = Find_Sign(x - x2);
        s3 = Find_Sign(x - x3);
        if (s1 != s2)
        {
            id3 = id2;
        }
        else if (s2 != s3)
        {
            id1 = id2;
        }
        if ((s1 == s3) || (count > 50))
        {
            idx = -1;
            Stop = 1;
        }
    }
    return idx;
}

double Interp_2D(double Ek, double z, double *Tab)
{
    /* 
    Interpolate (in log) from a 2D table
    -- inputs --
        Ek : kinetic energy
        z : redshift
        Tab : EFF table
     */
    int eid1, eid2, zid1, zid2;
    double le, le1, le2, lz, lz1, lz2, f11, f12, f21, f22, f1, f2, F1, F2, f;
    eid1 = Find_Index_1D(Ek_GeV_axis, Ek, Ek_axis_Size);
    zid1 = Find_Index_1D(Redshift_axis, z, Redshift_Size);
    if (zid1 == -1)
    {
        return 0.0;
    }
    if (eid1 == -1)
    {
        printf("Error: Ek not in range, exitting, debugging info:.\n");
        printf("Ek1 = %E, Ek = %E, Ek2 = %E\n", Ek_GeV_axis[eid1], Ek, Ek_GeV_axis[eid1 + 1]);
        printf("min(Ek) = %E, Ek = %E, max(Ek) = %E\n", Ek_GeV_axis[0], Ek, Ek_GeV_axis[Ek_axis_Size - 1]);
        exit(1);
    }

    eid2 = eid1 + 1;
    zid2 = zid1 + 1;

    le = log10(Ek);
    le1 = log10(Ek_GeV_axis[eid1]);
    le2 = log10(Ek_GeV_axis[eid2]);
    lz = log10(1+z);
    lz1 = log10(1+Redshift_axis[zid1]);
    lz2 = log10(1+Redshift_axis[zid2]);

    f11 = Read_EFF_Tab(zid1, eid1, Tab, 1);
    f12 = Read_EFF_Tab(zid1, eid2, Tab, 1);
    f21 = Read_EFF_Tab(zid2, eid1, Tab, 1);
    f22 = Read_EFF_Tab(zid2, eid2, Tab, 1);

    // fix E1
    f1 = f11;
    f2 = f21;
    F1 = (f2 - f1) / (lz2 - lz1) * (lz - lz1) + f1;
    // fix E2
    f1 = f12;
    f2 = f22;
    F2 = (f2 - f1) / (lz2 - lz1) * (lz - lz1) + f1;

    f = (F2 - F1) / (le2 - le1) * (le - le1) + F1;
    f = exp(f);
    return f;
}

double Interp_EFF(double z, double m, int particle, int channel, int HMF, int PS, int UseBoost)
{
    /*
    Interpolate DM deposition efficiencies
    -- inputs --
    z : z
    m : DM mass in GeV
    particle : particle that DM annihilate into
        0 - gamma gamma
        1 - electron + positron
    channel : deposition channel
        0 - HIon
        1 - LyA
        2 - Heat
    HMF : HMF model, follows UserParams
    PS : POWER_SPECTRUM in UserParams
    UseBoost : Whether to use Halo Boost
        0 - no, ignore HMF and PS
        1 - yes
    */
    double Ek, r;
    double EFF_Tab[Redshift_Size*Ek_axis_Size];
    CopyEFF(EFF_Tab, particle, channel, HMF, PS, UseBoost);
    Ek = (particle==1) ? m-0.5109989E-3 : m;
    if ((HMF < 0) || (HMF > 3))
    {
        printf("Error: unsupported HMF, HMF=%d\n", HMF);
        exit(1);
    }
    if ((PS < 0) || (PS > 4))
    {
        printf("Error: unsupported PS, PS=%d\n", PS);
        exit(1);
    }
    if (Ek < 0.0)
    {
        printf("Error: Ek is negative, Ek=%E\n", Ek);
        exit(1);
    }
    r = Interp_2D(Ek, z, EFF_Tab);
    return r;
}
