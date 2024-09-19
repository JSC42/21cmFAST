// DM/PBH Deposition efficiency interpolation module

// include these files if compiled independently
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "EFF_Tables.h"
#include "Useful_Functions.h"

double Interp_EFF_DM_Annihilation(double Mdm, double z, int Dep_Channel, int Reaction_Channel)
{
  /* Interpolate deposition efficiencies for DM annihilation
  -- inputs --
  Mdm : DM mass in GeV
  z: redshfit
  Dep_Channel: Deposition Channel
          1 - HIon
          3 - LyA
          4 - Heat
  Reaction_Channel: Product that DM decay/annihilate into
          1 - Photon
          2 - Electron
          3 - Higgs
          4 - Muon
          5 - Tau
          6 - Q
          7 - CHarm
          8 - Bottom
          9 - Top
          10 - W
          11 - Z
          12 - Gluon
  */

  double Mdm_eV, x_in, EFF;
  Mdm_eV = Mdm * 1.0E9;
  // Validate inputs
  if ((Dep_Channel != 1) && (Dep_Channel != 3) && (Dep_Channel != 4))
  {
    printf("Error: wrong choice of deposition channels.\n");
    exit(1);
  }

  // ---- Gamma ----
  if (Reaction_Channel == 1)
  {
    x_in = Mdm_eV;
    if (Dep_Channel == 1)
    {
      EFF = Interp_2D(x_in, 1 + z, Kinetic_Energy_Axis, Redshift_Axis, EFF_Phot_HIon_Ann_HMG, Kinetic_Energy_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 3)
    {
      EFF = Interp_2D(x_in, 1 + z, Kinetic_Energy_Axis, Redshift_Axis, EFF_Phot_LyA_Ann_HMG, Kinetic_Energy_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 4)
    {
      EFF = Interp_2D(x_in, 1 + z, Kinetic_Energy_Axis, Redshift_Axis, EFF_Phot_Heat_Ann_HMG, Kinetic_Energy_Axis_Size, Redshift_Size, 1);
    }
  }
  // ---- Elec ----
  else if (Reaction_Channel == 2)
  {
    x_in = Mdm_eV - Electron_Mass_eV;
    if (x_in < 0.0)
    {
      printf("Error: DM mass too small to annihilate into electron+positron.\n");
      exit(1);
    }
    if (Dep_Channel == 1)
    {
      EFF = Interp_2D(x_in, 1 + z, Kinetic_Energy_Axis, Redshift_Axis, EFF_Elec_HIon_Ann_HMG, Kinetic_Energy_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 3)
    {
      EFF = Interp_2D(x_in, 1 + z, Kinetic_Energy_Axis, Redshift_Axis, EFF_Elec_LyA_Ann_HMG, Kinetic_Energy_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 4)
    {
      EFF = Interp_2D(x_in, 1 + z, Kinetic_Energy_Axis, Redshift_Axis, EFF_Elec_Heat_Ann_HMG, Kinetic_Energy_Axis_Size, Redshift_Size, 1);
    }
  }
  // ---- Higgs ----
  else if (Reaction_Channel == 3)
  {
    x_in = Mdm_eV;
    if (Dep_Channel == 1)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Annihilation, Redshift_Axis, EFF_Higgs_HIon_Ann_HMG, Mdm_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 3)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Annihilation, Redshift_Axis, EFF_Higgs_LyA_Ann_HMG, Mdm_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 4)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Annihilation, Redshift_Axis, EFF_Higgs_Heat_Ann_HMG, Mdm_Axis_Size, Redshift_Size, 1);
    }
  }
  // ---- Muon ----
  else if (Reaction_Channel == 4)
  {
    x_in = Mdm_eV;
    if (Dep_Channel == 1)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Annihilation, Redshift_Axis, EFF_Muon_HIon_Ann_HMG, Mdm_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 3)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Annihilation, Redshift_Axis, EFF_Muon_LyA_Ann_HMG, Mdm_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 4)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Annihilation, Redshift_Axis, EFF_Muon_Heat_Ann_HMG, Mdm_Axis_Size, Redshift_Size, 1);
    }
  }
  // ---- Tau ----
  else if (Reaction_Channel == 5)
  {
    x_in = Mdm_eV;
    if (Dep_Channel == 1)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Annihilation, Redshift_Axis, EFF_Tau_HIon_Ann_HMG, Mdm_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 3)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Annihilation, Redshift_Axis, EFF_Tau_LyA_Ann_HMG, Mdm_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 4)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Annihilation, Redshift_Axis, EFF_Tau_Heat_Ann_HMG, Mdm_Axis_Size, Redshift_Size, 1);
    }
  }
  // ---- Q ----
  else if (Reaction_Channel == 6)
  {
    x_in = Mdm_eV;
    if (Dep_Channel == 1)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Annihilation, Redshift_Axis, EFF_Q_HIon_Ann_HMG, Mdm_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 3)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Annihilation, Redshift_Axis, EFF_Q_LyA_Ann_HMG, Mdm_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 4)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Annihilation, Redshift_Axis, EFF_Q_Heat_Ann_HMG, Mdm_Axis_Size, Redshift_Size, 1);
    }
  }
  // ---- Charm ----
  else if (Reaction_Channel == 7)
  {
    x_in = Mdm_eV;
    if (Dep_Channel == 1)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Annihilation, Redshift_Axis, EFF_Charm_HIon_Ann_HMG, Mdm_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 3)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Annihilation, Redshift_Axis, EFF_Charm_LyA_Ann_HMG, Mdm_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 4)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Annihilation, Redshift_Axis, EFF_Charm_Heat_Ann_HMG, Mdm_Axis_Size, Redshift_Size, 1);
    }
  }
  // ---- Bottom ----
  else if (Reaction_Channel == 8)
  {
    x_in = Mdm_eV;
    if (Dep_Channel == 1)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Annihilation, Redshift_Axis, EFF_Bottom_HIon_Ann_HMG, Mdm_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 3)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Annihilation, Redshift_Axis, EFF_Bottom_LyA_Ann_HMG, Mdm_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 4)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Annihilation, Redshift_Axis, EFF_Bottom_Heat_Ann_HMG, Mdm_Axis_Size, Redshift_Size, 1);
    }
  }
  // ---- Top ----
  else if (Reaction_Channel == 9)
  {
    x_in = Mdm_eV;
    if (Dep_Channel == 1)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Annihilation, Redshift_Axis, EFF_Top_HIon_Ann_HMG, Mdm_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 3)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Annihilation, Redshift_Axis, EFF_Top_LyA_Ann_HMG, Mdm_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 4)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Annihilation, Redshift_Axis, EFF_Top_Heat_Ann_HMG, Mdm_Axis_Size, Redshift_Size, 1);
    }
  }
  // ---- W ----
  else if (Reaction_Channel == 10)
  {
    x_in = Mdm_eV;
    if (Dep_Channel == 1)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Annihilation, Redshift_Axis, EFF_W_HIon_Ann_HMG, Mdm_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 3)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Annihilation, Redshift_Axis, EFF_W_LyA_Ann_HMG, Mdm_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 4)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Annihilation, Redshift_Axis, EFF_W_Heat_Ann_HMG, Mdm_Axis_Size, Redshift_Size, 1);
    }
  }
  // ---- Z ----
  else if (Reaction_Channel == 11)
  {
    x_in = Mdm_eV;
    if (Dep_Channel == 1)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Annihilation, Redshift_Axis, EFF_Z_HIon_Ann_HMG, Mdm_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 3)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Annihilation, Redshift_Axis, EFF_Z_LyA_Ann_HMG, Mdm_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 4)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Annihilation, Redshift_Axis, EFF_Z_Heat_Ann_HMG, Mdm_Axis_Size, Redshift_Size, 1);
    }
  }
  // ---- Gluon ----
  else if (Reaction_Channel == 12)
  {
    x_in = Mdm_eV;
    if (Dep_Channel == 1)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Annihilation, Redshift_Axis, EFF_Gluon_HIon_Ann_HMG, Mdm_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 3)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Annihilation, Redshift_Axis, EFF_Gluon_LyA_Ann_HMG, Mdm_Axis_Size, Redshift_Size, 1);
    }
    else if (Dep_Channel == 4)
    {
      EFF = Interp_2D(x_in, 1 + z, Mdm_Axis_Annihilation, Redshift_Axis, EFF_Gluon_Heat_Ann_HMG, Mdm_Axis_Size, Redshift_Size, 1);
    }
  }
  else
  {
    printf("Unknown DM decay/annihilation channel.\n");
    exit(1);
  }
  return EFF;
}

/*
int main()
{
  int idx, nz, dep_channel, Particle;
  double z1, z2, dz, z, mdm, fc;
  
  mdm = 0.0010573715456674328;
  dep_channel = 1;
  Particle = 1;

  for (idx=0; idx < Redshift_Size; idx++)
  {
    z = Redshift_Axis[idx];
    fc = Interp_EFF_DM_Annihilation(mdm, z, dep_channel, Particle);
    printf("%f %f\n", z, fc);
    z +=dz;
  }
  
  return 0;
}
*/