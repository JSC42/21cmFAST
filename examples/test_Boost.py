'''
Example for running 21cmFAST simulation for inhomogeneous DM boost factor, default reproduces the IHM simulation of 2312.17499
contains some of my own packages, comments are very welcome, feel free to contact if cannot run this
Implicit dependencies:
1 - PyLab.cosmo_tools
    PyLab is my own archive for some useful functions, contains HyRec and Switch_xe_format modules, required for calling my own HyRec and switching xe format
    url: https://github.com/Junsong-Cang/PyLab
    GitHub version might be out of date because I am having some troubles doing git push with new machine, latest version available upon request
2 - DarkSide branch of 21cmFAST
    Forked from 21cmFAST repository, docs for inputs can be found in src/py21cmfast/inputs.py
    url: https://github.com/Junsong-Cang/21cmFAST
3 - DarkSide branch of Dark_CosmoMC, forked from CosmoMC repository
    Module needed for setting xe, Tk initial conditions for 21cmFAST, we only need the HyRec module
    url: https://github.com/Junsong-Cang/Dark_CosmoMC
'''

FileName = 'Boost.h5'
L_X = 40.5
cache_loc = 'Your cache loc' # need this to store coevals for Tb and density, can give smoother PS than lightcones

redshift = 4.0
Z_HEAT_MAX = 60
LC_Quantities = ('brightness_temp','Ts_box','xH_box','Tk_box','Boost_box', 'density') # Boost_box is the new output box
GLB_Quantities = ('brightness_temp','Ts_box','xH_box','Tk_box', 'Boost_box', 'density')

# New params:
USE_HALO_BOOST = True # use boost
INHOMO_HALO_BOOST = True # use INHOMOGENEOUS boost
Pann27 = 1.0 # Set <sigma v>/mdm, unit: 10^{-27} cm^3/s/GeV

from cosmo_tools import *
import py21cmfast as p21c

user_params = p21c.UserParams(
  HII_DIM = 300,
  N_THREADS = 1,
  USE_RELATIVE_VELOCITIES = False,
  USE_INTERPOLATION_TABLES = True,
  FAST_FCOLL_TABLES = False,
  HMF = 0,
  POWER_SPECTRUM = 2,
  BOX_LEN = 500)

astro_params = p21c.AstroParams(
  Pann27 = Pann27,
  L_X = L_X,
  F_STAR10 = -1.3,
  ALPHA_STAR = 0.5,
  F_ESC10 = -1.0,
  ALPHA_ESC = -0.5,
  M_TURN = 8.7,
  t_STAR = 0.5,
  NU_X_THRESH = 500.0)
  
flag_options = p21c.FlagOptions(
  USE_MINI_HALOS = False,
  USE_MASS_DEPENDENT_ZETA = True,
  INHOMO_RECO = True,
  USE_TS_FLUCT = True,
  USE_HALO_BOOST = USE_HALO_BOOST,
  INHOMO_HALO_BOOST = INHOMO_HALO_BOOST,
)

# ---- Initialise ----
start_time = time.time()
# Setting Initial conditions

InitialCondition = HyRec(Pann = Pann27*1E-27, Use_SSCK = 1) # if this fails due to package dependencies, can also run HyRec externally and read from file
z = InitialCondition['z'][::-1]
xe = InitialCondition['xe'][::-1]
Tk = InitialCondition['Tk'][::-1]

XION_at_Z_HEAT_MAX = np.interp(x = Z_HEAT_MAX, xp = z, fp = xe)
XION_at_Z_HEAT_MAX = Switch_xe_format(xe = XION_at_Z_HEAT_MAX, format = 0) # HyRec (max 1.15) and 21cmFAST (max 1) has slightly different definitions for xe

TK_at_Z_HEAT_MAX = np.interp(x = Z_HEAT_MAX, xp = z, fp = Tk)

with p21c.global_params.use(Z_HEAT_MAX = Z_HEAT_MAX, XION_at_Z_HEAT_MAX = XION_at_Z_HEAT_MAX, TK_at_Z_HEAT_MAX = TK_at_Z_HEAT_MAX):
  lc = p21c.run_lightcone(
    redshift=redshift, 
    max_redshift=Z_HEAT_MAX,
    astro_params=astro_params, 
    flag_options=flag_options,
    user_params = user_params,
    lightcone_quantities=LC_Quantities,
    global_quantities=GLB_Quantities,
    direc = cache_loc,
    write = True # write coeval boxes for better PS eveluation
    )
end_time = time.time()
print('--------RunTime--------')
print(end_time - start_time)

try:
  os.remove(FileName)
except:
  pass

lc.save(FileName)
