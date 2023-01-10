%General Settings

% -- Constants
G_=6.6740831313131E-11; % Gravitational Constant
h_=6.626070040818181818E-34; % Planck Constant
hb_=h_/(2*pi); % Reduced Planck Constant
Q_=1.602176634E-19; % Electrical Charge Quanta
Elec_Mass_eV_=510998.9461;
kB=1.38064852E-23; % Boltzmann constant
C_=299792458; % Speed of light in m/s

% -- Cosmology
Omega_M=0.3111;
Omega_R=9.065340263942517E-5;
Omega_L=1-Omega_M-Omega_R;
As=2.105209331337507e-09;

Obh2=0.02242;
Och2=0.11933;

hubble_=0.6766;
H0_=hubble_*(3.240755744239557e-18);
%H_=@(z) H0_*sqrt(Omega_L+Omega_M*(1+z).^3+Omega_R*(1+z).^4);
H_=@(z) Hubble(z)
dtdz_=@(z) 1./((1+z).*H_(z));

dt_=@(z) dtdz_(z).*z/1000;
TIME_=@(x) integral(dtdz_,x,Inf);

% -- Numbers
MpC_=3.086E22; % in m
Msun = 1.98847E30; % Solar mass in kg

TIME_(0)
