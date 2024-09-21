#include <stdio.h>
#include <math.h>

double PeeblesFactor(double z, double xe, double Tk)
{
    /*
    Compute Peebles C factor following Dodelson's Modern Cosmology II (see p108)
    */
    double E0, alpha_fsc, hbar, LightSpeed, me, np0, x, fx1, fx2, fx, beta2, pi, L2y, nH, H, La1, La2, La, result;
    alpha_fsc = 1/137.03599976;
    E0 = 13.6 * 1.602176634E-19;
    hbar = 1.0545718002693302e-34;
    LightSpeed = 299792458.0;
    me = 9.109383632044565e-31;
    np0 = 0.1901567053460595;// Number density of H nuclei today
    pi = 3.14159265358979323846264338327;
    
    x = E0/(Tk * 1.38064852E-23);
    fx1 = log(x);
    fx2 = x * exp(x/4.0);
    fx = fx1/fx2;
    beta2 = 9.78 * pow(E0/(2 * pi), 1.5) * pow(alpha_fsc, 2.0) /(pow(me, 0.5) * LightSpeed * hbar) * fx;
    L2y = 8.227;
    nH = np0 * pow(1.0+z, 3.0) * (1.0-xe);
    H = 2.192695336552484e-18 * pow(0.69026731839 + 0.30964168161*pow(1.0+z, 3.0) + 9.1E-5 * pow(1.0+z, 4.0), 0.5);
    // H = hubble(z);
    La1 = H * pow(3.0*E0, 3.0);
    La2 = nH * pow(8*pi, 2.0) * pow(hbar*LightSpeed, 3.0);
    La = La1/La2;
    result = (La + L2y)/(La + L2y + beta2);
    return result;
}
