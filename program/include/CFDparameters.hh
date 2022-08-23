#ifndef CFDparameters_h
#define CFDparameters_h 1


class CFDparameters{
public:
  CFDparameters(unsigned short a, unsigned short b,unsigned short c, unsigned short d, double l,double err_l, double h,double err_h, double delta1,double err_delta1, double delta2, double r, double err_r){
sed_delay = a;

sed_frac = b;
si_delay = c;
si_frac = d;
sigma_LF = l;
err_sigma_LF = err_l;
sigma_HF = h;
err_sigma_HF = err_h;
delta_tof = delta1;
err_delta_tof = err_delta1;
delta_tof_mt = delta2;
 integral_ratio = r;
 err_integral_ratio = err_r;
};
public:
unsigned short sed_delay;
unsigned short sed_frac;
unsigned short si_delay;
unsigned short si_frac;
double sigma_LF;
  double err_sigma_LF;
double sigma_HF;
  double err_sigma_HF;
double delta_tof;
  double err_delta_tof;
double delta_tof_mt;//measured - evaluated//SRIM corrected
  double integral_ratio;
  double err_integral_ratio;
};
#endif
