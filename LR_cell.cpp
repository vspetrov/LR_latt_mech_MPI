#include "LR_cell.h"
#include "LR_lattice.h"
#include "math.h"
#include "parallel.h"

const double G_Na=23.;
const double E_Na=54.4;
double G_K=0.282; //0.705
const double E_Kl=-87.2268;
const double E_k=-77.;
double G_Kl;
double G_si=0.09;//0.06

/* Mechanics constants */
static const double G_s_mech = 0.005;
static const double E_mech_threshold = 20;
static const double V_t = 10;
static const double varepsilon_1 = 0.1;
static const double varepsilon_2 = 0.01;
static const double V_o = -80;
static const double S_o = 0;

double p,q;
double i_Na;
double i_si;
double E_si;
double i_K;
double Xi;
double i_Kl;
double Kl_inf;
double alpha_Kl, betta_Kl;
double i_Kp;
double Kp;
double i_b;
double alpha, betta;

double U, m, h, j, d, f, X, Cai, s_tension, tension;

inline double sigmoid(double x){
    return 1./(1+exp(-x*10));
}

double VFunction(int &ii, int &jj) //find V for the Cell numbered CN_c, CN_r
{
	G_Kl = gk1[ii*N+jj];	
	G_K = Gk[ii*N+jj];
	G_si = Gsi[ii*N+jj];
	U=V[ii*N+jj];
	m=Cell[ii*N+jj].m;
	h=Cell[ii*N+jj].h;
	j=Cell[ii*N+jj].j;
	d=Cell[ii*N+jj].d;
	f=Cell[ii*N+jj].f;
	X=Cell[ii*N+jj].X;
	Cai=Cell[ii*N+jj].Cai;
    s_tension = Cell[ii*N+jj].s_tension;


	//find fast sodium current I_Na
	i_Na=G_Na*m*m*m*h*j*(U-E_Na);
	
	//find slow-inward current I_si
	E_si=7.7-13.0287*log(Cai);
	i_si=G_si*d*f*(U-E_si);

	//find time-dependent potassium current I_K
	if (U>-100.)
		Xi=2.837*(exp(0.04*(U+77.))-1.)/((U+77.)*exp(0.04*(U+35.)));
	else
		Xi=1.;
	i_K=G_K*X*Xi*(U-E_k);
	
	//find time-independent potassium current I_K1
	alpha_Kl=1.02/(1+exp(0.2385*(U-E_Kl-59.215)));
	betta_Kl=(0.49124*exp(0.08032*(U-E_Kl+5.476))+
			exp(0.06175*(U-E_Kl-594.31)))/(1.+exp(-0.5143*(U-E_Kl+4.753)));
	Kl_inf=alpha_Kl/(alpha_Kl+betta_Kl);
	i_Kl=G_Kl*Kl_inf*(U-E_Kl);
		
	//find platou potassium current I_Kp
	Kp=1./(1.+exp((7.488-U)/5.98));
	i_Kp=0.0183*Kp*(U-E_Kl);
		
	//find background current I_b
	i_b=0.03921*(U+59.87);

    tension = S_o - s_tension;
    Cell[ii*N+jj].tension = tension;
    double I_SAC =  G_s_mech * sigmoid(tension)*tension*(U - E_mech_threshold);

    return -(i_Na+i_si+i_K+i_Kl+i_Kp+i_b+I_ext[ii*N+jj] + I_SAC);
}

double mFunction(double &delta_t)
{
	alpha=0.32*(U+47.13)/(1.-exp(-0.1*(U+47.13)));
	betta=0.08*exp(-U/11.);
	p=alpha/(alpha+betta);
	q=1./(alpha+betta);
	return p+(m-p)*exp(-delta_t/q);
}

double hFunction(double &delta_t)
{
	if (U < -40.) 
		{
			alpha=0.135*exp((80.+U)/-6.8);
			betta=3.56*exp(0.079*U)+3.1e5*exp(0.35*U);
		}
	else
		{
			alpha=0.;
			betta=1./(0.13*(1.+exp((U+10.66)/-11.1)));
		}
	p=alpha/(alpha+betta);
	q=1./(alpha+betta);
	return p+(h-p)*exp(-delta_t/q);
}

double jFunction(double &delta_t)
{
	if (U < -40.) 
		{
			alpha=(-1.2714e5*exp(0.2444*U)-3.474e-5*exp(-0.04391*U))*
			(U+37.78)/(1+exp(0.311*(U+79.23)));
			betta=0.1212*exp(-0.01052*U)/(1.+exp(-0.1378*(U+40.14)));
		}
	else
		{
			alpha=0.;
			betta=0.3*exp(-2.535e-7*U)/(1.+exp(-0.1*(U+32.)));
		}
	p=alpha/(alpha+betta);
	q=1./(alpha+betta);
	return p+(j-p)*exp(-delta_t/q);
}

double dFunction(double &delta_t)
{
	alpha=0.095*exp(-0.01*(U-5.0))/(1.+exp(-0.072*(U-5.0)));
	betta=0.07*exp(-0.017*(U+44.0))/(1.+exp(0.05*(U+44.0)));
	p=alpha/(alpha+betta);
	q=1./(alpha+betta);
	return p+(d-p)*exp(-delta_t/q);
}

double fFunction(double &delta_t)
{
	alpha=0.012*exp(-0.008*(U+28.))/(1.+exp(0.15*(U+28.)));
	betta=0.0065*exp(-0.02*(U+30.))/(1.+exp(-0.2*(U+30.)));
	p=alpha/(alpha+betta);
	q=1./(alpha+betta);
	return p+(f-p)*exp(-delta_t/q);
}

double XFunction(double &delta_t)
{
	alpha=0.0005*exp(0.083*(U+50.))/(1.+exp(0.057*(U+50.)));
	betta=0.0013*exp(-0.06*(U+20.))/(1.+exp(-0.04*(U+20.)));
	p=alpha/(alpha+betta);
	q=1./(alpha+betta);
	return p+(X-p)*exp(-delta_t/q);
}

double CaiFunction()
{
	return -1e-4*i_si+0.07*(1e-4-Cai);
}

double s_tensionFunction(){
    double epsilon = varepsilon_1 + (varepsilon_2 - varepsilon_1)*sigmoid(U-V_o);
    return epsilon*(U/V_t-s_tension);
}
