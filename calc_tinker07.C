#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>

#include "proto.h"
#include "allvars.h"
#include "nrutil.h"
#include <gsl/gsl_sf.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multiroots.h>

using namespace std;

static char TransTable[2000];
static double xp[1000],yp[1000],yp2[1000];
static double tab_z[NPOINTS],tab_chi[NPOINTS],err_chi[NPOINTS];
static double scale_f[NPOINTS],GF[NPOINTS],GF2[NPOINTS];

static double tab_m[NPOINTS];
static double tab_z_mf[NPOINTS];
static double tab_mf_norm[NPOINTS],err_mf_norm[NPOINTS];
static double **tab_mf,**err_mf;
static double **tab_bh,**err_bh;
static double **tab_cv,**err_cv;
static double **tab_nuM, **err_nuM;

static double tab_R[NPOINTS],tab_dsdr[NPOINTS],tab_ds2dr2[NPOINTS],err_dsdr[NPOINTS],err_ds2dr2[NPOINTS];
static double tab_sig2[NPOINTS],err_sig2[NPOINTS];

static double tab_SC_dens[1000], tab_SC_velo[1000], err_SC_velo[1000];

static int WhichSpectrum, np, WhichWindow, OPT, OPT_fit, WhichGF;

static double bling;

static double r_tophat,Rsmooth,Delta_c,fdelta_c,Const_MF,st_norm;

static double AA,BB,CC;
static double B1,B2,B3,B4,B5;
static double nu, sigma, Omega_z;

static double UnitTime_in_s, UnitLength_in_cm, UnitMass_in_g, UnitVelocity_in_cm_per_s;
static double G, Hubble;

static double Norm, InitTime, calc_z, a_exp;
static double Dplus; /* growth factor */

static double Mh_min, Mh_max;

//HOD
struct HODparams{
  double logMmin;
  double sig_sq;
  double logM1;
  double alpha;
  double kappa;  
  double alpha_inc;
  double logMinc;
};
static struct HODparams hod;

static double thre;

struct get_mvir{
  double delta;
  double rd;
  double redshift;
};

static double zhalo, sigma_NL_mass;

static double tab_mod_ln_s2[NPOINTS];
static double tab_mod_ln_N[NPOINTS], tab_mod_ln_A[NPOINTS], tab_mod_ln_w2[NPOINTS];
static double err_mod_ln_N[NPOINTS], err_mod_ln_A[NPOINTS], err_mod_ln_w2[NPOINTS];

int main(int argc, char **argv){
  FILE *fout;
  int i,j;
	
  //input parameters======================================
  WhichSpectrum = 3;  // 1 Eisenstein-Hu, 2 Bond-Efstathiou, 3 CMBFAST table, 4 BBKS
  WhichGF = 2;        // 1 apploximation formula, 2 caluculate linear density growth eq.	
  OPT = 6;            // 0 JK  1 ST  2 PS 3 Yahagi 4 Reed et al (astro-ph/0607150v4) 5 Bhattacharya et al (astro-ph/1005.2239) 6 Tinker et al 2008, 2010
  OPT_fit = 2; // 1 smith et al 03, 2 takahashi et al 12

  if(argc!=8){
    printf( "usage: %s output Pk_lin_at_z=0 Om h0 ns HOD_param z\n", argv[0]);
    return 1;
  }
		
  sprintf(TransTable, "%s", argv[2]);
	
  //set cosmological model
  w  = -1.;	
  Omega = atof(argv[3]);
  OmegaLambda = 1.-Omega;
  HubbleParam = atof(argv[4]);
  ns = atof(argv[5]);
	
  FILE *fin;
  fin = fopen(argv[6], "r");
  if(fin == NULL){
    fprintf(stderr, "can not find %s.\n", argv[4]);
    exit(1);
  }
	
  fscanf(fin, "%lf %lf %lf %lf %lf\n",
	 &hod.logMmin, &hod.sig_sq, &hod.logM1, &hod.alpha, &hod.kappa);

  hod.logM1 = log10(hod.logM1);

  fclose(fin);
	
  Mh_min = 1e11;
  Mh_max = 1e16;
		
  zhalo = atof(argv[7]);
	
  if(WhichGF==2){growth_spline();}
	
  i=set_units();
  i=initialize_powerspectrum(WhichSpectrum);	
	
  stack_table_and_spline();	
  stack_data_and_spline_Pk();
		
  double rho_crit = 2.7754e11;
  double rhom_comv = rho_crit*Omega;

  //save_SC_param();

  char oname[256];
  sprintf(oname, "%s", argv[1]);

  fout = fopen(oname, "w");
  if(fout == NULL){
    printf("you can not make %s\n", oname);
    exit(1);
  }
  
  int nint  = 30;

  double vmin = -2500.0;
  double vmax = +2500.0;
  double tab_v[5*NPOINTS], tab_pdf[5*NPOINTS], err_pdf[5*NPOINTS];

  for(double r = 3.0; r <= 40.0; r += 3.0){

    sigma_NL_mass = TopHatSigma2_NL(r);
    
    // halo exclusion effect, see http://articles.adsabs.harvard.edu/pdf/2007MNRAS.374..477T
    double Rlim1 = r - r_delta(zhalo, Mh_min, 200.0) * (1.0+zhalo); // comoving
    double Mlim1 = 4.*M_PI*rhom_comv*200.0*pow(Rlim1, 3)/3.0;

    if(Mlim1 >= Mh_max) Mlim1 = Mh_max;
    
    double norm = 0;
    for(i=0;i<nint;i++){

      double xmin1 = log10(Mh_min);
      double xmax1 = log10(Mlim1);
      double dx1 = (xmax1-xmin1)/nint;
      double x1 = xmin1 + (i+0.5)*dx1;
      double mf1 = dndlogm_fast(x1, zhalo);
      
      for(j=0;j<nint;j++){
	
	double xmin2 = log10(Mh_min);
	double Rlim2 = r - r_delta(zhalo, pow(10., x1), 200.0) * (1.0+zhalo); // comoving 
	double Mlim2 = 4.*M_PI*rhom_comv*200.0*pow(Rlim2, 3)/3.0;
	if(Mlim2 >= Mh_max) Mlim2 = Mh_max;
	double xmax2 = log10(Mlim2);
	double dx2 = (xmax2-xmin2)/nint;
	double x2 = xmin2 +(j+0.5)*dx2;
	double mf2 = dndlogm_fast(x2, zhalo);

	double mass1 = pow(10., x1);
	double mass2 = pow(10., x2);

	double hod1 = HOD_gal(mass1, 0) + HOD_gal(mass1, 1);
	double hod2 = HOD_gal(mass2, 0) + HOD_gal(mass2, 1);

	norm += dx1*dx2*mf1*mf2*hod1*hod2;
		
      }
    }
    
    double mean_infall = 0.0;
    double sigr = 0.0;
    
    for(int iv = 0; iv<5*NPOINTS; iv++){
      tab_v[iv] = vmin + (vmax-vmin)*(iv)/(5*NPOINTS-1);
      tab_pdf[iv] = 0.0;
    }

    for(i=0;i<nint;i++){

      double xmin1 = log10(Mh_min);
      double xmax1 = log10(Mlim1);
      double dx1 = (xmax1-xmin1)/nint;
      double x1 = xmin1 + (i+0.5)*dx1;
      double mass1 = pow(10., x1);
      double rd1 = r_delta(zhalo, mass1, 200.0) * (1.0+zhalo);

      double mf1 = dndlogm_fast(log10(mass1), zhalo);
	
      fprintf(stdout, "#"); fflush(stdout);

      for(j=0;j<nint;j++){
	
	double xmin2 = log10(Mh_min);
	double Rlim2 = r - r_delta(zhalo, pow(10., x1), 200.0) * (1.0+zhalo); // comoving 
	double Mlim2 = 4.*M_PI*rhom_comv*200.0*pow(Rlim2, 3)/3.0;
	if(Mlim2 >= Mh_max) Mlim2 = Mh_max;

	double xmax2 = log10(Mlim2);
	double dx2 = (xmax2-xmin2)/nint;
	double x2 = xmin2 +(j+0.5)*dx2;
	double mass2 = pow(10., x2);
	
	double hod1 = HOD_gal(mass1, 0) + HOD_gal(mass1, 1);
	double hod2 = HOD_gal(mass2, 0) + HOD_gal(mass2, 1);
	
	double mf2 = dndlogm_fast(log10(mass2), zhalo);

	for(int iv = 0; iv<5*NPOINTS; iv++){

	  double pdf = PDF_halo(tab_v[iv], r, mass1, mass2);	  
	  if(isnan(pdf)){
	    printf("Nan found at v=%e r=%e\n", tab_v[iv], r);	    
	    exit(1);
	  }	  
	  tab_pdf[iv] += dx1*dx2*mf1*mf2*hod1*hod2*pdf;	  
	}

      }
    }

    for(int iv=0;iv<5*NPOINTS;iv++){
      tab_pdf[iv] = log10(tab_pdf[iv] *2 / (norm));
    }

    double yp1 = 1.e31;
    double ypn = 1.e31;
    spline(tab_v-1, tab_pdf-1, 5*NPOINTS, yp1, ypn, err_pdf-1);

    int Nv = 1000;
    double dv = (vmax-vmin)/Nv;
    double res0, res1, res2;
    
    res0=0;    
    for(int iv=0; iv<Nv/2; iv++){
      double vl = vmin + (2*iv+0)*dv;
      double vm = vmin + (2*iv+1)*dv;
      double vh = vmin + (2*iv+2)*dv;
	  
      double yl, ym, yh;
      splint(tab_v-1, tab_pdf-1, err_pdf-1, 5*NPOINTS, vl, &yl);
      splint(tab_v-1, tab_pdf-1, err_pdf-1, 5*NPOINTS, vm, &ym);
      splint(tab_v-1, tab_pdf-1, err_pdf-1, 5*NPOINTS, vh, &yh);

      yl = pow(10.0, yl);
      ym = pow(10.0, ym);
      yh = pow(10.0, yh);

      res0 += dv * (yl+4*ym+yh) / 3;
    }    

    res1=0;
    for(int iv=0; iv<Nv/2; iv++){
      double vl = vmin + (2*iv+0)*dv;
      double vm = vmin + (2*iv+1)*dv;
      double vh = vmin + (2*iv+2)*dv;
	  
      double yl, ym, yh;
      splint(tab_v-1, tab_pdf-1, err_pdf-1, 5*NPOINTS, vl, &yl);
      splint(tab_v-1, tab_pdf-1, err_pdf-1, 5*NPOINTS, vm, &ym);
      splint(tab_v-1, tab_pdf-1, err_pdf-1, 5*NPOINTS, vh, &yh);

      yl = pow(10.0, yl) * vl / res0;
      ym = pow(10.0, ym) * vm / res0;
      yh = pow(10.0, yh) * vh / res0;

      res1 += dv * (yl+4*ym+yh) / 3;
    }

    mean_infall = res1;

    res2=0;
    for(int iv=0; iv<Nv/2; iv++){
      double vl = vmin + (2*iv+0)*dv;
      double vm = vmin + (2*iv+1)*dv;
      double vh = vmin + (2*iv+2)*dv;
	  
      double yl, ym, yh;
      splint(tab_v-1, tab_pdf-1, err_pdf-1, 5*NPOINTS, vl, &yl);
      splint(tab_v-1, tab_pdf-1, err_pdf-1, 5*NPOINTS, vm, &ym);
      splint(tab_v-1, tab_pdf-1, err_pdf-1, 5*NPOINTS, vh, &yh);

      yl = pow(10.0, yl) * pow(vl-res1, 2) / res0;
      ym = pow(10.0, ym) * pow(vm-res1, 2) / res0;
      yh = pow(10.0, yh) * pow(vh-res1, 2) / res0;

      res2 += dv * (yl+4*ym+yh) / 3;
    }

    sigr = res2;
    
    fprintf(fout, "%e %e %e\n", r, mean_infall, sigr);
    fprintf(stdout, "%e %e %e\n", r, mean_infall, sigr);    
    

    /*
      char fname_each[256];
      sprintf(fname_each, "%s_R%3.2f_2halo_pdf.dat", argv[1], r);

      FILE *fout2;
      fout2 = fopen(fname_each, "w");
      if(fout2==NULL){
      fprintf(stderr, "cannot make %s.\n", fout2);
      exit(1);
      }

      for(int iv=0; iv<Nv; iv++){
      double v = vmin+(iv+0.5)*dv;
      double y;
      splint(tab_v-1, tab_pdf-1, err_pdf-1, NPOINTS, v, &y);
      fprintf(fout2, "%e %e\n", v, pow(10., y)/res0);
      }
      fclose(fout2);
      exit(1);
    */

  }
    
  fclose(fout);

  free_dmatrix(tab_mf,1,NPOINTS,1,NPOINTS);
  free_dmatrix(err_mf,1,NPOINTS,1,NPOINTS);
  free_dmatrix(tab_bh,1,NPOINTS,1,NPOINTS);
  free_dmatrix(err_bh,1,NPOINTS,1,NPOINTS);
  free_dmatrix(tab_cv,1,NPOINTS,1,NPOINTS);
  free_dmatrix(err_cv,1,NPOINTS,1,NPOINTS);
  free_dmatrix(tab_nuM,1,NPOINTS,1,NPOINTS);
  free_dmatrix(err_nuM,1,NPOINTS,1,NPOINTS);

  return 0;
}

 //==================================================
int initialize_powerspectrum(int Spectrum)
{
  double res;
  int i;

  for(i=0;i<1000;i++)
    xp[i]=yp[i]=yp2[i]=0;

  if(WhichSpectrum==3){
    //fprintf(stdout,"initialising...\n");
    readCMB_and_do_spline();
  }

  a_exp=1/(1+calc_z);

  AA=6.4/Gamma; 
  BB=3.0/Gamma; 
  CC=1.7/Gamma;  
  nu=1.13;


  B1=2.34/Gamma; 
  B2=3.89/Gamma; 
  B3=16.1/Gamma; 
  B4=5.46/Gamma; 
  B5=6.71/Gamma; 


  Norm = 1.0;
  Sigma8 = sqrt(TopHatSigma2(8.));
  //Norm=Sigma8*Sigma8/res;
  fprintf(stdout,"Sigma8 = %g \n",Sigma8);

  Dplus= GrowthFactor(a_exp, 1.0);

  return i;
}

 //==================================================
int set_units()    /* ... set some units */
{
  UnitLength_in_cm= 3.085678e21; /* 1.0 kpc */
  UnitMass_in_g=    1.989e43;    /* 1.0e10 solar masses */ 
  UnitVelocity_in_cm_per_s=1e5;  /* 1 km/sec */

  UnitTime_in_s= UnitLength_in_cm / UnitVelocity_in_cm_per_s;

  G=GRAVITY/pow(UnitLength_in_cm,3)*UnitMass_in_g*pow(UnitTime_in_s,2);
  Hubble = HUBBLE * UnitTime_in_s;

  //  cout<<"UnitLength_in_cm        ="<<UnitLength_in_cm<<endl;
  //  cout<<"UnitMass_in_g           ="<<UnitMass_in_g<<endl;
  //  cout<<"UnitVelocity_in_cm_per_s="<<UnitVelocity_in_cm_per_s<<endl;
  //  cout<<"UnitTime_in_s           ="<<UnitTime_in_s<<endl;
  //  cout<<" "<<endl;

  return 1;
}



 //==================================================
double PowerSpec(double kmag)
{

  switch(WhichSpectrum)
    {
    case 1:
      return PowerSpec_EH(kmag);
      break;
    case 2:
      return PowerSpec_Efstathiou(kmag);
      break;
    case 3:
      return PowerSpec_CMBFAST(kmag);
      break;
    case 4:
      return PowerSpec_BBKS(kmag);
      break;
    default:
      fprintf(stdout,"Not supported\n");  
    }

}

//==================================================
inline double PowerSpec_Efstathiou(double k)
{
  return Norm*pow(k,ns) / pow(1+pow(AA*k+pow(BB*k,1.5)+CC*CC*k*k,nu),2/nu);
}




//==================================================
inline double PowerSpec_BBKS(double k)
{
  return Norm*pow(k,ns) * pow(log(1.0+B1*k)/(B1*k),2)/ 
    pow(1+ B2*k + B3*B3*k*k + pow(B4*k,3) + pow(B5*k,4),0.5);
}

//==================================================
double   tk_eh(double k)  
{
  double   q,theta,ommh2,a,s,gamma,L0,C0;
  double   tmp;
  double   omegam, ombh2, hubble;

  /* other input parameters */
  hubble= HubbleParam;

  omegam= Omega;
  ombh2=  OmegaBaryon*HubbleParam*HubbleParam;

  //k*= 1000.0;    /* convert to h/Mpc */
  /*k*= HubbleParam;*/  /* convert to 1/Mpc */

  theta = 2.728/2.7;
  ommh2 = omegam*hubble*hubble;
  s     = 44.5*log(9.83/ommh2)/sqrt( 1.+10.*exp(0.75*log(ombh2)) )*hubble;      
  a     = 1.-0.328*log(431.*ommh2)*ombh2/ommh2
    +0.380*log(22.3*ommh2)*(ombh2/ommh2)*(ombh2/ommh2);
  gamma = a+(1.-a)/(1.+exp(4*log(0.43*k*s)) );
  gamma*= omegam*hubble;
  q     = k*theta*theta/gamma;
  L0    = log( 2.*exp(1.)+1.8*q );
  C0    = 14.2 + 731./(1.+62.5*q);
  tmp   = L0/(L0+C0*q*q);
  return(tmp);
}

 //==================================================
 inline double PowerSpec_EH(double k)
 {
	 return Norm*pow(k,ns)*pow( tk_eh(k), 2);
 }



 //==================================================
 inline double PowerSpec_CMBFAST(double k)
 {
	 //return Norm*pow(k,ns)*pow(transfunc_cmbfast(k), 2);
	 return Norm *transfunc_cmbfast(k);
 }


 //==================================================
 //==================================================
 double transfunc_cmbfast(double k)
 {
	 int i;
	 double lk;
	 double pow_index;

	 //k *= 1000.0; /* convert to h/Mpc */

	 lk=log10(k);

	 if(lk < xp[0]){
		 double dummy = (yp[2]-yp[1])/(xp[2]-xp[1])*(lk-xp[1])+yp[1];
		 return pow(10.,dummy);
	 } /* usually should not be executed */
	 if(lk > xp[np-1]){
		 double dummy = (yp[np-2]-yp[np-1])/(xp[np-2]-xp[np-1])*(lk-xp[np-1])+yp[np-1];
		 return pow(10.,dummy);
		 //return pow(10.,yp[np-1])*pow(k/pow(10.,xp[np-1]),-2.);
		 //return pow(10.,yp[np-1]);
	 }

	 splint(xp-1, yp-1, yp2-1, np, lk, &pow_index);

	 return pow(10.0,pow_index);
 }

 //==================================================
void readCMB_and_do_spline()
{
  int i,iline;
  double yp1,ypn;
  char tfunc_table[100];
  int  errorFlag=0;
  FILE *fd;

  sprintf(tfunc_table,"%s",TransTable);

  fprintf(stdout,"Reading %s .\n",tfunc_table);
  iline=0;
  double dummy1,dummy2,dummy3,dummy4,dummy5;
  fd=fopen(tfunc_table,"r");
  if(fd != NULL){
    while(!feof(fd)){
      fscanf(fd,"%lf %lf\n",&xp[iline],&yp[iline]);
      //fprintf(stdout,"%g %g \n",xp[iline],yp[iline]);
      iline++; 
    }
    fclose(fd);
  }
  else{
    fprintf(stdout,"transfer function file %s not found.\n",tfunc_table);
    exit(1);
  }

  fprintf(stdout,"read in %d data points \n",iline);

  np=iline;

  for(i=0;i<iline;i++){
    xp[i]=log10(xp[i]);
    yp[i]=log10(yp[i]);
  }

  yp1 = 1.e31;
  ypn = 1.e31;

  spline(xp-1, yp-1, iline, yp1, ypn, yp2-1);

  //for(i=0;i<iline;i++)
  //fprintf(stdout,"%g %g \n",xp[i],yp[i]);

}

//==================================================
inline double TopHatSigma2(double R)
{
  r_tophat= R;

  return qromb(sigma2_int, 0, 1000/R);
}

//==================================================
double sigma2_int(double k)
{
  double kr,kr3,kr2,wf,x;

  kr=r_tophat*k;
  kr2=kr*kr;
  kr3=kr2*kr;

  if(kr<1e-8) return 0;

  wf=3*( sin(kr)/kr3-cos(kr)/kr2 ); 
  x=4*PI*k*k*wf*wf*PowerSpec(k)/pow(2*PI,3.);

  return x;
}


//==================================================
double GrowthFactor(double astart, double aend)
{
  return growth(aend)/growth(astart);
}


inline double Hubble_a(double a)
{
  double res;
	 
  res= sqrt(Omega/(a*a*a) + (1-Omega-OmegaLambda)/(a*a) + OmegaLambda/pow(a,3.0*(1.0+w)));

  return res;
}

double Omega_de(double a){
  double res;
  res = OmegaLambda/(Omega*pow(a,3*w)+OmegaLambda);
  return res;
}	

double coeff1(double a){
  double res;
  res = 0.5*(5-3*w*Omega_de(a));
  return res;
}

double coeff2(double a){
  double res;
  res = 1.5*(1-w)*Omega_de(a);
  return res;
}

double growth(double a)
{
  double hubble_a;

  if(w != -1){
    hubble_a= sqrt(Omega/(a*a*a) + (1-Omega-OmegaLambda)/(a*a) + OmegaLambda/pow(a,3.0*(1.0+w)));
  }else{
    hubble_a= sqrt(Omega/(a*a*a) + (1-Omega-OmegaLambda)/(a*a) + OmegaLambda);
  }

  switch(WhichGF){
  case 1:
    return hubble_a*qromb(growth_int, 0, a);
    break;
  case 2:
    return growth_for_any_w(a);
    break;
  default:
    fprintf(stdout,"Not supported\n"); 
  }
}

inline double growth_int(double a)
{
  return pow(a / (Omega + (1-Omega-OmegaLambda)*a + OmegaLambda*a*a*a), 1.5);
}

inline double F_Omega(double a)
{
  double omega_a;

  omega_a= Omega/(Omega + a*(1-Omega-OmegaLambda) + a*a*a*OmegaLambda/pow(a,3.0*(1.0+w)));

  return pow(omega_a, 0.55);
}

inline double delta_c_func()
{ 
  double res;
  res = 3.0*pow(12*PI,0.666666666)*(1.0-0.00123*log10(Omega_z))/20.0;
  return res;
}

inline double efnn(double x, double rsphere)
{
  double rk,res;

  rk=exp(x);
  res=rk*var4(rk, rsphere);

  return res;
}


inline double var4(double x, double rsphere)
{
  double rk,res,pspec,xk;

  res=ddweight(x, rsphere)*PowerSpec(x)*bling*bling;

  return res;
}

inline double efn(double x, double rsphere)
{
  double rk,res;

  rk=exp(x);
  res=rk*var3(rk, rsphere);

  return res;
}

inline double var3(double x, double rsphere)
{
  double rk,res,pspec,xk;
  res=dweight(x, rsphere)*PowerSpec(x)*bling*bling;

  return res;
}

inline double var2(double x, double rsphere)
{
  double res,pspec,xk;
  res=weight(x, rsphere)*PowerSpec(x)*bling*bling;

  return res;
}

inline double evar2(double x, double rsphere)
{
  double rk,res;

  rk = exp(x);
  res = var2(rk, rsphere)*rk;

  return res;
}

inline double weight(double x, double rsphere) 
{
  // Tophat filter * k^2
  // there appears to be appreciable differences between C-version and fortran
  double y,res,yinv,yinv2;

  y=rsphere*x;
  yinv=1.0/y;
  yinv2=yinv*yinv;

  res=36.0*M_PI*x*x* pow((sin(y)/y-cos(y))*yinv2, 2)/pow(2*PI,3);

  if(sin(y)/y - cos(y) < TINY)
    res=0;

  return res;
}


inline double dweight(double x, double rsphere)  // derivative of weight
{
  double y,yinv,yinv2,yinv3,yinv4,res;

  y=rsphere*x;
  yinv=1.0/y;
  yinv2=yinv*yinv;
  yinv3=yinv*yinv2;
  yinv4=yinv*yinv3;


  res=M_PI*x*x*x*72.0*(3.0*cos(y)*yinv3 - 3.0*sin(y)*yinv4+sin(y)*yinv2)
    *(sin(y)*yinv3-cos(y)*yinv2)/pow(2*PI,3);

  return res;
}

inline double ddweight(double x, double rsphere)  // 2nd derivative of weight
{
  double y,yinv,yinv2,yinv3,yinv4,yinv5,res;

  y=rsphere*x;
  yinv=1.0/y;
  yinv2=yinv*yinv;
  yinv3=yinv*yinv2;
  yinv4=yinv*yinv3;
  yinv5=yinv*yinv4;


  res = M_PI*x*x*x*x*72.0*(cos(y)*yinv2-5.0*sin(y)*yinv3-12.0*cos(y)*yinv4+12.0*sin(y)*yinv5)
    *(sin(y)*yinv3-cos(y)*yinv2)/pow(2*PI,3);
  res+= M_PI*x*x*x*x*72.0*(3.0*cos(y)*yinv3 - 3.0*sin(y)*yinv4+sin(y)*yinv2)
    *(3.0*cos(y)*yinv3 - 3.0*sin(y)*yinv4+sin(y)*yinv2)/pow(2*PI,3);

  return res;
}

double dndlogm(double logm)
{
  // Number density per interval in log_{10} in mass

  double rho_crit=2.7754e11;

  double result,rm,rmass,sig,ans,fract;
  double r, res;


  rm = pow(10,logm);
  rmass = rm/rho_crit/Omega;

  sig = sigma_m(rmass, &r);

  ans = dlogdsigma(rmass, r, sig);

  fract = fract_mass(sig);
  if(OPT==4){
    double nu=delta_c_func()/sig;
    double n_eff=-6*ans-3;

    fract=fract*exp(-0.03/pow(n_eff+3,2)/pow(nu,0.6));
  }
  res = Rlog_e10*(rho_crit*Omega/rm*fabs(ans)*fract);

  //cout<<"rm="<<rm<<" ans="<<ans<<" sig="<<sig<<" fract="<<fract<<endl;
  return res;
}

double sigma_m(double m, double *rsphere_return)
{
  //   Use unit of mass where 1h^{-1}Mpc^3 has mass 1

  double res,rsphere,ling,aend;

  rsphere = pow((3.*m/4./M_PI), 0.33333333333);
  *rsphere_return = rsphere;

  res = sqrt(unnsigma(rsphere));

  return res;
}


double dlogdsigma(double mass, double rsphere, double sigma)
{
  double rho_crit=2.7754e11;
  double usig,ling,res,aend;

  res=sigdsigdr(rsphere);

  res= res * (mass/sigma/sigma)/4.0/M_PI/rsphere/rsphere;

  return res;
}

double sigdsigdr(double rsphere)
{
  int i;
  double dxk=0.005,xk;
  double sum=0,res;

  for(i=0;i<=8000;i++){
    xk=-20.0 + dxk*i;
    sum += efn(xk, rsphere)*dxk;
  }

  res = 0.5*sum;
  return res;
}

double dsig2dr2(double rsphere)
{
  int i;
  double dxk=0.0025,xk;
  double sum=0,res;

  for(i=0;i<=16000;i++){
    xk=-20.0 + dxk*i;
    sum += efn(xk, rsphere)*dxk;
  }

  res = sum;
  return res;
}

double fract_mass(double sig)  //  f(\ln\sigma^{-1})
{
  double delta=1.68647;
  double sqrt_two_over_pi=0.79788456;
  double sqrt_two=1.414213562;
  double nu_prime,s,fract,fdeltac;
  int iselect, ivalid;

  fdeltac=delta_c_func();

  double Nu = fdeltac/sig;

  switch(OPT)
    {
    case 0:{ // Jenkins et al 2000 fitting formula
      //fract=Jenkins(sig);

      //s = log(1./sig);
      //fract = 0.315*exp(-pow(fabs(s+0.61),3.8));      

      //if (s > -1.2 &&  s < 1.05)
      //	ivalid = 1;
      //else
      //	ivalid=0;
      //break;

      //Jenkins et al for FOF0.164 LCDM
      s = log(1./sig);
      fract = 0.301*exp(-pow(fabs(s+0.64),3.88));      

      if (s > -1.2 &&  s < 1.05)
	ivalid = 1;
      else
	ivalid=0;
      break;
    }
    case 1:{ // Sheth-Tormen
      //fract=Sheth_Tormen(nu);

      nu_prime = sqrt(0.707)*Nu;
      fract = 0.3222*sqrt_two_over_pi*nu_prime*exp(-nu_prime*nu_prime/2.)
	*(1.+ 1./pow(nu_prime,0.6));

      ivalid = 1;
      break;
    }
    case 2:{ // Press-Schechter
      //fract=Press-Schechter(nu);

      fract = sqrt_two_over_pi * Nu * exp(-Nu*Nu/2.);
      ivalid = 1; 
      break;
    }
    case 3:{ // Yahagi et al;

      fract = 0.298*(1.0+pow(0.893*Nu/sqrt_two, 1.39))*pow(Nu,0.408)
	*exp(-0.893*0.893*Nu*Nu/2.);

      ivalid = 1; 
      break;
    }
    case 4:{ // Reed et al 2007 fitting formula

      s=log(1./sig);
      nu_prime=sqrt(0.707)*Nu;
      double G1=exp(-pow((s-0.40)/0.6,2.0)/2);
      double G2=exp(-pow((s-0.75)/0.2,2.0)/2);
      fract= 0.3222*sqrt(2/PI)*nu_prime*exp(-Nu*Nu*0.746/2.)
	*(1.+1./pow(nu_prime,0.6)+0.6*G1+0.4*G2);

      if(s>-1.2 && s<1.05)
	ivalid=1;
      else
	ivalid=0;
      break;
    }
    case 5:{ // Bhattacharya et al 2011 fitting formula
      double A,B,p,q;
      A = 0.333/pow(1+calc_z,0.11);
      B = 0.788/pow(1+calc_z,0.01);
      p = 0.807;
      q = 1.795;

      fract = A*sqrt(2/PI)*exp(-0.5*B*Nu*Nu)
	*(1+pow(1/B/Nu/Nu,p))*pow(sqrt(B)*Nu,q);
      ivalid = 1;
      break;
    }
    case 6:{ 
      // Tinker et al 2008

      double A,a,b,c;
      double Delta=200.0;
      double alpha = pow(10.,-pow(0.75/(log10(Delta/75.)), 1.2));
      A = 0.186*pow(1.+calc_z, -0.14);
      a = 1.47*pow(1.+calc_z,-0.06);
      b = 2.57*pow(1.+calc_z, -alpha);
      c = 1.19;
      fract = A*(1+pow(sig/b,-a))*exp(-c/sig/sig);

      // Tinker et al 2010 (http://arxiv.org/pdf/1001.3162v2.pdf)
      /*
	Nu = 1.69/sig;
	double alp, bet, gam, phi, eta;
	alp = 0.368;
	bet = 0.589;
	gam = 0.864;
	phi = -0.729;
	eta = -0.243;

	bet = bet * pow(1.+calc_z, 0.20);
	gam = gam * pow(1.+calc_z, -0.01);
	phi = phi * pow(1.+calc_z, -0.08);
	eta = eta * pow(1.+calc_z, 0.27);

	fract = alp * Nu * (1+pow(bet*Nu, -2*phi))*pow(Nu, 2*eta)*exp(-0.5*gam*Nu*Nu);
      */
      ivalid = 1;
      break;
    }
  }

  return fract;

}

double unnsigma(double rsphere)
{
  int i,j,k;
  double dxk=0.01,xk;
  double sum=0;

  for(i=0;i<=4000;i++){
    xk=-20.0 + dxk*i;
    sum += evar2(xk, rsphere)*dxk;
  }

  return sum;
}

double H_z(double z){
  double res;
  res=100*HubbleParam*sqrt(Omega*pow((1+z),3)+OmegaLambda*pow(1+z,3*(1+w)));
  return res;
}

double H_z_1(double z,void *params){
  double res;
  res =1/H_z(z);
  return res;
}

double chi(double z){
  double result,abserr;
  double param=0;
  gsl_integration_workspace *w_gsl = gsl_integration_workspace_alloc(1000);
  gsl_function F;

  F.function=&H_z_1;
  F.params=&param;
  if(z>1e-10){
    gsl_integration_qag(&F,0,z,0,1e-7,1000,1,w_gsl,&result, &abserr);
  }else{
    result = 0;
  }

  gsl_integration_workspace_free(w_gsl);
  return C*result*HubbleParam; //Mpc/h
}	


double delta_v(double z){
  //double x=pow((1.0/Omega-1.0),0.3333333)*pow(1+z,-1.0);
  //double Omega_m=1.0/(1.0+pow(x,3.0));
  //return 18.0*pow(PI,2.0)*(1.0+0.4093*pow(x,2.7152));
  double Omz = (Omega*pow(1.0+z,3.0))/(Omega*pow(1.0+z,3.0)+OmegaLambda);
  return 18.0*M_PI*M_PI+82.*(Omz-1)-39.0*(Omz-1)*(Omz-1);
}

double r_vir(double z,double M){
  double rho_crit,logr;
  rho_crit=2.7754e11*(Omega*pow(1.0+z,3.0) + OmegaLambda);
  logr=(log10(3*M)-log10(4.0*PI*rho_crit*delta_v(z)))/3.0;
  return pow(10.0,logr);
}

double r_delta(double z, double Mass, double delta){
  double rho_crit,logr;

  //This difinition is the same as that in Tinker et al
  double rhom=2.7754e11*Omega*pow(1+z,3);
  logr=(log10(3*Mass)-log10(4.0*PI*rhom*delta))/3.0;
  return pow(10.0,logr);

  //This difinition is the same as that in Prada et al
  //rho_crit = 2.7754e11*(Omega*pow(1+z,3)+OmegaLambda);
  //logr=(log10(3*Mass)-log10(4.0*PI*rho_crit*delta))/3.0;
  //return pow(10.0,logr);

}

double halo_bias(double logm){ //based on peak-background split
  double delta=1.68647;
  double sqrt_two_over_pi=0.79788456;
  double sqrt_two=1.414213562;
  double fdeltac;
  double rho_crit=2.7754e11;
  double rm,rmass,sig,ans,fract;
  double r, res;

  fdeltac=delta_c_func();

  rm = pow(10,logm);
  rmass = rm/rho_crit/Omega;

  sig = sigma_m(rmass, &r);

  double Nu = fdeltac/sig;

  //For press-schechter
  //res = 1.+ (nu*nu-1.)/fdeltac;

  if(OPT==6){
    //Tinker et al 2010 arXiv:1001.3162
    Nu = 1.686/sig;
    double a1,a2,b1,b2,c1,c2;
    double Delta=200;
    double y = log10(Delta);
    a1 = 1.0+0.24*y*exp(-pow(4./y,4.));
    a2 = 0.44*y-0.88;
    b1 = 0.183;
    b2 = 1.5;
    c1 = 0.019+0.107*y+0.19*exp(-pow(4./y,4.));
    c2 = 2.4;
    res = 1.0-a1*(pow(Nu,a2))/(pow(Nu,a2)+pow(1.686,a2))+b1*pow(Nu,b2)+c1*pow(Nu,c2);
     
    //Tinker et al 2005 https://arxiv.org/abs/astro-ph/0411777
    /*
      double a = 0.707;
      double b = 0.35;
      double c = 0.80;
      Nu = delta/sig;
      res = 1.0 + 1./sqrt(a)/delta * (sqrt(a)*(a*Nu*Nu)+sqrt(a)*b*pow(a*Nu*Nu, 1-c) - pow(a*Nu*Nu,c)/(pow(a*Nu*Nu,c)+b*(1-c)*(1-0.5*c)));
    */

  }
  if(OPT==5){

    Nu = fdeltac/sig;
    double A,B,p,q;
    A = 0.333/pow(1+calc_z,0.11);
    B = 0.788/pow(1+calc_z,0.01);
    p = 0.807;
    q = 1.795;

    double X=B*Nu*Nu;

    res = 1.0 + 2*X/fdeltac*(0.5 - 1/X*(q/2+(-p+q/2)*pow(X,-p))/(1+pow(X,-p)));

  }
  if(OPT==1){
    //For sheth-tormen
    Nu = fdeltac/sig;
    double a = 0.75,p=0.3;
    res = 1.0 + (a*Nu*Nu-1.)/fdeltac + (2.*p)/fdeltac/(1.+pow(a*Nu*Nu,p));
  }

  return res;
}

double halo_bias_fast(double logm, double z){
  double lz = log10(z);
  double pow_index;

  if( z   < pow(10., tab_z_mf[0]) || lz   > pow(10., tab_z_mf[NPOINTS-1]) ){
    calc_z = z;
    bling = growth(1./(1.+calc_z))/growth(1.);
    return halo_bias(logm);
  }
  if(pow(10., logm) < pow(10.,  tab_m[0]) || pow(10., logm) > pow(10.,  tab_m[NPOINTS-1])){
    calc_z = z;
    bling = growth(1./(1.+calc_z))/growth(1.);
    return halo_bias(logm);
  }
  splin2(tab_z_mf-1, tab_m-1, tab_bh, err_bh, NPOINTS, NPOINTS, lz, logm, &pow_index);
  if(isnan(pow_index)){
    printf("fail to spline dndlogm at m=%e z=%e\n",pow(10.,logm),z);
    exit(1);
  }
  return pow(10.0,pow_index);
}

double RungeKutta(double a_in,double a){
  // u=D/a ,initial condition---> du/dlna=0,u=1
  int i,j;
  double h=(log(a)-log(a_in))/10000;
  double x=log(a_in);
  double u=1;
  double dudlna=0;
  double k0[2],k1[2],k2[2],k3[2];

  if(a_in==0){
    printf("you cannot solve calculate linear density growth eq.");
  }if(a == a_in){
    u=1;
  }else{
    for(i=0;i<10000;i++){

      k0[0]=h*dudlna;
      k0[1]=h*(-coeff1(exp(x))*dudlna-coeff2(exp(x))*(u));

      k1[0]=h*(dudlna+k0[1]/2);
      k1[1]=h*(-coeff1(exp(x+h/2))*(dudlna+k0[1]/2)-coeff2(exp(x+h/2))*(u+k0[0]/2));

      k2[0]=h*(dudlna+k1[1]/2);
      k2[1]=h*(-coeff1(exp(x+h/2))*(dudlna+k1[1]/2)-coeff2(exp(x+h/2))*(u+k1[0]/2));

      k3[0]=h*(dudlna+k2[1]);
      k3[1]=h*(-coeff1(exp(x+h))*(dudlna+k2[1])-coeff2(exp(x+h))*(u+k2[0]));

      u = u + (k0[0]+2*k1[0]+2*k2[0]+k3[0])/6;
      dudlna = dudlna + (k0[1]+2*k1[1]+2*k2[1]+k3[1])/6;
      x = x+h;
    }
  }

  return a*u;
}

void growth_spline(){
  int i;
  double yp1,ypn;
  double da = (1.0-(1.0/(1.0+1088.2)))/NPOINTS;

  for(i=0;i<NPOINTS;i++){
    scale_f[i] = (double)(i+1)*da + 1.0/(1.0+1088.2);
    GF[i] = RungeKutta(1.0/(1.0+1088.2),scale_f[i]);
  }

  for(i=0;i<NPOINTS;i++){
    scale_f[i] = log10(scale_f[i]);
    GF[i] = log10(GF[i]);
  }

  yp1 = 1.e31;
  ypn = 1.e31;

  spline(scale_f-1, GF-1, NPOINTS, yp1, ypn, GF2-1);

}

double growth_for_any_w(double a){
  double la;
  double pow_index;

  la=log10(a);

  if(la < scale_f[0])
    return a; /* usually should not be executed */
  if(la > scale_f[NPOINTS-1]){
    double dummy = (GF[NPOINTS-2]-GF[NPOINTS-1])/(GF[NPOINTS-2]-GF[NPOINTS-1])*(la-scale_f[NPOINTS-1])+GF[NPOINTS-1];
    return pow(10.,dummy);
    //return pow(10.,yp[np-1]);
  }

  splint(scale_f-1, GF-1, GF2-1, NPOINTS, la, &pow_index);
  return pow(10.0,pow_index);
}

void stack_table_and_spline(){
  int i,j,k;
  double yp1,ypn;
  double dlogz = (log10(5.0)-(-3.))/(NPOINTS-1);
  double dlogz2 = (log10(0.9)-log10(0.1))/(NPOINTS-1);
  double dlogm = (17.-10.)/(NPOINTS-1);

  yp1 = 1.e31;
  ypn = 1.e31;

  for(i=0;i<NPOINTS;i++){
    tab_z[i] = i*dlogz + (-3.);
    tab_z_mf[i] = i*dlogz2 + log10(0.1);
    tab_chi[i] = chi(pow(10,tab_z[i]));
  }	

  for(i=0;i<NPOINTS;i++){
    tab_chi[i] = log10(tab_chi[i]);
  }

  printf("stock chi(z) data...\n");

  tab_mf=dmatrix(1,NPOINTS,1,NPOINTS);
  err_mf=dmatrix(1,NPOINTS,1,NPOINTS);
  tab_bh=dmatrix(1,NPOINTS,1,NPOINTS);
  err_bh=dmatrix(1,NPOINTS,1,NPOINTS);
  tab_nuM=dmatrix(1,NPOINTS,1,NPOINTS);
  err_nuM=dmatrix(1,NPOINTS,1,NPOINTS);
  tab_cv=dmatrix(1,NPOINTS,1,NPOINTS);
  err_cv=dmatrix(1,NPOINTS,1,NPOINTS);

  for(i=0;i<NPOINTS;i++){
    a_exp=1.0/(1.0+pow(10,tab_z_mf[i])); // expansion parameter at which MF is computed
    calc_z=1.0/a_exp -1.0;

    Delta_c=1.68647;

    Omega_z = Omega*pow(1.0+calc_z,3.0)/(Omega*pow(1.0+calc_z,3.0) + OmegaLambda);
    fdelta_c = delta_c_func();

    bling = growth(a_exp)/growth(1.);

    for(j=0;j<NPOINTS;j++){
      tab_m[j] = j*dlogm + 10.;
      tab_mf[i+1][j+1] = dndlogm(tab_m[j]);
      tab_bh[i+1][j+1] = halo_bias(tab_m[j]);
      //if(tab_mf[i+1][j+1] <1e-30)cout << calc_z << " " << pow(10.,tab_m[j]) << " " << tab_mf[i+1][j+1] << " " << tab_bh[i+1][j+1] << endl;
      tab_nuM[i+1][j+1] = nu_M(calc_z, pow(10., tab_m[j]));
    }
  }
  printf("stock dn/dlogm data...\n");
  printf("stock halo bias data...\n");

  for(i=1;i<=NPOINTS;i++){
    for(j=1;j<=NPOINTS;j++){
      tab_mf[i][j] = log10(tab_mf[i][j]);
      tab_bh[i][j] = log10(tab_bh[i][j]);
      tab_nuM[i][j] = log10(tab_nuM[i][j]);
    }
  }
  spline(tab_z-1, tab_chi-1, NPOINTS, yp1, ypn, err_chi-1);
  printf("set spline chi(z)...\n");
  splie2(tab_z_mf-1, tab_m-1,tab_mf,NPOINTS,NPOINTS,err_mf);
  printf("set spline dn/dlogm...\n");
  splie2(tab_z_mf-1, tab_m-1,tab_bh,NPOINTS,NPOINTS,err_bh);
  printf("set spline halo bias...\n");
  splie2(tab_z_mf-1, tab_m-1,tab_nuM,NPOINTS,NPOINTS,err_nuM);
  printf("set spline nu(z, M)...\n");

  for(i=0;i<NPOINTS;i++){
    for(j=0;j<NPOINTS;j++){
      tab_cv[i+1][j+1] = c_200c_DK15(pow(10,tab_z_mf[i]), pow(10., tab_m[j]));
      //cout << pow(10,tab_z_mf[i]) << " " << pow(10., tab_m[j]) << " " << tab_cv[i+1][j+1] << endl;
      tab_cv[i+1][j+1] = log10(tab_cv[i+1][j+1]);
    }
  }

  splie2(tab_z_mf-1, tab_m-1,tab_cv,NPOINTS,NPOINTS,err_cv);
  printf("set spline c200c(z, M)...\n");

  for(i=0;i<NPOINTS;i++){
    tab_mf_norm[i] = mean_ngal(pow(10,tab_z_mf[i]));
    tab_mf_norm[i] = log10(tab_mf_norm[i]);
  }
  spline(tab_z_mf-1, tab_mf_norm-1, NPOINTS, yp1, ypn, err_mf_norm-1);

  //exit(1);

}

void splie2(double x1a[], double x2a[], double **ya, int m, int n, double **y2a){
  int j;
  for(j=1;j<=m;j++){
    spline(x2a,ya[j],n,1.e31,1.e31,y2a[j]);
  }
}
void splin2(double x1a[], double x2a[], double **ya, double **y2a, int m, int n, double x1, double x2, double *y){
  int j;
  double *ytmp,*yytmp;

  ytmp=dvector(1,n);
  yytmp=dvector(1,n);

  for(j=1;j<=m;j++){
    splint(x2a,ya[j],y2a[j],n,x2,&yytmp[j]);
  }
  spline(x1a,yytmp,m,1.e31,1.e31,ytmp);
  splint(x1a,yytmp,ytmp,m,x1,y);
  free_dvector(yytmp,1,n);
  free_dvector(ytmp,1,n);
}

void splie3(double x1a[], double x2a[], double x3a[], double ***ya, int l, int m, int n, double ***y2a){
  int j;
  for(j=1;j<=l;j++){
    splie2(x2a,x3a,ya[j],m,n,y2a[j]);
  }
}
void splin3(double x1a[], double x2a[], double x3a[], double ***ya, double ***y2a, int l, int m, int n, double x1, double x2, double x3, double *y){
  int j;
  double *ytmp,*yytmp;

  ytmp=dvector(1,l);
  yytmp=dvector(1,l);

  for(j=1;j<=l;j++){
    splin2(x2a,x3a,ya[j],y2a[j],m,n,x2,x3,&yytmp[j]);
  }
  spline(x1a,yytmp,l,1.e31,1.e31,ytmp);
  splint(x1a,yytmp,ytmp,l,x1,y);
  free_dvector(yytmp,1,l);
  free_dvector(ytmp,1,l);
}

double chi_fast(double z){
  double lz;
  double pow_index;

  lz=log10(z);

  if(lz < tab_z[0])
    return 1e-30; /* usually should not be executed */
  if(lz > tab_z[NPOINTS-1]){
    return chi(z);
  }

  splint(tab_z-1, tab_chi-1, err_chi-1, NPOINTS, lz, &pow_index);
  //printf("OK.spline chi(z)\n");
  return pow(10.0,pow_index);
}

double dndlogm_fast(double logm, double z){
  double lz = log10(z);
  double pow_index;

  if( z   < pow(10., tab_z_mf[0]) || lz   > pow(10., tab_z_mf[NPOINTS-1]) ){
    calc_z = z;
    bling = growth(1./(1.+calc_z))/growth(1.);
    return dndlogm(logm);
  }
  if(pow(10., logm) < pow(10.,  tab_m[0]) || pow(10., logm) > pow(10.,  tab_m[NPOINTS-1])){
    calc_z = z;
    bling = growth(1./(1.+calc_z))/growth(1.);
    return dndlogm(logm);
  }

  splin2(tab_z_mf-1, tab_m-1, tab_mf, err_mf, NPOINTS, NPOINTS, lz, logm, &pow_index);
  //printf("OK.spline dn/dlogm\n");
  if(isnan(pow_index)){
    printf("fail to spline dndlogm at m=%e z=%e\n",pow(10.,logm),z);
    return 1;
  }
  return pow(10.0,pow_index);
}

double P_nonlinear(double z, double k){

  double ksig ,n_eff, C_halo;
  double inv_D = (growth(1.0)/growth(1./(1+z)))*(growth(1.0)/growth(1./(1+z)));
  double param[3];
  set_halofit_param(z, param);
  ksig = param[0];
  n_eff = param[1];
  C_halo = param[2];
  double a_n,b_n,c_n,alpha_n,gamma_n,beta_n,nu_n,mu_n;
  if(OPT_fit==1){
    a_n = pow(10., 1.4861 + 1.8369*n_eff + 1.6762*n_eff*n_eff + 0.7940*n_eff*n_eff*n_eff + 0.1670*n_eff*n_eff*n_eff*n_eff -0.6206*C_halo);
    b_n = pow(10.,0.9463 + 0.9466*n_eff + 0.3084*n_eff*n_eff - 0.9400*C_halo);
    c_n = pow(10.,-0.2807 + 0.6669*n_eff + 0.3214*n_eff*n_eff - 0.0793*C_halo);
    gamma_n = 0.8649 + 0.2989*n_eff + 0.1631*C_halo;
    alpha_n = 1.3884 + 0.3700*n_eff - 0.1452*n_eff*n_eff;
    beta_n = 0.8291 + 0.9854*n_eff + 0.3401*n_eff*n_eff;
    mu_n = pow(10.,-3.5442 + 0.1908*n_eff);
    nu_n = pow(10., 0.9589 + 1.2857*n_eff);
  }
  if(OPT_fit==2){
    double om_w = OmegaLambda*pow(1+z,3*w)/(Omega+OmegaLambda*pow(1+z,3*w)); 
    a_n = pow(10.,1.5222 + 2.8553*n_eff + 2.3706*n_eff*n_eff + 0.9903*n_eff*n_eff*n_eff + 0.2250*n_eff*n_eff*n_eff*n_eff -0.6038*C_halo +0.1749*om_w*(1.+w));
    b_n = pow(10.,-0.5642 + 0.5864*n_eff + 0.5716*n_eff*n_eff - 1.5474*C_halo +0.2279*om_w*(1.+w));
    c_n = pow(10.,0.3698 + 2.0404*n_eff + 0.8161*n_eff*n_eff + 0.5869*C_halo);
    gamma_n = 0.1971 - 0.0843*n_eff + 0.8460*C_halo;
    alpha_n = fabs(6.0835 + 1.3373*n_eff - 0.1959*n_eff*n_eff - 5.5274*C_halo);
    beta_n = 2.0379 - 0.7354*n_eff + 0.3157*n_eff*n_eff + 1.2490*n_eff*n_eff*n_eff + 0.3980*n_eff*n_eff*n_eff*n_eff - 0.1682*C_halo;
    mu_n = 0;
    nu_n = pow(10., 5.2105 + 3.6902*n_eff);
  }
  //printf("%e %e %e %e %e %e %e\n",a_n,b_n,c_n,alpha_n,gamma_n,beta_n,nu_n,mu_n);


  double delta_L,delta_H,delta_Q;
  //only support flat universe
  double omz;
  omz = Omega/(Omega+OmegaLambda*pow(1+z,3*w));

  double f1 = pow(omz,-0.0307);
  double f2 = pow(omz,-0.0585);
  double f3 = pow(omz,0.0743);

  double fac =4*PI*k*k*k/(2*PI)/(2*PI)/(2*PI);
  delta_L = fac*PowerSpec(k)/inv_D;

  double y = k/ksig;
  delta_Q = delta_L*(pow(1.+delta_L,beta_n)/(1.+alpha_n*delta_L))*exp(-y/4-y*y/8);

  delta_H = a_n*pow(y,3*f1)/(1+b_n*pow(y,f2)+pow(c_n*f3*y,3-gamma_n));
  delta_H = delta_H/(1+mu_n/y+nu_n/y/y);

  return (delta_Q+delta_H)/fac;
}

void set_halofit_param(double z, double *param){
  int i;
  double inv_D = (growth(1.0)/growth(1./(1+z)))*(growth(1.0)/growth(1./(1+z)));
  thre = inv_D;
  //double init = pow(10.,-1.065+4.332e-1*(1.+z)-2.516e-2*pow(1.+z,2)+9.069e-4*pow(1.+z,3));
  //double R0=solver(z);
  //printf("z=%e\n",z);
  double R0=solver(z);
  //check the solver
  //double dummy=0;
  //printf("check_solver %e\n",sigma2_gauss(log(R0),&dummy));

  param[0] = 1./R0; //k_sig
  param[1] = neff(R0);
  param[2] = C_halofit(R0);

}

double solver(double z){
  int status; 
  int iter = 0, max_iter = 100;
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;
  double r = 0.;
  double x_lo = -5.0, x_hi = 5.0;//for gravitino or WDM, we should set x_lo=-7.0, otherwise x_lo=-5. to calculate faster
  gsl_function F; 
  double params = z;
  F.function = &sigma2_gauss; 
  F.params = &params; 
  T = gsl_root_fsolver_brent; 
  s = gsl_root_fsolver_alloc(T);

  gsl_root_fsolver_set(s, &F, x_lo, x_hi); 

  do {
    iter++;
    status = gsl_root_fsolver_iterate(s); 
    r = gsl_root_fsolver_root(s); 
    x_lo = gsl_root_fsolver_x_lower(s); 
    x_hi = gsl_root_fsolver_x_upper(s); 
    status = gsl_root_test_interval (x_lo, x_hi, 0, 1e-7); 
    //if (status == GSL_SUCCESS) printf ("Converged:\n");
    //printf("%5d [%.7f, %.7f] %.7f %.7f\n",iter, x_lo, x_hi, r,  x_hi-x_lo);

  } while (status == GSL_CONTINUE && iter < max_iter);

  gsl_root_fsolver_free(s);

  double R0 = exp(r);

  //printf("Rsig = %e\n",R0);
  return R0;
}

double get_delta_k(double k){
  double lk;
  double pow_index;
  double fac = k*k*k*4*PI/(2*PI)/(2*PI)/(2*PI);

  return fac*PowerSpec(k);

}

double sigma2_gauss_int(double lnk, void *params){
  double res;
  double R = *(double *)params;

  double k = exp(lnk);

  res = get_delta_k(k)*exp(-k*k*R*R);
  return res;
}

double sigma2_gauss(double lnR, void *params){
  double result,abserr;
  double params_int=exp(lnR);
  double z = *(double *)params;
  double inv_D = (growth(1.0)/growth(1./(1+z)))*(growth(1.0)/growth(1./(1+z)));
  size_t neval;
  gsl_function F;

  //gsl_integration_workspace *wgsl = gsl_integration_workspace_alloc(1000);

  //F.function=&sigma2_gauss_int;
  //F.params=&params_int;
  //gsl_integration_qng(&F,-3,3,0,1e-6, &result, &abserr,&neval);
  //gsl_integration_qag(&F,log(calc_z),log(zLSS),0,1e-7,neval,6,wgsl,&result, &abserr);
  //gsl_integration_workspace_free(wgsl);

  /*int Nint = 5000;

    double dlnk = 10.0/Nint;
    result = 0;
    for(int i=0;i<Nint;i++){
    result += dlnk*sigma2_gauss_int(i*dlnk-5.0,&params_int);
    }*/
  //res = qromb(integral_Pkappa,log(calc_z),log(zLSS));

  return sigma2_gauss_fast(params_int)-inv_D;
}

double dsig_dR_int(double lnk, void *params){
  double res;
  double R = *(double *)params;

  double k = exp(lnk);

  res = get_delta_k(k)*exp(-k*k*R*R)*(-2*k*k*R);
  return res;
}
double dsig_dR(double R){
  double result,abserr;
  double params=R;
  size_t neval;
  gsl_function F;

  //gsl_integration_workspace *wgsl = gsl_integration_workspace_alloc(1000);

  //F.function=&dsig_dR_int;
  //F.params=&params;
  //gsl_integration_qng(&F,-3,3,0,1e-6, &result, &abserr,&neval);
  //gsl_integration_qag(&F,log(calc_z),log(zLSS),0,1e-7,neval,6,wgsl,&result, &abserr);
  //gsl_integration_workspace_free(wgsl);

  int Nint = 5000;

  double dlnk = 20.0/Nint;
  result = 0;
  for(int i=0;i<Nint;i++){
    result += dlnk*dsig_dR_int(i*dlnk-10.0,&params);
  }
  //res = qromb(integral_Pkappa,log(calc_z),log(zLSS));

  return result;
}
double d2sig_dR2_int(double lnk, void *params){
  double res;
  double R = *(double *)params;

  double k = exp(lnk);
  double y = k*R;

  res = get_delta_k(k)*exp(-y*y)*(y*y-y*y*y*y);

  return res;
}
double d2sig_dR2(double R){
  double result,abserr;
  double params=R;
  size_t neval;
  gsl_function F;

  //gsl_integration_workspace *wgsl = gsl_integration_workspace_alloc(1000);

  //F.function=&d2sig_dR2_int;
  //F.params=&params;
  //gsl_integration_qng(&F,-3,3,0,1e-6, &result, &abserr,&neval);
  //gsl_integration_qag(&F,log(calc_z),log(zLSS),0,1e-7,neval,6,wgsl,&result, &abserr);
  //gsl_integration_workspace_free(wgsl);

  int Nint = 5000;

  double dlnk = 20.0/Nint;
  result = 0;
  for(int i=0;i<Nint;i++){
    result += dlnk*d2sig_dR2_int(i*dlnk-10.0,&params);
  }
  //res = qromb(integral_Pkappa,log(calc_z),log(zLSS));
  //if(result < 0){printf("R=%e d^2sig/dR^2 = %e\n",R,result);}

  return result;
}
double neff(double R){
  double sig2 = thre;
  double res = -R*dsig_dR_fast(R)/sig2;

  return res -3.0;
}

double C_halofit(double R){
  double sig2 = thre;
  double n_eff = neff(R);

  double res = (3.+n_eff)*(3+n_eff) +4./thre*d2sig_dR2_fast(R); //cf .smith et al

  return res;

}

void stack_data_and_spline_Pk(){
  double logR = (3.-(-3.))/NPOINTS;
  int i;

  for(i=0;i<NPOINTS;i++){
    tab_R[i] = (double)(i)*logR -3.;
    tab_dsdr[i] = dsig_dR(pow(10.,tab_R[i]));
    tab_ds2dr2[i] = d2sig_dR2(pow(10.,tab_R[i]));
    //printf("%e %e %e\n",tab_R[i],tab_dsdr[i],tab_ds2dr2[i]);
    tab_dsdr[i] = log10(-tab_dsdr[i]);
    tab_ds2dr2[i] = log10(10+tab_ds2dr2[i]);
    //printf("%e %e %e\n",tab_R[i],tab_dsdr[i],tab_ds2dr2[i]);
  }

  int Nint = 5000;

  double dlnk = 20.0/Nint;
  double params_int;
  for(i=0;i<NPOINTS;i++){
    tab_sig2[i]=0.;
    params_int = pow(10.,tab_R[i]);
    for(int j=0;j<Nint;j++){
      tab_sig2[i] += dlnk*sigma2_gauss_int(j*dlnk-10.0,&params_int);
    }
    tab_sig2[i] = log10(tab_sig2[i]);
  }

  double yp1 = 1.e31;
  double ypn = 1.e31;

  spline(tab_R-1, tab_sig2-1, NPOINTS, yp1, ypn, err_sig2-1);
  printf("spline sigma2 ... \n");
  spline(tab_R-1, tab_dsdr-1, NPOINTS, yp1, ypn, err_dsdr-1);
  printf("spline dsig/dr ... \n");
  spline(tab_R-1, tab_ds2dr2-1, NPOINTS, yp1, ypn, err_ds2dr2-1);
  printf("spline d2sig/dr2 ... \n");
}

double dsig_dR_fast(double R){
  double lR;
  double pow_index;

  lR=log10(R);

  if(lR < tab_R[0] || lR > tab_R[NPOINTS-1]){
    return dsig_dR(R);
  }else{
    splint(tab_R-1, tab_dsdr-1, err_dsdr-1, NPOINTS, lR, &pow_index);
    return -pow(10.,pow_index);
  }
}
double d2sig_dR2_fast(double R){
  double lR;
  double pow_index;

  lR=log10(R);

  if(lR < tab_R[0] || lR > tab_R[NPOINTS-1]){
    return d2sig_dR2(R);
  }else{
    splint(tab_R-1, tab_ds2dr2-1, err_ds2dr2-1, NPOINTS, lR, &pow_index);
    return pow(10.,pow_index)-10.;
  }
}

double sigma2_gauss_fast(double R){
  double lR;
  double pow_index;

  lR=log10(R);

  if(lR < tab_R[0] || lR > tab_R[NPOINTS-1]){
    int Nint = 5000;
    double dlnk = 20.0/Nint;
    double result=0.;

    for(int i=0;i<Nint;i++){
      result += dlnk*sigma2_gauss_int(i*dlnk-10.0,&R);
    }
    return result;
  }else{
    splint(tab_R-1, tab_sig2-1, err_sig2-1, NPOINTS, lR, &pow_index);
    return pow(10.,pow_index);
  }
}

double Redshift(double x,double z_in){
  //x,x(z_in) [Mpc/h]  dx [Mpc]  
  int i;
  int nint = 5000;
  double dx=(x-chi(z_in))/HubbleParam/nint;
  double z=z_in,dz;

  for(i=0;i<nint;i++){
    dz = H_z(z)/C*dx;
    z += dz;
  }
  return z;
}

double HOD_gal(double M, int cen_or_sat){

  //https://arxiv.org/pdf/1407.1856.pdf
  double res;
  double mass=M;

  double logMmin = hod.logMmin;
  double sig = sqrt(hod.sig_sq);
  double fac = (log10(mass)-logMmin)/sig;

  res = 0.5*(1.+erf(fac));

  if(cen_or_sat == 0){//central
    return res;
  }else{
    double Mcut = hod.kappa * pow(10., logMmin);
    double M1 = pow(10., hod.logM1);
    double alpha = hod.alpha;
    if(mass < Mcut){return 0;}

    res = res*pow((mass-Mcut)/M1, alpha);
    return res;
  }
}

double mean_ngal_int(double logM, double z){
  double M= pow(10., logM);
  double mf = dndlogm_fast(logM, z);

  return (HOD_gal(M, 0)+HOD_gal(M, 1))*mf;	
}

double mean_ngal(double z){
  double res,abserr;
  double x_lo = log10(Mh_min);
  double x_hi = log10(Mh_max);

  int Nint = 10;
  //Romberg
  int i,j;
  double h[Nint];
  double s[Nint][Nint];
  for(i=1;i<=Nint;i++){h[i-1] = (x_hi-x_lo)/pow(2.,i-1);}

  s[0][0] = 0.5*h[0]*(mean_ngal_int(x_hi, z)+mean_ngal_int(x_lo, z));
  for(i=2;i<=Nint;i++){
    s[i-1][0] = s[i-2][0];
    for(j=1;j<=pow(2.,i-2);j++){
      s[i-1][0] += h[i-2]*mean_ngal_int(x_lo+(2*j-1)*h[i-1], z);
    }
    s[i-1][0] = 0.5*s[i-1][0];
  }

  for(i=2;i<=Nint;i++){
    for(j=2;j<=i;j++){
      s[i-1][j-1] = s[i-1][j-2]+(s[i-1][j-2]-s[i-2][j-2])/(pow(4.,j-1)-1);
    }
  }

  res = s[Nint-1][Nint-1];

  return res;
}

double mean_ngal_fast(double z){
  double lz;
  double pow_index;

  lz=log10(z);

  if(lz < tab_z_mf[0])
    return 0.0; /* usually should not be executed */
  if(lz > tab_z_mf[NPOINTS-1]){
    return 0.0;//window_kappa_dSigma(z);
  }

  splint(tab_z_mf-1, tab_mf_norm-1, err_mf_norm-1, NPOINTS, lz, &pow_index);
  //printf("OK.spline chi(z)\n");
  return pow(10.0,pow_index);
}

double dens_sat_int(double logM, double z){
  double M= pow(10., logM);
  double mf = dndlogm_fast(logM, z);

  return HOD_gal(M, 1)*mf;	
}

double dens_sat(double z){
  double res,abserr;
  double x_lo = log10(Mh_min);
  double x_hi = log10(Mh_max);

  int Nint = 12;
  //Romberg
  int i,j;
  double h[Nint];
  double s[Nint][Nint];
  for(i=1;i<=Nint;i++){h[i-1] = (x_hi-x_lo)/pow(2.,i-1);}

  s[0][0] = 0.5*h[0]*(dens_sat_int(x_hi, z)+dens_sat_int(x_lo, z));
  for(i=2;i<=Nint;i++){
    s[i-1][0] = s[i-2][0];
    for(j=1;j<=pow(2.,i-2);j++){
      s[i-1][0] += h[i-2]*dens_sat_int(x_lo+(2*j-1)*h[i-1], z);
    }
    s[i-1][0] = 0.5*s[i-1][0];
  }

  for(i=2;i<=Nint;i++){
    for(j=2;j<=i;j++){
      s[i-1][j-1] = s[i-1][j-2]+(s[i-1][j-2]-s[i-2][j-2])/(pow(4.,j-1)-1);
    }
  }

  res = s[Nint-1][Nint-1];

  return res;
}

double dlnP_dlnk(double lnk){
  double dlnk = 0.01;
  double k1 = exp(lnk-dlnk);
  double k2 = exp(lnk+dlnk);
  double P1 = PowerSpec_EH(k1);
  double P2 = PowerSpec_EH(k2);
  return log(P2/P1)/dlnk/2;
}

double nu_M(double z, double M){

  double rho_crit=2.7754e11;
  double r, rmass = M/rho_crit/Omega;
  bling=GrowthFactor(1.,1./(1.+z));
  double sig = sigma_m(rmass, &r);

  return 1.686/sig;
}

double nu_M_fast(double z, double M){
  double logm = log10(M);
  double lz = log10(z);
  double pow_index;

  splin2(tab_z_mf-1, tab_m-1, tab_nuM, err_nuM, NPOINTS, NPOINTS, lz, logm, &pow_index);
  //printf("OK.spline dn/dlogm\n");
  if(isnan(pow_index)){
    printf("fail to spline nu at m=%e z=%e\n",pow(10.,logm),z);
    exit(1);
  }
  return pow(10.0,pow_index);
}

double c_200c_DK15(double z, double M){ //https://arxiv.org/pdf/1407.4730.pdf
  double nu = nu_M_fast(z, M);
  double rho_crit=2.7754e11;
  double R = pow(3*M/rho_crit/Omega/(4*M_PI), 0.3333);

  //paramaters
  double kap = 0.69;
  double phi0 = 7.14;
  double phi1 = 1.60;
  double eta0 = 4.10;
  double eta1 = 0.75;
  double alpha = 1.40;
  double beta = 0.67;

  double n = dlnP_dlnk(log(2*M_PI/R*kap));

  double cmin = phi0 + phi1 *n;
  double numin = eta0 + eta1 *n;

  double res = cmin/2 * (pow(nu/numin, -alpha) + pow(nu/numin, beta));
  //cout << n << " " << res << endl;
  return res;
}

double c_200c_DK15_fast(double z, double M){
  double logm = log10(M);
  double lz = log10(z);
  double pow_index;

  splin2(tab_z_mf-1, tab_m-1, tab_cv, err_cv, NPOINTS, NPOINTS, lz, logm, &pow_index);
  //printf("OK.spline dn/dlogm\n");
  if(isnan(pow_index)){
    printf("fail to spline c200c at m=%e z=%e\n",pow(10.,logm),z);
    exit(1);
  }
  return pow(10.0,pow_index);
}

double get_eq_for_M_200c(double x, void *params){
  struct get_mvir *p = (struct get_mvir *)params;

  double delta = p->delta;
  double rd = p->rd; 
  double z = p->redshift;

  double rhom =2.7754e11*Omega*pow(1+z,3);
  double rhoc =2.7754e11*(Omega*pow(1+z,3)+OmegaLambda);

  //solve 3 rho_s (z, m200c) rs(z, m200c)^3 f(rd/rs(z, m200c)) = rho_m * delta * rd^3 where f(x) = ln(1+x)-x/(1+x)	
  double m200c = pow(10., x);
  double r200c = pow(10., (log10(3*m200c)-log10(4.0*PI*rhoc*200.))/3.0);
  double c200c = c_200c_DK15_fast(z, m200c);
  double rs = r200c/c200c;
  double rhos = m200c/(log(1.+c200c)-c200c/(1.+c200c))/(4*PI*rs*rs*rs);
  double cd  = rd/rs;

  double rhs = 3 * rhos * rs* rs* rs* (log(1+cd) -cd/(1+cd));
  double lfs = rhom * delta * rd * rd * rd;
  return rhs-lfs;
}

double M_delta_to_M_200c(double z, double mdelta, double delta){
  double rd = r_delta(z, mdelta, delta);

  //get mvir
  int status; 
  int iter = 0, max_iter = 100;
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;
  double r = 1.0;
  double x_lo = 5.0, x_hi = 20.0;
  gsl_function F;
  struct get_mvir params={delta, rd, z};

  //printf("%e %e\n",get_eq_for_M_vir(x_lo, &params), get_eq_for_M_vir(x_hi, &params));

  F.function = &get_eq_for_M_200c;
  F.params = &params; 
  T = gsl_root_fsolver_brent; 
  s = gsl_root_fsolver_alloc(T);

  gsl_root_fsolver_set(s, &F, x_lo, x_hi); 

  do {
    iter++;
    status = gsl_root_fsolver_iterate(s); 
    r = gsl_root_fsolver_root(s); 
    x_lo = gsl_root_fsolver_x_lower(s); 
    x_hi = gsl_root_fsolver_x_upper(s); 
    status = gsl_root_test_interval (x_lo, x_hi, 0, 1e-7); 
    //if (status == GSL_SUCCESS) printf ("Converged:\n");
    //printf("%5d [%.7f, %.7f] %.7f %.7f\n",iter, x_lo, x_hi, r,  x_hi-x_lo);

  } while (status == GSL_CONTINUE && iter < max_iter);

  gsl_root_fsolver_free(s);

  double m200c = pow(10., r);
  return m200c;
}

void save_SC_param(){

  /*FILE *fd;
    fd = fopen("SC_dens_velo.dat", "w");

    int i, nplot=1000;
    for(i=0;i<nplot;i++){
    double gamma = 0.0+(2*M_PI-0.0)/nplot*(i+0.5);
    double dens = SC_param_density(gamma);
    double velo = SC_param_velocity(gamma);
    fprintf(fd, "%e %e %e\n", gamma, dens, velo);
    }
    fclose(fd);
  */

  int i;
  for(i=0;i<1000;i++){
    double gamma = 0.0+(2*M_PI-0.0)/1000.0*(i+0.5);
    double dens = SC_param_density(gamma);
    double velo = SC_param_velocity(gamma);
    tab_SC_dens[i] = log10(dens);
    tab_SC_velo[i] = log10(1-velo);
  }
  
  double yp1 = 1.e31;
  double ypn = 1.e31;

  spline(tab_SC_dens-1, tab_SC_velo-1, 1000, yp1, ypn, err_SC_velo-1);
  
}

double SC_param_density(double gamma){
  return 9.0*pow(gamma-sin(gamma), 2)/2.0/pow(1-cos(gamma),3)-1;
}

double SC_param_velocity(double gamma){
  return 3.0*sin(gamma)*(gamma-sin(gamma))/pow(1-cos(gamma),2)-1;
}
double SC_param_velocity_fast(double density){
  int i;
  double lk;
  double pow_index;
		
  if(density < 0){
    return 1.0;
  }

  lk=log10(density);
	
  if(lk < tab_SC_dens[0]){
    pow_index = (tab_SC_velo[1]-tab_SC_velo[0])/(tab_SC_dens[1]-tab_SC_dens[0])*(lk-tab_SC_dens[0])+tab_SC_velo[0];
    return 1.0-pow(10.0, pow_index);
  } /* usually should not be executed */
  if(lk > tab_SC_dens[1000-1]){
    pow_index = (tab_SC_velo[1000-1]-tab_SC_velo[1000-2])/(tab_SC_dens[1000-1]-tab_SC_dens[1000-2])*(lk-tab_SC_dens[1000-2])+tab_SC_velo[1000-2];
    return 1.0-pow(10.0, pow_index);
  }
	
  splint(tab_SC_dens-1, tab_SC_velo-1, err_SC_velo-1, 1000, lk, &pow_index);
	
  return 1.0-pow(10.0,pow_index);
}

double mean_infall_velocity(double r, double density){
    
    
  /*
  double mu_lin = -H_z(zhalo)/HubbleParam/(1+zhalo) * r * F_Omega(1./(1.+zhalo)) * density / 3.0; // km/s  
  double mu_sc = H_z(zhalo)/HubbleParam/(1+zhalo) * r * F_Omega(1./(1.+zhalo)) * SC_param_velocity_fast(density) * exp(-pow(4.5/r/(1+density), 2)); // km/s

  double weight;
  if(r <= 4) weight = 1.0;
  if(4 <= r && r < 20.0) weight = (0.0-1.0)/(log(20.0)-log(4.0))*(log(r)-log(4.0)) + 1.0;
  if(20 <= r) weight = 0.0;

  return mu_lin * (1-weight) + mu_sc * weight;  
  */
  
  double deltac = 1.686;
  // Eq (B2) in https://arxiv.org/pdf/1305.5548.pdf
  double mu_sc = -H_z(zhalo)/HubbleParam/(1+zhalo)* r * F_Omega(1./(1.+zhalo)) / 3 * deltac * (pow(1+density, 1./deltac)-1.0);
  return mu_sc;  
  
}

double sigma_r_velocity(double r, double density, double mass1, double mass2){
  
  double rd1 = r_delta(zhalo, mass1, 200.0) * (1.0+zhalo);
  double rd2 = r_delta(zhalo, mass2, 200.0) * (1.0+zhalo);
  double rd0 = rd1 + rd2;

  double rho_scale = pow(r/(5.0)/sqrt(rd1), -4.0) + pow(r/(11.5)/sqrt(rd0), -1.3) + 0.50;
  double alpha = pow(r/35.0, 0.1);
  
  double f = F_Omega(1./(1.+zhalo)) / F_Omega(1.0);
  bling = growth(1./(1.+zhalo))/growth(1.);
  
  return 200.0 * pow(Omega/0.3, 0.6) * (Sigma8/0.8) * pow((1+density)/rho_scale, alpha) * f * bling /(1+zhalo) * H_z(zhalo) / H_z(0.0); 

}

//==================================================
double TopHatSigma2_NL(double R)
{
  r_tophat= R;

  double res,abserr;
  double x_lo = 0.0;
  double x_hi = 1000.0/R;

  int Nint = 10;
  //Romberg
  int i,j;
  double h[Nint];
  double s[Nint][Nint];
  for(i=1;i<=Nint;i++){h[i-1] = (x_hi-x_lo)/pow(2.,i-1);}

  s[0][0] = 0.5*h[0]*(sigma2_NL_int(x_hi, zhalo)+sigma2_NL_int(x_lo, zhalo));
  for(i=2;i<=Nint;i++){
    s[i-1][0] = s[i-2][0];
    for(j=1;j<=pow(2.,i-2);j++){
      s[i-1][0] += h[i-2]*sigma2_NL_int(x_lo+(2*j-1)*h[i-1], zhalo);
    }
    s[i-1][0] = 0.5*s[i-1][0];
  }

  for(i=2;i<=Nint;i++){
    for(j=2;j<=i;j++){
      s[i-1][j-1] = s[i-1][j-2]+(s[i-1][j-2]-s[i-2][j-2])/(pow(4.,j-1)-1);
    }
  }

  res = s[Nint-1][Nint-1];

  return res;
  
}

//==================================================
double sigma2_NL_int(double k, double z)
{
  double kr,kr3,kr2,wf,x;

  kr=r_tophat*k;
  kr2=kr*kr;
  kr3=kr2*kr;

  if(kr<1e-8) return 0;

  wf=3*( sin(kr)/kr3-cos(kr)/kr2 ); 
  x=4*PI*k*k*wf*wf*P_nonlinear(z, k)/pow(2*PI,3.);

  return x;
}

double PDF_mass_NL(double density, double r, double sigm, double mass1, double mass2){

  double bh1 = halo_bias_fast(log10(mass1), zhalo);
  double bh2 = halo_bias_fast(log10(mass2), zhalo);
  
  double sig1 = log(1+sigm);  
  double rd1 = r_delta(zhalo, mass1, 200.0) * (1.0+zhalo);  
  double pdf = 1.0/sqrt(2*M_PI*sig1)*exp(-pow(log(1+density)+sig1/2,2)/2/sig1);

  /*
  double out[3];
  get_modified_ln_params(sigm, out);
  double N = out[0];
  double A = out[1];
  double w2 = out[2];
  double pdf = N * exp(-pow(log(1+density)+w2*0.5,2)*(1.0+A/(1+density))/2/w2);
  */

  double rho_scale = 1.41*(bh1+bh2) + pow(r/(9.4*rd1), -2.2);
  bling = growth(1./(1.+zhalo))/growth(1.);

  return exp(-rho_scale/(1+density))*pdf;  

}

double PDF_halo_int_n(double density, double r, double sigm, double mass1, double mass2){  
  double delta = exp(density) -1; // density = log(1+delta) 
  double cond_pdf = PDF_mass_NL(delta, r, sigm, mass1, mass2);

  return cond_pdf;
}

double PDF_halo_int_d(double velocity, double density, double r, double sigm, double mass1, double mass2){

  double v = velocity;
  double delta = exp(density) -1; // density = log(1+delta)
  
  double mean = mean_infall_velocity(r, delta);
  double sig = sigma_r_velocity(r, delta, mass1, mass2);
  double gauss_pdf = exp(-pow(v-mean, 2)/2/(sig*sig))/sqrt(2*M_PI*sig*sig);
  double cond_pdf = PDF_mass_NL(delta, r, sigm, mass1, mass2);

  if(isnan(gauss_pdf)){
    fprintf(stdout, "found Nan in PDF_halo_int_d\n");
    fprintf(stdout, "v=%e, r=%e, delta=%e, M1=%e, M2=%e\n", v, r, delta, mass1, mass2);
    fprintf(stdout, "mean = %e, sig = %e\n", mean, sig);
    exit(1);
  }

  return gauss_pdf * cond_pdf;

}

double mean_overdesnity_ln(double r, double sigm, double mass1, double mass2){

  double sig1 = log(1+sigm);

  double res,abserr;
  double x_lo = -10.0*sqrt(sig1);
  double x_hi = +10.0*sqrt(sig1);

  int Nint = 10;
  //Romberg
  int i,j;
  double h[Nint];
  double s[Nint][Nint];
  for(i=1;i<=Nint;i++){h[i-1] = (x_hi-x_lo)/pow(2.,i-1);}

  s[0][0] = 0.5*h[0]*((exp(x_hi)-1)*PDF_halo_int_n(x_hi, r, sigm, mass1, mass2)+(exp(x_lo)-1)*PDF_halo_int_n(x_lo, r, sigm, mass1, mass2));
  for(i=2;i<=Nint;i++){
    s[i-1][0] = s[i-2][0];
    for(j=1;j<=pow(2.,i-2);j++){
      s[i-1][0] += h[i-2]*(exp(x_lo+(2*j-1)*h[i-1])-1)*PDF_halo_int_n(x_lo+(2*j-1)*h[i-1], r, sigm, mass1, mass2);
    }
    s[i-1][0] = 0.5*s[i-1][0];
  }

  for(i=2;i<=Nint;i++){
    for(j=2;j<=i;j++){
      s[i-1][j-1] = s[i-1][j-2]+(s[i-1][j-2]-s[i-2][j-2])/(pow(4.,j-1)-1);
    }
  }

  res = s[Nint-1][Nint-1];

  double norm;

  s[0][0] = 0.5*h[0]*(PDF_halo_int_n(x_hi, r, sigm, mass1, mass2)+PDF_halo_int_n(x_lo, r, sigm, mass1, mass2));
  for(i=2;i<=Nint;i++){
    s[i-1][0] = s[i-2][0];
    for(j=1;j<=pow(2.,i-2);j++){
      s[i-1][0] += h[i-2]*PDF_halo_int_n(x_lo+(2*j-1)*h[i-1], r, sigm, mass1, mass2);
    }
    s[i-1][0] = 0.5*s[i-1][0];
  }

  for(i=2;i<=Nint;i++){
    for(j=2;j<=i;j++){
      s[i-1][j-1] = s[i-1][j-2]+(s[i-1][j-2]-s[i-2][j-2])/(pow(4.,j-1)-1);
    }
  }
  
  norm = s[Nint-1][Nint-1];

  return res / norm;

}

double PDF_halo(double velocity, double r, double mass1, double mass2){

  double sigm = sigma_NL_mass;
  double sig1 = log(1+sigm);
  double rd1 = r_delta(zhalo, mass1, 200) * (1.0+zhalo);

  double res,abserr;
  //double x_lo = log(1e-3);
  //double x_hi = log(1e+3);
  double x_lo = -10.0 * sqrt(sig1);
  double x_hi = +10.0 * sqrt(sig1);
  
  int nint=100;
  double dx=(x_hi-x_lo)/nint;
  double tab1[nint+1], tab2[nint+1];

  for(int i=0;i<=nint;i++){
    double x = x_lo + (i+0.0)*dx;
    tab1[i] =  PDF_halo_int_d(velocity, x, r, sigm, mass1, mass2);
    tab2[i] = PDF_halo_int_n(x, r, sigm, mass1, mass2);
  }

  res=0;      
  for(int i=0; i<nint/2; i++){
    double xl = x_lo + (2*i+0)*dx;
    double xm = x_lo + (2*i+1)*dx;
    double xh = x_lo + (2*i+2)*dx;
    
    double yl, ym, yh;
    
    yl = tab1[2*i+0];
    ym = tab1[2*i+1];
    yh = tab1[2*i+2];

    res += dx * (yl+4*ym+yh) / 3;
  }    

  double norm =0;
  for(int i=0; i<nint/2; i++){
    double xl = x_lo + (2*i+0)*dx;
    double xm = x_lo + (2*i+1)*dx;
    double xh = x_lo + (2*i+2)*dx;
    
    double yl, ym, yh;
    
    yl = tab2[2*i+0];
    ym = tab2[2*i+1];
    yh = tab2[2*i+2];

    norm += dx * (yl+4*ym+yh) / 3;
  }    

  return res/norm;
  
}

double modified_ln_pdf(double x, double *p, int opt){
  double N = p[0];
  double A = p[1];
  double w2 = p[2];

  if(opt==-1){
    return N * exp(-pow(x+w2*0.5,2)/2/w2*(1.0+A*exp(-x)));
  }
  if(opt==0){
    return exp(-pow(x+w2*0.5,2)/2/w2*(1.0+A*exp(-x)));
  }
  if(opt==1){
    return N * exp(-pow(x+w2*0.5,2)/2/w2*(1.0+A*exp(-x))) * (-pow(x+w2*0.5,2)/2/w2) * exp(-x);
  }
  if(opt==2){
    return N * exp(-pow(x+w2*0.5,2)/2/w2*(1.0+A*exp(-x))) * (1.0+A*exp(-x)) * (-(x+w2*0.5)/2/w2) * (1.0+1/w2*(x+w2*0.5));
  }
}

double modified_ln_pdf_norm_int(double x, double *p, int opt){
  return modified_ln_pdf(x, p, opt);
}

double modified_ln_pdf_norm(double *p, int opt){

  double w2 = p[2];

  double res,abserr;
  double x_lo = -10.0 * sqrt(w2);
  double x_hi = +10.0 * sqrt(w2);

  int Nint = 10;
  //Romberg
  int i,j;
  double h[Nint];
  double s[Nint][Nint];
  for(i=1;i<=Nint;i++){h[i-1] = (x_hi-x_lo)/pow(2.,i-1);}

  s[0][0] = 0.5*h[0]*(modified_ln_pdf_norm_int(x_hi,p,opt)+modified_ln_pdf_norm_int(x_lo, p,opt));
  for(i=2;i<=Nint;i++){
    s[i-1][0] = s[i-2][0];
    for(j=1;j<=pow(2.,i-2);j++){
      s[i-1][0] += h[i-2]*modified_ln_pdf_norm_int(x_lo+(2*j-1)*h[i-1], p, opt);
    }
    s[i-1][0] = 0.5*s[i-1][0];
  }

  for(i=2;i<=Nint;i++){
    for(j=2;j<=i;j++){
      s[i-1][j-1] = s[i-1][j-2]+(s[i-1][j-2]-s[i-2][j-2])/(pow(4.,j-1)-1);
    }
  }

  res = s[Nint-1][Nint-1];

  return res;
}

double modified_ln_pdf_mean_int(double x, double *p, int opt){
  return exp(x) * modified_ln_pdf(x, p, opt);
}

double modified_ln_pdf_mean(double *p, int opt){
    double w2 = p[2];

  double res,abserr;
  double x_lo = -10.0 * sqrt(w2);
  double x_hi = +10.0 * sqrt(w2);

  int Nint = 10;
  //Romberg
  int i,j;
  double h[Nint];
  double s[Nint][Nint];
  for(i=1;i<=Nint;i++){h[i-1] = (x_hi-x_lo)/pow(2.,i-1);}

  s[0][0] = 0.5*h[0]*(modified_ln_pdf_mean_int(x_hi,p,opt)+modified_ln_pdf_mean_int(x_lo, p, opt));
  for(i=2;i<=Nint;i++){
    s[i-1][0] = s[i-2][0];
    for(j=1;j<=pow(2.,i-2);j++){
      s[i-1][0] += h[i-2]*modified_ln_pdf_mean_int(x_lo+(2*j-1)*h[i-1], p, opt);
    }
    s[i-1][0] = 0.5*s[i-1][0];
  }

  for(i=2;i<=Nint;i++){
    for(j=2;j<=i;j++){
      s[i-1][j-1] = s[i-1][j-2]+(s[i-1][j-2]-s[i-2][j-2])/(pow(4.,j-1)-1);
    }
  }

  res = s[Nint-1][Nint-1];

  return res;
}

double modified_ln_pdf_var_int(double x, double *p, int opt){
  return (exp(x)-1) * (exp(x)-1) * modified_ln_pdf(x, p, opt);
}

double modified_ln_pdf_var(double *p, int opt){
  double w2 = p[2];

  double res,abserr;
  double x_lo = -10.0 * sqrt(w2);
  double x_hi = +10.0 * sqrt(w2);

  int Nint = 10;
  //Romberg
  int i,j;
  double h[Nint];
  double s[Nint][Nint];
  for(i=1;i<=Nint;i++){h[i-1] = (x_hi-x_lo)/pow(2.,i-1);}

  s[0][0] = 0.5*h[0]*(modified_ln_pdf_var_int(x_hi,p,opt)+modified_ln_pdf_var_int(x_lo, p, opt));
  for(i=2;i<=Nint;i++){
    s[i-1][0] = s[i-2][0];
    for(j=1;j<=pow(2.,i-2);j++){
      s[i-1][0] += h[i-2]*modified_ln_pdf_var_int(x_lo+(2*j-1)*h[i-1], p, opt);
    }
    s[i-1][0] = 0.5*s[i-1][0];
  }

  for(i=2;i<=Nint;i++){
    for(j=2;j<=i;j++){
      s[i-1][j-1] = s[i-1][j-2]+(s[i-1][j-2]-s[i-2][j-2])/(pow(4.,j-1)-1);
    }
  }

  res = s[Nint-1][Nint-1];

  return res;
}

int modified_ln_mod_constraints(const gsl_vector * x, void *p, gsl_vector * f){

  double *pin = (double *)p;
  double s2 = pin[0];
  
  double x0 = gsl_vector_get(x, 0);
  double x1 = gsl_vector_get(x, 1);
  double x2 = gsl_vector_get(x, 2);

  double params[3] = {exp(x0), exp(x1), exp(x2)};

  gsl_vector_set (f, 0, modified_ln_pdf_norm(params, -1)-1.0);
  gsl_vector_set (f, 1, modified_ln_pdf_mean(params, -1)-1.0);
  gsl_vector_set (f, 2, modified_ln_pdf_var(params, -1)-s2);
  
  return GSL_SUCCESS;

}

int modified_ln_mod_constraints_df(const gsl_vector * x, void *p, gsl_matrix *J){

  double *pin = (double *)p;
  double s2 = pin[0];
  
  double x0 = gsl_vector_get(x, 0);
  double x1 = gsl_vector_get(x, 1);
  double x2 = gsl_vector_get(x, 2);

  double params[3] = {exp(x0), exp(x1), exp(x2)};

  gsl_matrix_set(J, 0, 0, modified_ln_pdf_norm(params, 0)*exp(x0));
  gsl_matrix_set(J, 0, 1, modified_ln_pdf_norm(params, 1)*exp(x1));
  gsl_matrix_set(J, 0, 2, modified_ln_pdf_norm(params, 2)*exp(x2));

  gsl_matrix_set(J, 1, 0, modified_ln_pdf_mean(params, 0)*exp(x0));
  gsl_matrix_set(J, 1, 1, modified_ln_pdf_mean(params, 1)*exp(x1));
  gsl_matrix_set(J, 1, 2, modified_ln_pdf_mean(params, 2)*exp(x2));

  gsl_matrix_set(J, 2, 0, modified_ln_pdf_var(params, 0)*exp(x0));
  gsl_matrix_set(J, 2, 1, modified_ln_pdf_var(params, 1)*exp(x1));
  gsl_matrix_set(J, 2, 2, modified_ln_pdf_var(params, 2)*exp(x2));
  
  return GSL_SUCCESS;
  
}

int modified_ln_mod_constraints_fdf(const gsl_vector * x, void *p, gsl_vector * f, gsl_matrix *J){

  modified_ln_mod_constraints(x,p,f);
  modified_ln_mod_constraints_df(x,p,J);
  
  return GSL_SUCCESS;
}

int solve_modified_ln (double s2, double tolerance, double *pinit, double *out){

  const gsl_multiroot_fdfsolver_type *T;
  gsl_multiroot_fdfsolver *s;
  
  gsl_vector *v;
  size_t i, iter = 0, max_iter=1000;
  const size_t ndim = 3;

  int status;
  double size;

  double p[3] = {s2, 0, 0};
  
  gsl_multiroot_function_fdf func = {&modified_ln_mod_constraints,
                                     &modified_ln_mod_constraints_df,
                                     &modified_ln_mod_constraints_fdf, ndim, &p};

  /* Starting point */
  v = gsl_vector_alloc (ndim);
  gsl_vector_set (v, 0, log(pinit[0]));
  gsl_vector_set (v, 1, log(pinit[1]));
  gsl_vector_set (v, 2, log(pinit[2]));
  
  /* Initialize method and iterate */
  T = gsl_multiroot_fdfsolver_hybridsj;
  s = gsl_multiroot_fdfsolver_alloc (T,ndim);
  gsl_multiroot_fdfsolver_set (s, &func, v);

  do {
    iter++;
    status = gsl_multiroot_fdfsolver_iterate(s);
    if (status) break;
    status = gsl_multiroot_test_residual (s->f, tolerance);
    //if ((status == GSL_SUCCESS)) {
      //printf ("converged to minimum\n");
      //printf ("Niter= %5d, lnN= %e, lnA= %e, lnw2 = %e, cond1=%e, cond2=%e, cond3=%e\n",
      //      (int)iter,
      //      gsl_vector_get(s->x, 0), gsl_vector_get(s->x,1), gsl_vector_get(s->x, 2),
      //      gsl_vector_get(s->f, 0), gsl_vector_get(s->f,1), gsl_vector_get(s->f, 2)); 
    //}

  }while(status == GSL_CONTINUE && iter < max_iter);
  
  double u0 = gsl_vector_get (s->x, 0);
  double u1 = gsl_vector_get (s->x, 1);
  double u2 = gsl_vector_get (s->x, 2);

  out[0] = exp(u0);
  out[1] = exp(u1);
  out[2] = exp(u2);

  gsl_vector_free(v);
  gsl_multiroot_fdfsolver_free(s);
  
  return status;

}

void get_modified_ln_params(double s2, double *out){

  int i;
  double lk;
  double pow_index;

  lk=log10(s2);

  if(s2 > pow(10., tab_mod_ln_s2[0]) || s2 < pow(10., tab_mod_ln_s2[NPOINTS-1])){
    fprintf(stderr, "can not get the parameters for modified log normal at s2 = %e\n", s2);
    exit(1);
  }

  splint(tab_mod_ln_s2-1, tab_mod_ln_N-1, err_mod_ln_N-1, NPOINTS, lk, &pow_index);
  out[0] = pow(10.0,pow_index);

  splint(tab_mod_ln_s2-1, tab_mod_ln_A-1, err_mod_ln_A-1, NPOINTS, lk, &pow_index);
  out[1] = pow(10.0,pow_index);

  splint(tab_mod_ln_s2-1, tab_mod_ln_w2-1, err_mod_ln_w2-1, NPOINTS, lk, &pow_index);
  out[2] = pow(10.0,pow_index);
  
}

void set_modified_ln_params(){
  int i;
  double pinit[3] = {0.95, 4.0, 2.0};
  double p_mod_ln[3];
  double logs2_max = +1.0;
  double logs2_min = -2.0;
  double dlogs2 = (logs2_min-logs2_max)/(NPOINTS-1);

  for(i=0;i<NPOINTS;i++){
    double logs2 = logs2_max + dlogs2 * i;
    double s2 = pow(10.0, logs2);
    solve_modified_ln (s2, 1e-5, pinit, p_mod_ln);
    pinit[0] = p_mod_ln[0];
    pinit[1] = p_mod_ln[1];
    pinit[2] = p_mod_ln[2];
    //fprintf(stdout, "%e %e %e %e\n", s2, p_mod_ln[0], p_mod_ln[1], p_mod_ln[2]);
    tab_mod_ln_s2[i] = log10(s2);
    tab_mod_ln_N[i] = log10(p_mod_ln[0]);
    tab_mod_ln_A[i] = log10(p_mod_ln[1]);
    tab_mod_ln_w2[i] = log10(p_mod_ln[2]);
  }

  double yp1 = 1.e31;
  double ypn = 1.e31;
  spline(tab_mod_ln_s2-1, tab_mod_ln_N-1, NPOINTS, yp1, ypn, err_mod_ln_N-1);
  spline(tab_mod_ln_s2-1, tab_mod_ln_A-1, NPOINTS, yp1, ypn, err_mod_ln_A-1);
  spline(tab_mod_ln_s2-1, tab_mod_ln_w2-1, NPOINTS, yp1, ypn, err_mod_ln_w2-1);
    
}
