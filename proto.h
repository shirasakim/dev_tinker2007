#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

double PowerSpec(double k);
double GrowthFactor(double astart, double aend);
double growth(double a);
inline double growth_int(double);
double qromb(double (*func)(double), double a, double b);
double sigma2_int(double k);
double xi_int(double k);
inline double TopHatSigma2(double R);
inline double PowerSpec_Efstathiou(double k);
inline double PowerSpec_EH(double k);
inline double PowerSpec_BBKS(double k);
inline double PowerSpec_CMBFAST(double k);

int initialize_powerspectrum(int Spectrum);
int set_units();

double tk_eh(double k);
double transfunc_cmbfast(double k);
double transfunc_WDM(double k);
inline double F_Omega(double a);
inline double Hubble_a(double a);

inline double delta_c_func();
inline double efnn(double x, double);
inline double var4(double x, double);
inline double efn(double x, double);
inline double var3(double x, double);
inline double ddweight(double x, double);
inline double dweight(double x, double);
inline double weight(double x, double);
double dlogdsigma(double mass, double, double);

double sigma_m(double m, double *rsphere_return);
double fract_mass(double sig);
double sigdsigdr(double);
double dsig2dr2(double);
double unnsigma(double);
inline double evar2(double x, double);
inline double var2(double x, double);
double dndlogm(double logm);

void readCMB_and_do_spline();

void spline(double x[], double y[], int n, double yp1, double ypn, double y2[]);
void splint(double xa[], double ya[], double y2a[], int n, double x, double *y);

void splie2(double x1a[], double x2a[], double **ya, int m, int n, double **y2a);
void splin2(double x1a[], double x2a[], double **ya, double **y2a, int m, int n, double x1, double x2, double *y);
void splie3(double x1a[], double x2a[], double x3a[], double ***ya, int l, int m, int n, double ***y2a);
void splin3(double x1a[], double x2a[], double x3a[], double ***ya, double ***y2a, int l, int m, int n, double x1, double x2, double x3, double *y);

double H_z(double z);
double H_z_1(double z,void *params);
double chi(double z);
double r_vir(double z, double M);
double delta_c(double z,double M);
double delta_v(double z);
double r_delta(double z, double Mass, double delta);
double get_eq_for_M_vir(double x, void *params);
double M_delta_to_M_vir(double z, double mdelta, double delta);

double get_eq_for_M_delta(double x, void *params);
double M_vir_to_M_delta(double z, double mvir, double delta);

double halo_bias(double logm);
double halo_bias_fast(double logm, double z);

double Omega_de(double a);
double coeff1(double a);
double coeff2(double a);//u"+coeff1(a)u'+coeff2(a)u=0, u=D/a, '=d/dlna
double RungeKutta(double a_in,double a); //calculate linear density growth eq.

void stack_table_and_spline();

double chi_fast(double z);
void growth_spline();
double growth_for_any_w(double a);
double dndlogm_fast(double logm, double z);
double Redshift(double x,double z_in);

//non-linear matter Pk
double P_nonlinear(double z, double k);
void set_halofit_param(double z, double *param);
double solver(double z);
double get_delta_k(double k);
double sigma2_gauss_int(double lnk, void *params);
double sigma2_gauss(double lnR, void *params);
double dsig_dR_int(double lnk, void *params);
double dsig_dR(double R);
double d2sig_dR2_int(double lnk, void *params);
double d2sig_dR2(double R);
double neff(double R);
double C_halofit(double R);

void stack_data_and_spline_Pk();

double dsig_dR_fast(double R);
double d2sig_dR2_fast(double R);
double sigma2_gauss_fast(double lnR);

//galaxy HOD
double HOD_gal(double M, int cen_or_sat);
double mean_ngal_int(double logM, double z);
double mean_ngal(double z);
double mean_ngal_fast(double z);

double dens_sat_int(double logM, double z);
double dens_sat(double z);

//halo concentration by Diemer & Kravtosov 2015, https://arxiv.org/abs/1407.4730
double dlnP_dlnk(double lnk);
double nu_M(double z, double M);
double nu_M_fast(double z, double M);
double c_200c_DK15(double z, double M);
double c_200c_DK15_fast(double z, double M);
double get_eq_for_M_200c(double x, void *params);
double M_delta_to_M_200c(double z, double mdelta, double delta);

//two-halo term in tinker 2007 http://articles.adsabs.harvard.edu/pdf/2007MNRAS.374..477T
void save_SC_param();
double SC_param_density(double gamma);
double SC_param_velocity(double gamma);
double SC_param_velocity_fast(double density);

double mean_infall_velocity(double r, double density);
double sigma_r_velocity(double r, double density, double mass1, double mass2);
double TopHatSigma2_NL(double R);
double sigma2_NL_int(double k, double z);
double PDF_mass_NL(double density, double r, double sigm, double mass1, double mass2);
double PDF_halo_int_n(double density, double r, double sigm, double mass1, double mass2);
double mean_overdesnity_ln(double r, double sigm, double mass1, double mass2);
double PDF_halo_int_d(double velocity, double density, double r, double sigm, double mass1, double mass2);
double PDF_halo(double velocity, double r, double mass1, double mass2);

//modified log normal proposed in https://iopscience.iop.org/article/10.1086/504032/pdf
void set_modified_ln_params();
double modified_ln_pdf(double x, double *p);

double modified_ln_pdf_norm_int(double x, double *p, int opt);
double modified_ln_pdf_norm(double *p, int opt);
double modified_ln_pdf_mean_int(double x, double *p, int opt);
double modified_ln_pdf_mean(double *p, int opt);
double modified_ln_pdf_var_int(double x, double *p, int opt);
double modified_ln_pdf_var(double *p, int opt);

int modified_ln_mod_constraints(const gsl_vector * x, void *p, gsl_vector * f);
int modified_ln_mod_constraints_df(const gsl_vector * x, void *p, gsl_matrix * J);
int modified_ln_mod_constraints_fdf(const gsl_vector * x, void *p, gsl_vector * f, gsl_matrix * J);

int solve_modified_ln (double s2, double tolerance, double *pinit, double *out);
void get_modified_ln_params(double s2, double *out);
