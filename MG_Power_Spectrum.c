#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#define PI 3.141592
#define Nphi 500
#define Ndel 200

/*Structure with the information about the cosmolgy*/
typedef struct Cosmology {
	float H0;
	float w;
	float Ob;
	float Odm;
	float Om;
	float Ol;
	float Ok;
	float Onu;
	float A;
	float ns;
	float Yhe;
	float T;
} COSMOLOGY;

/*Structure with the information about Hu-Sawicki theory parameters*/
typedef struct HuSawicki {
	float n;
	float fr0;
} HS;

/*Structure with the information about Symmetron theory parameters*/
typedef struct Symmetron {
	float L;
	float zssb;
	float beta;
} SYM;

/*Global variables contaning the information about the cosmoloy and modified gravity*/
int model;
float *u;
COSMOLOGY cosmo;
HS hs;
SYM sym;

/*The Hubble parameters over H0 (E(a))*/
float E(float a){
	float resp;
	resp = sqrt(cosmo.Om*pow(a,-3) + cosmo.Ok*pow(a,-2) + cosmo.Ol*pow(a,-3.0*(1.0 + cosmo.w)));
	return resp;
}

/*Derivative of equation above*/
float dE(float a){
	float resp;
	resp = -0.5*(3.0*cosmo.Om*pow(a,-4) + 2.0*cosmo.Ok*pow(a,-3))/E(a);
	return resp;
}

/*Equation for \dot{v} = \ddot{u}*/
float ddu(float a, float u, float v, float mu2, float assb){
	float resp;
	resp = -(4.0/a + dE(a)/E(a))*v - mu2/(a*a*E(a)*E(a))*((pow(assb/a,3) - 1.0)*u + pow(u,3));
	return resp;
}

/*Equation for \dot{u}*/
float du(float a, float u, float v){
	return v;
}

/*Solve the equation for the symmetron field*/
void phi(float aini, float assb, float mu2, float u[], int N){
	float *v, h, a;
	float q1, q2, q3, q4, k1, k2, k3, k4;
	int i, j;
	
	/*Allocate the velocity vector*/
	v = (float*)malloc(N*sizeof(float));	

	/*Initial conditions*/
	u[0] = 1e-3;	
	v[0] = 1e-3;

	/*Solve the EDO*/
	h = (1.0 - aini)/(N-1);

	/*Solve the EDO for each k using Runge-Kutta 4*/
	a = aini;
	for(i=0;i<N-1;i++){
		q1 = du(a, u[i], v[i]);
		k1 = ddu(a, u[i], v[i], mu2, assb);
		q2 = du(a+h/2.0, u[i]+h/2.0*q1, v[i]+h/2.0*k1);
		k2 = ddu(a+h/2.0, u[i]+h/2.0*q1, v[i]+h/2.0*k1, mu2, assb);
		q3 = du(a+h/2.0, u[i]+h/2.0*q2, v[i]+h/2.0*k2);
		k3 = ddu(a+h/2.0, u[i]+h/2.0*q2, v[i]+h/2.0*k2, mu2, assb);
		q4 = du(a+h, u[i]+h*q3, v[i]+h*k3);		
		k4 = ddu(a+h, u[i]+h*q3, v[i]+h*k3, mu2, assb);

		u[i+1] = u[i] + h/6.0*(q1 + 2.0*q2 + 2.0*q3 + q4);	
		v[i+1] = v[i] + h/6.0*(k1 + 2.0*k2 + 2.0*k3 + k4);
		a = a + h;
	}
}

/*Mass of the symmetron field*/
float Mass_sym(float a, float assb, float mu2){
	float m2;	
	if(a < assb)
		m2 = mu2*(pow(assb/a, 3.0) - 1.0);
	else
		m2 = 2.0*mu2*(1.0 - pow(assb/a, 3.0));
	return m2;
}

/*Wave lenght of the symmetron field*/
float Lamb_sym(float a, float assb){
	float lamb2;	
	if(a < assb)
		lamb2 = 2.0*pow(sym.L,2)/(pow(assb/a, 3.0) - 1.0);
	else
		lamb2 = pow(sym.L,2)/(1.0 - pow(assb/a, 3.0));
	return lamb2;
}

/*The value of the potential minimum*/
float min(float a, float assb){
	float resp;
	if(a<assb)
		resp = 0.0;
	else
		resp = sqrt(1.0 - pow(assb/a, 3.0));	
	return resp;
}

/*Function that sumarize the modied Poisson's equation in f(R)*/
float Mu(float k, float a){
	float resp;
	if(model == 1){
		float m0, m;
		m0 = (100.0/299792.0)*sqrt((cosmo.Om + 4.0*cosmo.Ol)/((hs.n + 1.0)*hs.fr0));
		m = m0*pow((cosmo.Om*pow(a,-3)+4.0*cosmo.Ol)/(cosmo.Om + 4.0*cosmo.Ol),(hs.n+2.0)/2.0);
		resp = ((4.0/3.0)*pow(k,2) + pow(m*a,2))/(pow(k,2) + pow(m*a,2));
	}
	if(model == 2){
		float lamb2, assb;
		int i;
		assb = 1.0/(1.0 + sym.zssb);
		lamb2 = Lamb_sym(a,assb);
		i = (int)floor((a - 0.01)/(1.0 - 0.01)*(Nphi*Ndel - 1.0));
		/*resp = 1.0 + 2.0*pow(min(a,assb)*sym.beta,2)/(1.0 + pow(a/k,2)/lamb2);*/
		if(a<=assb)	resp = 1.0;
		else resp = 1.0 + (2.0*pow(sym.beta*sym.L*k, 2.0)*(1.0 - pow(assb/a,3.0)))/(pow(k*sym.L,2.0) + a*a - pow(assb,3)/a);
	}
	return resp;
}

/*Equation for \dot{v} = \ddot{\delta}*/
float dv(float k, float a, float d, float v){
	float resp;
	resp = -(3.0/a + dE(a)/E(a))*v + 3.0/2.0*cosmo.Om*Mu(k, a)/(pow(E(a),2)*pow(a,5))*d;
	return resp;
}

/*Equation for \dot{\delta}*/
float dd(float a, float d, float v){
	return v;
}

/*Create the param file for input in CAMB*/
void Create_Params_ini(float Rmin, float zini){
	FILE *camb;
	char *cambfile = "params_mg.ini";
	float kmax;

	kmax = 2.0*PI/Rmin;

	camb = fopen(cambfile, "w");
	if (camb == NULL) {
		printf("Unable to open %s\n",cambfile);
		exit(0);
	}

	fprintf(camb, "output_root = test\n"
		"get_scalar_cls = F\nget_vector_cls = F\nget_tensor_cls = F\nget_transfer   = T\n"
		"do_lensing     = T\n"
		"do_nonlinear = 0\n"
		"l_max_scalar      = 2200\nk_eta_max_scalar  = 4000\n"
		"l_max_tensor      = 1500\nk_eta_max_tensor  = 3000\n"
		"use_physical   = F\n#ombh2          = 0.0226\n#omch2          = 0.112\n#omnuh2         = 0\nomk            = %f\nhubble         = %f\n"
		"w              = %f\ncs2_lam        = 1\n"
		"omega_baryon   = %f\nomega_cdm      = %f\nomega_lambda   = %f\nomega_neutrino = %f\n"
		"temp_cmb           = %f\nhelium_fraction    = %f\n"
		"massless_neutrinos = 3.046\nmassive_neutrinos  = 0\n"
		"nu_mass_eigenstates = 0\nnu_mass_degeneracies = 0\nnu_mass_fractions = 1\n"
		"initial_power_num         = 1\npivot_scalar              = 0.05\npivot_tensor              = 0.05\nscalar_amp(1)             = %e\nscalar_spectral_index(1)  = %f\nscalar_nrun(1)            = 0\ntensor_spectral_index(1)  = 0\ninitial_ratio(1)          = 1\n"
		"reionization         = T\n"
		"re_use_optical_depth = T\nre_optical_depth     = 0.09\nre_redshift          = 11\nre_delta_redshift    = 1.5\nre_ionization_frac   = -1\n"
		"RECFAST_fudge = 1.14\nRECFAST_fudge_He = 0.86\nRECFAST_Heswitch = 6\nRECFAST_Hswitch  = T\n"
		"initial_condition   = 1\ninitial_vector = -1 0 0 0 0\n"
		"vector_mode = 0\n"
		"COBE_normalize = F\nCMB_outputscale = 7.42835025e12\n"
		"transfer_high_precision = F\ntransfer_kmax           = %f\ntransfer_k_per_logint   = 0\ntransfer_num_redshifts  = 2\ntransfer_interp_matterpower = T\ntransfer_redshift(1)    = %f\ntransfer_redshift(2)    = %f\ntransfer_filename(1)    = transfer_out100.dat\ntransfer_filename(2)    = transfer_out99.dat\ntransfer_matterpower(1) = matterpower100.dat\ntransfer_matterpower(2) = matterpower99.dat\n"
		"scalar_output_file = scalCls.dat\nvector_output_file = vecCls.dat\ntensor_output_file = tensCls.dat\ntotal_output_file  = totCls.dat\nlensed_output_file = lensedCls.dat\nlensed_total_output_file  =lensedtotCls.dat\nlens_potential_output_file = lenspote	ntialCls.dat\nFITS_filename      = scalCls.fits\n"
		"do_lensing_bispectrum = F\ndo_primordial_bispectrum = F\n"
		"bispectrum_nfields = 1\nbispectrum_slice_base_L = 0\nbispectrum_ndelta=3\nbispectrum_delta(1)=0\nbispectrum_delta(2)=2\nbispectrum_delta(3)=4\nbispectrum_do_fisher= F\nbispectrum_fisher_noise=0\nbispectrum_fisher_noise_pol=0\nbispectrum_fisher_fwhm_arcmin=7\nbispec	trum_full_output_file=\nbispectrum_full_output_sparse=F\nbispectrum_export_alpha_beta=F\n"
		"feedback_level = 1\n"
		"derived_parameters = F\n"
		"lensing_method = 1\naccurate_BB = F\n"
		"massive_nu_approx = 1\n"
		"accurate_polarization   = T\n"
		"accurate_reionization   = T\n"
		"do_tensor_neutrinos     = T\n"
		"do_late_rad_truncation   = T\n"
		"number_of_threads       = 0\n"
		"high_accuracy_default=T\n"
		"accuracy_boost          = 1\n"
		"l_accuracy_boost        = 1\n"
		"l_sample_boost          = 1\n", cosmo.Ok, cosmo.H0, cosmo.w, cosmo.Ob, cosmo.Odm, cosmo.Ol, cosmo.Onu, cosmo.T, cosmo.Yhe, cosmo.A, cosmo.ns, kmax, zini + 1.0, zini);

	fclose(camb);
}


/*Define the window function for the sigma's calculation*/
float W(float k, float R){
	float resp;
	resp = 3.0/(pow(k*R,2))*(sin(k*R)/(k*R) - cos(k*R));
	return resp;
}

/*Evaluate the square root of matter variance*/
float calc_sigma(float *k, float *P, int cont, float R){
	int i;
	float resp;
	resp = 0.0;
	for(i=0;i<cont-2;i++)
		resp += (k[i+1] - k[i])/2.0*(P[i]*pow(k[i],2)*pow(W(k[i],R),2) + P[i+1]*pow(k[i+1],2)*pow(W(k[i+1],R),2));
	return resp/(2.0*PI*PI);
}

int main(int argc,char *argv[])
{
FILE *param, *camb;
char paramfile[20];
float Rmin, Rmax, n, fr0, zini;
float k[1000], Pini[1000], Paux[1000], d[1000], v[1000], P[1000], aini, aaux, h, a;
float q1, q2, q3, q4, k1, k2, k3, k4, trash;
float sigma[100], R[100];
int cont, i, j, N;

if (argc != 2){
	printf("Wrong number of arguments.\n");
	printf("arg1: The name of parameters file.\n\n");	
	exit(0);
}

sprintf(paramfile,"%s", argv[1]);

/*Open the oparameters file*/
param = fopen(paramfile, "r");
if (param == NULL) {
	printf("Unable to open %s\n",paramfile);
	exit(0);
}

/*Read the parameters*/
fscanf(param,"%d", &model);
fscanf(param,"%f %f", &Rmin, &Rmax);
if(model == 1)
	fscanf(param,"%f %f %f", &hs.n, &hs.fr0, &trash);
if(model == 2)
	fscanf(param,"%f %f %f", &sym.zssb, &sym.beta, &sym.L);
fscanf(param,"%f %f %f %f %f %f %f %f %f %f %f", &cosmo.H0, &cosmo.w, &cosmo.Odm, &cosmo.Ob, &cosmo.Ol, &cosmo.Ok, &cosmo.Onu, &cosmo.Yhe, &cosmo.A, &cosmo.ns, &cosmo.T);
cosmo.Om = cosmo.Ob + cosmo.Odm;
zini = 99.0;

fclose(param);

printf("%d %f %f %f\n", model, sym.zssb, sym.beta, sym.L);

/*Create the params.ini file for camb*/
Create_Params_ini(Rmin, zini);

/*Run the camb in the zini*/
system("~/Documentos/CAMB/camb/camb params_mg.ini");

/*Read the CAMB output (z=99)*/
camb = fopen("test_matterpower99.dat", "r");
if (camb == NULL) {
	printf("Unable to open %s\n", "test_matterpower99.dat");
	exit(0);
}

cont = 0;
while(!feof(camb)){
	fscanf(camb,"%f %f", &k[cont], &Pini[cont]);
	cont ++;
}
fclose(camb);

/*Read the aux power spectrum (z=100)*/
camb = fopen("test_matterpower100.dat", "r");
if (camb == NULL) {
	printf("Unable to open %s\n", "test_matterpower100.dat");
	exit(0);
}

cont = 0;
while(!feof(camb)){
	fscanf(camb,"%f %f", &k[cont], &Paux[cont]);
	cont ++;
}
fclose(camb);

/*Delete the crated files*/
system("rm test_params.ini test_transfer_out99.dat test_transfer_out100.dat test_matterpower99.dat test_matterpower100.dat");

/*Construct the initial condictions*/
aini = 1.0/(1.0 + zini);
aaux = 1.0/(2.0 + zini);
for(i=0;i<cont;i++){
	d[i] = sqrt(Pini[i]);
	v[i] = (d[i] - sqrt(Paux[i]))/(aini - aaux);
}

/*Solve the EDO*/
N = Ndel;	/*number of steps*/
h = (1.0 - aini)/(N-1);

/*Evaluate the background symmetron evolution
if(model == 2){
	u = (float *)malloc((Nphi*N)*sizeof(float));
	phi(aini, 1.0/(1.0 + sym.zssb), pow(2997.92/(sqrt(2.0)*sym.L), 2.0), u, Nphi*N);
}

for(i=0;i<Nphi*N;i++)
	printf("%f\n", u[i]);
*/
omp_set_num_threads(6);
#pragma omp parallel for private(i, q1, q2, q3, q4, k1, k2, k3, k4, a)
for(j=0;j<cont;j++){
/*Solve the EDO for each k using Runge-Kutta 4*/
	a = aini;
	for(i=0;i<N;i++){
		q1 = dd(a, d[j], v[j]);
		k1 = dv(k[j], a, d[j], v[j]);
		q2 = dd(a+h/2.0, d[j]+h/2.0*q1, v[j]+h/2.0*k1);
		k2 = dv(k[j], a+h/2.0, d[j]+h/2.0*q1, v[j]+h/2.0*k1);
		q3 = dd(a+h/2.0, d[j]+h/2.0*q2, v[j]+h/2.0*k2);
		k3 = dv(k[j], a+h/2.0, d[j]+h/2.0*q2, v[j]+h/2.0*k2);
		q4 = dd(a+h, d[j]+h*q3, v[j]+h*k3);		
		k4 = dv(k[j], a+h, d[j]+h*q3, v[j]+h*k3);

		d[j] = d[j] + h/6.0*(q1 + 2.0*q2 + 2.0*q3 + q4);	
		v[j] = v[j] + h/6.0*(k1 + 2.0*k2 + 2.0*k3 + k4);
		a = a + h;
	}
}

for(i=0;i<cont;i++)
	P[i] = d[i]*d[i];

camb = fopen("Matter_Power_MG.dat", "w");
if (camb == NULL) {
	printf("Unable to open %s\n", "Matter_Power_MG.dat");
	exit(0);
}

for(i=0;i<cont-1;i++)
	fprintf(camb,"%f %f\n", k[i], P[i]);
fclose(camb);

/*Evaluating the radius array in log10 space*/ 
h = (log10(Rmax) - log10(Rmin))/99.0;
for(i=0;i<100;i++)
	R[i] = pow(10, log10(Rmin) + i*h);

/*Evaluating the square root of variance (\sigma)*/
omp_set_num_threads(6);
#pragma omp parallel for private(i)
for(i=0;i<100;i++)
	sigma[i] = sqrt(calc_sigma(k, P, cont, R[i]));

/*Print the variance*/
camb = fopen("Sigma_MG.dat", "w");
if (camb == NULL) {
	printf("Unable to open %s\n", "Sigma_MG.dat");
	exit(0);
}

for(i=0;i<100;i++)
	fprintf(camb,"%f %f\n", R[i], sigma[i]);
fclose(camb);

/*Print the seigma_8*/
printf("Sigma_8 = %f\n", sqrt(calc_sigma(k, P, cont, 8.0)));

return 0;
}
