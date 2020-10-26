#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define PI 3.141592

/*Structure with the information about the cosmolgy*/
typedef struct Cosmology {
	double H0;
	double w;
	double Om;
	double Ol;
	double Ok;
} COSMOLOGY;

/*Structure with the information about Hu-Sawicki theory parameters*/
typedef struct HuSawicki {
	double n;
	double fr0;
} HS;

/*Global variables contaning the information about the cosmoloy and modified gravity*/
COSMOLOGY cosmo;
HS hs;

/*The Hubble parameters over H0 (E(a))*/
double E(double n){
	double resp;
	resp = sqrt(cosmo.Om*exp(-3.0*n) + cosmo.Ok*exp(-2.0*n) + cosmo.Ol*exp(-3.0*(1.0 + cosmo.w)*n));
	return resp;
}

/*Derivative of equation above*/
double dE(double n){
	double resp;
	resp = -0.5*(3.0*cosmo.Om*exp(-3.0*n) + 2.0*cosmo.Ok*exp(-2.0*n))/E(n);
	return resp;
}

/*Define the mass of the scalar field in the background*/
double Mass(double n){
	double resp, m0;

	m0 = cosmo.H0*sqrt((cosmo.Om + 4.0*cosmo.Ol)/((hs.n + 1.0)*hs.fr0));
	resp = m0*pow((cosmo.Om*exp(-3.0*n)+4.0*cosmo.Ol)/(cosmo.Om + 4.0*cosmo.Ol),(hs.n+2.0)/2.0);

	return resp;
}

/*Function that defines the change in gravity at linear order*/
double Eps(double k, double n){
	double resp, b, m;

	b = 1.0/sqrt(6);
	m = Mass(n);

	resp = 2.0*pow(b, 2.0)*pow(k, 2.0)/(pow(k, 2.0) + pow(m, 2.0)*exp(2.0*n));

	return resp;
}


/*Integrate using the Simpson's method*/
double Integrate(double *f, double *x, int N, double h){
	double resp;
	int i;

	resp = f[0];

	for(i=1;i<N/2;i++)
		resp += 2.0*f[2*i];

	for(i=1;i<N/2 + 1;i++)
		resp += 4.0*f[2*i - 1];

	resp += f[N];

	return h*resp/3.0;
}

/*Integrate using the Romberg's method*/
double Romb(double *f, double *x, int N, int n, int m){
	double resp, h, tmp;
	int k, r, Np = pow(2, N);

	resp = 0.0;

	h = (x[Np] - x[0])/pow(2, n);
	
	if(n==0 && m==0) 
		resp = (f[0] + f[Np])/h/2.0;

	else if(m == 0){
		r = N - n;
		tmp = 0.0;
		for(k=1; k<=(int) pow(2, n-1); k++)
			tmp += f[(2*k-1)*((int) pow(2, r))];

		resp = 0.5*Romb(f, x, N, n-1, 0) + h*tmp;
	}

	else
		resp = 1.0/(pow(4, m) - 1.0)*(pow(4.0, m)*Romb(f, x, N, n, m-1) - Romb(f, x, N, n-1, m-1));

	return resp;
}


/*Evaluate the Fourier transform of the delta*/
void Fourier_delta(double *delta, double *deltak, double *r, double *k, int nr, int nk){
	int i, j, Nr = pow(2, nr), Nk = pow(2, nk);
	double *tmpr, a = 1.0e-3;

	tmpr = (double *)malloc((Nr+1)*sizeof(double));

	for(j=0;j<Nr+1;j++)
		tmpr[j] = r[j]*delta[j]*r[j];

	deltak[0] = Romb(tmpr, r, nr, nr, 1)*4.0*PI;

	for(i=1;i<Nk+1; i++){
		for(j=0;j<Nr+1;j++)
			tmpr[j] = r[j]*delta[j]*sin(k[i]*r[j]);

		deltak[i] = Romb(tmpr, r, nr, nr, 1)*4.0*PI/k[i];
	}
	free(tmpr);
}

/*Test the current range in k*/
double Fourier_test(double *delta, double *r, double k, int nr){
	int i, j, Nr = pow(2, nr);
	double *tmpr, resp, a = 1.0e-3;

	tmpr = (double *)malloc((Nr+1)*sizeof(double));

	if(k == 0.0){
		for(j=0;j<Nr+1;j++)
			tmpr[j] = r[j]*delta[j]*r[j];

		resp = Romb(tmpr, r, nr, nr, 1)*4.0*PI;
	}	
	
	else{
		for(j=0;j<Nr+1;j++)
			tmpr[j] = r[j]*delta[j]*sin(k*r[j]);

		resp = Romb(tmpr, r, nr, nr, 1)*4.0*PI/k;
	}
	
	free(tmpr);

	return resp;
}

/*Evaluate the MG contribution to the Poisson's equation*/
void Int_Poisson(double *Int, double *deltak, double *r, double *k, int nr, int nk, double n){
	int i, j, Nr = pow(2, nr), Nk = pow(2, nk);
	double *tmpk;

	tmpk = (double *)malloc((Nk+1)*sizeof(double));

	for(j=0;j<Nk+1;j++)
		tmpk[j] = k[j]*deltak[j]*k[j]*Eps(k[j], n);

	Int[0] = Romb(tmpk, k, nk, nk, 3)/(2.0*PI*PI);
	

	for(i=1;i<Nr+1; i++){
		for(j=0;j<Nk+1;j++)
			tmpk[j] = k[j]*deltak[j]*sin(k[j]*r[i])*Eps(k[j], n);

		Int[i] = Romb(tmpk, k, nk, nk, 3)/(2.0*PI*PI*r[i]);
	}
	free(tmpk);
}

/*Evaluate the inverse Fourier transform*/
void InvFourier_delta(double *Int, double *deltak, double *r, double *k, int nr, int nk){
	int i, j, Nr = pow(2, nr), Nk = pow(2, nk);
	double *tmpk;

	tmpk = (double *)malloc((Nk+1)*sizeof(double));

	for(j=0;j<Nk+1;j++)
		tmpk[j] = k[j]*deltak[j]*k[j];

	Int[0] = Romb(tmpk, k, nk, nk, 3)/(2.0*PI*PI);

	for(i=1;i<Nr+1; i++){
		for(j=0;j<Nk+1;j++)
			tmpk[j] = k[j]*deltak[j]*sin(k[j]*r[i]);

		Int[i] = Romb(tmpk, k, nk, nk, 3)/(2.0*PI*PI*r[i]);
	}
	free(tmpk);
}

/*Non-linear equation for \dot{v} = \ddot{\delta}*/
double dv(double n, double d, double v, double Int){
	double resp;
	resp = -(2.0 + dE(n)/E(n))*v + 4.0/3.0*(1.0/(1.0+d))*v*v + 3.0/2.0*cosmo.Om/pow(E(n),2)*exp(-3.0*n)*(1.0+d)*(d + Int);
	return resp;
}

/*Non-linear equation for \dot{v} = \ddot{\delta}*/
double dv_lin(double n, double d, double v, double Int){
	double resp;
	resp = -(2.0 + dE(n)/E(n))*v + 3.0/2.0*cosmo.Om/pow(E(n),2)*exp(-3.0*n)*(d + Int);
	return resp;
}

/*Equation for \dot{\delta}*/
double dd(double n, double d, double v){
	return v;
}

/*Evolve the non-linear EDO one steep in time using the Runge-Kutta 4*/
void Evolve(double n, double *r, double *k, double *delta, double *vel, double h, int nr, int nk){
	int i, Nr = pow(2, nr), Nk = pow(2, nk);
	double q1[Nr+1], q2[Nr+1], q3[Nr+1], q4[Nr+1], k1[Nr+1], k2[Nr+1], k3[Nr+1], k4[Nr+1], delta_tmp[Nr+1], deltak[Nk+1], Int[Nr+1], errF[Nr+1];

	/*Compute the error in the integrals*/
	Fourier_delta(delta, deltak, r, k, nr, nk);
	InvFourier_delta(Int, deltak, r, k, nr, nk);

	for(i=0;i<=Nr;i++){
		errF[i] = delta[i]/Int[i]; printf("%lf ", errF[i]); }
	printf("\n");
	

	/*Compute q1 and k1 of RK4*/
	Int_Poisson(Int, deltak, r, k, nr, nk, n);
	for(i=0;i<=Nr;i++)
		Int[i] = Int[i]*errF[i];

	for(i=0;i<Nr+1;i++){
		q1[i] = dd(n, delta[i], vel[i]);
		k1[i] = dv(n, delta[i], vel[i], Int[i]);
	}

	/*Compute q2 and k2 of RK4*/
	for(i=0;i<Nr+1;i++)	delta_tmp[i] = delta[i] + h/2.0*q1[i];

	if(delta[0] > 100.0){
		Fourier_delta(delta_tmp, deltak, r, k, nr, nk);
		InvFourier_delta(Int, deltak, r, k, nr, nk);

		for(i=0;i<=Nr;i++)
			errF[i] = delta_tmp[i]/Int[i];
	}

	Fourier_delta(delta_tmp, deltak, r, k, nr, nk);
	Int_Poisson(Int, deltak, r, k, nr, nk, n);
	for(i=0;i<=Nr;i++)
		Int[i] = Int[i]*errF[i];

	for(i=0;i<Nr+1;i++){
		q2[i] = dd(n+h/2.0, delta[i] + h/2.0*q1[i], vel[i] + h/2.0*k1[i]);
		k2[i] = dv(n+h/2.0, delta[i] + h/2.0*q1[i], vel[i] + h/2.0*k1[i], Int[i]);
	}

	/*Compute q3 and k3 of RK4*/
	for(i=0;i<Nr+1;i++)	delta_tmp[i] = delta[i] + h/2.0*q2[i];

	if(delta[0] > 100.0){
		Fourier_delta(delta_tmp, deltak, r, k, nr, nk);
		InvFourier_delta(Int, deltak, r, k, nr, nk);

		for(i=0;i<=Nr;i++)
			errF[i] = delta_tmp[i]/Int[i];
	}

	Fourier_delta(delta_tmp, deltak, r, k, nr, nk);
	Int_Poisson(Int, deltak, r, k, nr, nk, n);
	for(i=0;i<=Nr;i++)
		Int[i] = Int[i]*errF[i];

	for(i=0;i<Nr+1;i++){
		q3[i] = dd(n+h/2.0, delta[i] + h/2.0*q2[i], vel[i] + h/2.0*k2[i]);
		k3[i] = dv(n+h/2.0, delta[i] + h/2.0*q2[i], vel[i] + h/2.0*k2[i], Int[i]);
	}

	/*Compute q4 and k4 of RK4*/
	for(i=0;i<Nr+1;i++)	delta_tmp[i] = delta[i] + h*q3[i];

	if(delta[0] > 100.0){
		Fourier_delta(delta_tmp, deltak, r, k, nr, nk);
		InvFourier_delta(Int, deltak, r, k, nr, nk);

		for(i=0;i<=Nr;i++)
			errF[i] = delta_tmp[i]/Int[i];
	}

	Fourier_delta(delta_tmp, deltak, r, k, nr, nk);
	Int_Poisson(Int, deltak, r, k, nr, nk, n);
	for(i=0;i<=Nr;i++)
		Int[i] = Int[i]*errF[i];

	for(i=0;i<Nr+1;i++){
		q4[i] = dd(n+h, delta[i] + h*q3[i], vel[i] + h*k3[i]);
		k4[i] = dv(n+h, delta[i] + h*q3[i], vel[i] + h*k3[i], Int[i]);
	}

	/*Evolve delta and vel one step in time*/
	for(i=0;i<Nr+1;i++){
		delta[i] = delta[i] + h/6.0*(q1[i] + 2.0*q2[i] + 2.0*q3[i] + q4[i]);
		vel[i] = vel[i] + h/6.0*(k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]);
	}
}

/*Evolve the linear EDO one steep in time using the Runge-Kutta 4*/
void Evolve_Linear(double n, double *r, double *k, double *delta, double *vel, double h, int nr, int nk){
	int i, Nr = pow(2, nr), Nk = pow(2, nk);
	double q1[Nr+1], q2[Nr+1], q3[Nr+1], q4[Nr+1], k1[Nr+1], k2[Nr+1], k3[Nr+1], k4[Nr+1], delta_tmp[Nr+1], deltak[Nk+1], Int[Nr+1], errF[Nr+1];

	/*Compute the error in the integrals*/
	Fourier_delta(delta, deltak, r, k, nr, nk);
	InvFourier_delta(Int, deltak, r, k, nr, nk);

	for(i=0;i<=Nr;i++)
		errF[i] = delta[i]/Int[i];

	/*Compute q1 and k1 of RK4*/
	Int_Poisson(Int, deltak, r, k, nr, nk, n);
	for(i=0;i<=Nr;i++)
		Int[i] = Int[i]*errF[i];

	for(i=0;i<Nr+1;i++){
		q1[i] = dd(n, delta[i], vel[i]);
		k1[i] = dv_lin(n, delta[i], vel[i], Int[i]);
	}

	/*Compute q2 and k2 of RK4*/
	for(i=0;i<Nr+1;i++)	delta_tmp[i] = delta[i] + h/2.0*q1[i];

	Fourier_delta(delta_tmp, deltak, r, k, nr, nk);
	Int_Poisson(Int, deltak, r, k, nr, nk, n);
	for(i=0;i<=Nr;i++)
		Int[i] = Int[i]*errF[i];

	for(i=0;i<Nr+1;i++){
		q2[i] = dd(n+h/2.0, delta[i] + h/2.0*q1[i], vel[i] + h/2.0*k1[i]);
		k2[i] = dv_lin(n+h/2.0, delta[i] + h/2.0*q1[i], vel[i] + h/2.0*k1[i], Int[i]);
	}

	/*Compute q3 and k3 of RK4*/
	for(i=0;i<Nr+1;i++)	delta_tmp[i] = delta[i] + h/2.0*q2[i];

	Fourier_delta(delta_tmp, deltak, r, k, nr, nk);
	Int_Poisson(Int, deltak, r, k, nr, nk, n);
	for(i=0;i<=Nr;i++)
		Int[i] = Int[i]*errF[i];

	for(i=0;i<Nr+1;i++){
		q3[i] = dd(n+h/2.0, delta[i] + h/2.0*q2[i], vel[i] + h/2.0*k2[i]);
		k3[i] = dv_lin(n+h/2.0, delta[i] + h/2.0*q2[i], vel[i] + h/2.0*k2[i], Int[i]);
	}

	/*Compute q4 and k4 of RK4*/
	for(i=0;i<Nr+1;i++)	delta_tmp[i] = delta[i] + h*q3[i];

	Fourier_delta(delta_tmp, deltak, r, k, nr, nk);
	Int_Poisson(Int, deltak, r, k, nr, nk, n);
	for(i=0;i<=Nr;i++)
		Int[i] = Int[i]*errF[i];

	for(i=0;i<Nr+1;i++){
		q4[i] = dd(n+h, delta[i] + h*q3[i], vel[i] + h*k3[i]);
		k4[i] = dv_lin(n+h, delta[i] + h*q3[i], vel[i] + h*k3[i], Int[i]);
	}

	/*Evolve delta and vel one step in time*/
	for(i=0;i<Nr+1;i++){
		delta[i] = delta[i] + h/6.0*(q1[i] + 2.0*q2[i] + 2.0*q3[i] + q4[i]);
		vel[i] = vel[i] + h/6.0*(k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]);
	}
}


int main(int argc,char *argv[])
{
char profile[100], velfile[100], parmfile[100], linearfile[100];
FILE *prof, *outv, *outp, *outl;
double h, nini, rb, s, d0, n1, n2;
double n, d, v, nmax,  *nt, *nc, *dt, dc, dt_lin, errnt, errdt, errnc, errk, errk2, Kcut, Kini;
double *delta, *vel, *r, *k, *velk, d_tmp, adj, *delta0;
int i, j, Nr, Nk, *flag, r_min, nr, nk;


/*Confere the correct number of parameters*/
if (argc != 7){
	printf("Wrong number of arguments.\n");
	printf("arg1: f_{R0}.\n");
	printf("arg2: n.\n");
	printf("arg3: Profile file name.\n");
	printf("arg4: d_{0}.\n");
	printf("arg5: r_{b}.\n");	
	printf("arg6: Output prefix.\n\n");
	exit(0);
}

/*Give the name to the output files*/
sprintf(parmfile,"%s.dat", argv[6]);
sprintf(profile,"%s", argv[3]);

/*Set the cosmology*/
cosmo.H0 = 1.0/2997.92;	/*Hubble's constant/c in units of 1/Mpc*/
cosmo.w = -1.0;
cosmo.Om = 0.31;
cosmo.Ol = 1.0 - cosmo.Om;
cosmo.Ok = 1.0 - cosmo.Om - cosmo.Ol;

/*Set the MG parameters*/
hs.fr0 = atof(argv[1]);
hs.n = atof(argv[2]);

/*EDO numerical parameters*/
d0 = atof(argv[4]);	/*Initial profile amplitude*/
rb = atof(argv[5]);	/*The radius of teh perturbation*/
nini = log(1.0/501.0);	/*Initial time*/
nmax = 1.6;		/*The biggest time*/
h = 0.01;		/*Initial time step*/
r_min = 0;		/*Smallest radius evolved*/
nk = 6;
Nk = pow(2, nk);

/*Open the profile file*/
prof = fopen(profile, "r");
if (prof == NULL) {
	printf("Unable to open %s\n", profile);
	exit(0);
}

/*Read the initial profile*/
fscanf(prof, "%d", &Nr);
nr = (int) log(Nr)/log(2.0);

/*Alloc the main vectors*/
r = (double *)malloc((Nr+1)*sizeof(double));
delta = (double *)malloc((Nr+1)*sizeof(double));
delta0 = (double *)malloc((Nr+1)*sizeof(double));
vel = (double *)malloc((Nr+1)*sizeof(double));

nt = (double *)malloc((Nr+1)*sizeof(double));
nc = (double *)malloc((Nr+1)*sizeof(double));
dt = (double *)malloc((Nr+1)*sizeof(double));
flag = (int *)malloc((Nr+1)*sizeof(int));

k = (double *)malloc((Nk+1)*sizeof(double));

for(i=0;i<Nr+1;i++){
	fscanf(prof, "%lf %lf", &r[i], &delta0[i]);
	delta[i] = d0*delta0[i]/delta0[0];
	vel[i] = exp(nini)*delta[i];
	flag[i] = 0;
	dt[i] = 0.0;
	nt[i] = nmax;
	nc[i] = nmax;
}
fclose(prof);


/*Parameters for the numerical integrations*/
Kcut = 0.15;
for(i=0;i<Nk+1;i++)
	k[i] = Kcut*i/Nk;

n = nini;
/*Solve the EDO using Runge-Kutta 4*/
while(flag[0] == 0 && n<nmax){
	/*printf("eta = %.5f and delta[0] = %lf\n", n, delta[0]);*/

	/*Prepare the correct k array*/
	errk = fabs(Fourier_test(delta, r, Kcut,  nr)/Fourier_test(delta, r, 0.0,  nr));
	errk2 = fabs(Fourier_test(delta, r, Kcut*1.01,  nr)/Fourier_test(delta, r, 0.0,  nr));
	adj = 1.0;
	while((errk > 1.0e-5 || errk2 > errk) && (adj<1.1 || delta[0]<0.1)){
		Kcut *= 1.01;
		adj *= 1.01;
		errk = fabs(Fourier_test(delta, r, Kcut, nr)/Fourier_test(delta, r, 0.0, nr));
		errk2 = fabs(Fourier_test(delta, r, Kcut*1.01, nr)/Fourier_test(delta, r, 0.0, nr));
	}

	/*printf("Kcut = %lf, adj = %lf\n", Kcut, adj);*/

	for(i=0;i<Nk+1;i++)
		k[i] = 1.0*Kcut*i/Nk;

	/*Evolve all shells one step in time*/
	Evolve(n, r, k, delta, vel, h, nr, nk);
	n = n + h;

	for(j=1;j<Nr;j++)
		if(delta[j] > delta[j-1])
			delta[j] = (delta[j-1] + delta[j+1])/2.0;

	if(delta[0] <= d_tmp) break;

	/*See if each shell is in turn-around or collapse and take the main information*/
	for(j=r_min;j<Nr+1;j++){
		if(vel[j] - 3.0*(1.0 + delta[j]) > 0.0 && flag[j] == 0){
			nt[j] = n - h/2.0;
			errnt = h/2.0;
			dt[j] = (delta[j] + d_tmp)/2.0;
			errdt = (delta[j] - d_tmp)/2.0;
			flag[j] = 1;
		}

		if(delta[j] > 200.0){
			nc[j] = n;
			errnc = h/2.0;
			r_min = j + 1;
		}
	}
	d_tmp = delta[0];

	/*Decrease the time step in non-linear regime*/
	if(delta[0] > 1.0 && h == 0.01){
		n1 = n;
		h = 0.001;
	}

	if(delta[0] > 100.0 && h == 0.001){
		n2 = n;
		h = 0.0001;
	}
}

/*free the memory*/
free(r);
free(delta);
free(vel);
free(k);

/*Print some usefull parameters of this run*/
outp = fopen(parmfile, "w");
if (outp == NULL) {
	printf("Unable to open %s\n", parmfile);
	exit(0);
}

fprintf(outp, "%lf %lf %lf %lf", exp(nt[0]), exp(nt[0])*errnt, dt[0], errdt);

fclose(outp);

return 0;
}
