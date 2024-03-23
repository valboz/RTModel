// bumper.cpp 
// Implementation of all methods in bumper.h

#define _CRT_SECURE_NO_WARNINGS
#define _USE_MATH_DEFINES
#include "bumper.h"
#include <math.h>
#include <malloc.h>

const double epsilon = 1.e-99;
	
bumper::bumper(double *pr,int ns){
	nps=ns;
	p0=(double *)malloc(sizeof(double)*nps);
	dp=(double *)malloc(sizeof(double)*nps);
	for(int i=0;i<nps;i++){
		p0[i]=pr[i];
	}
	curv=cov=0;

	next=0;
}

bumper::~bumper(){
	free(p0);
	free(dp);
	if(curv) free(curv);
	if (cov) free(cov);
}

void bumper::SetCurvature(double *sr,double nrm){
	curv=(double *)malloc(sizeof(double)*nps*nps);
	for(int i=0;i<nps*nps;i++){
		curv[i]=sr[i]*nrm;
	}
}

void bumper::SetCovariance(double* sr, double nrm) {
	cov = (double*)malloc(sizeof(double) * nps * nps);
	for (int i = 0; i < nps * nps; i++) {
		cov[i] = sr[i] / nrm;
	}
	//curv = (double*)malloc(sizeof(double) * nps * nps);
	//Inverse(sr, curv, nps);
}

void bumper::UpdateCurvature(double bumperpower){
	for(int i=0;i<nps*nps;i++){
		curv[i]/=(bumperpower*bumperpower);
	}
}

double bumper::distance(double *pr){
	double dist;
	for(int i=0;i<nps;i++){
		dp[i]=pr[i]-p0[i]+epsilon;
	}
	dist=0;
	for(int i=0;i<nps;i++){
		dist+=dp[i]*curv[i*nps+i]*dp[i];
		for(int j=0;j<i;j++){
			dist+=2*dp[i]*curv[i*nps+j]*dp[j];
		}
	}
	return dist;
}

void bumper::signCovariance(int j) {
	for (int i = 0; i < nps; i++) {
		cov[j * nps + i] = -cov[j * nps + i];
		cov[i * nps + j] = -cov[i * nps + j];
	}
}

void bumper::flipCovariance(int i,int j) {
	double sc;
	for (int k = 0; k < nps; k++) {
		sc = cov[k * nps + i];
		cov[k * nps + i] = cov[k * nps + j];
		cov[k * nps + j] = sc;
	}
	for (int k = 0; k < nps; k++) {
		sc = cov[i * nps + k];
		cov[i * nps + k] = cov[j * nps + k];
		cov[j * nps + k] = sc;
	}
}

double Determinant(double *A,int n){
	double *U,L,M;
	double det=1.;
	int iM;

	U=(double *)malloc(sizeof(double)*n*n);

	for(int i=0;i<n;i++){
		for(int j=0;j<n;j++){
			U[i*n+j]=A[i*n+j];
		}
	}

	for(int j=0;j<n-1;j++){
		for(int i=j+1;i<n;i++){
// find row with largest pivot
			M=0.;
			iM=j;
			for(int i1=j;i1<n;i1++){
				if (M < fabs(U[i1 * n + j])) {
					M = fabs(U[i1 * n + j]);
					iM = i1;
				}
			}
			if(M==0.){ 
				free(U); 
				return 0.;
			}
// select row with largest pivot as the next to be processed
			for(int j1=j;j1<n;j1++){
				L=U[j*n+j1];
				U[j*n+j1]=U[iM*n+j1];
				U[iM*n+j1]=-L;
			}

			L = (fabs(U[j * n + j])>1.e-100) ? U[i * n + j] / U[j * n + j] : 1.0;
			for(int k=j;k<n;k++){
				U[i*n+k]-=L*U[j*n+k];
			}
		}
	}

	for(int i=0;i<n;i++){
		det*=U[i*n+i];
	}

	free(U);
	return det;
}


void Inverse(double* source, double* dest, int n) {
	double p1;
	int j1, k1;
	double* minor;
	minor = (double*)malloc(sizeof(double) * (n - 1) * (n - 1));
	p1 = Determinant(source, n)+1.e-100;
	for (int i = 0; i < n; i++) {
		for (int i2 = 0; i2 < n; i2++) {
			for (int j = 0; j < n - 1; j++) {
				for (int k = 0; k < n - 1; k++) {
					j1 = (j < i2) ? j : j + 1;
					k1 = (k < i) ? k : k + 1;
					minor[j * (n - 1) + k] = source[j1 * n + k1];
				}
			}
			dest[i * n + i2] = Determinant(minor, n - 1) / p1 * (((i + i2) % 2) ? 1 : -1);
		}
	}
	free(minor);
}

void CombineCovariances(bumper* scanbumper, bumper* scanbumper2,double *Cov, double *Curv, int nps) {
	for (int i = 0; i < nps; i++) {
		for (int j = 0; j < nps; j++) {
			Cov[i * nps + j] = scanbumper->cov[i * nps + j] + scanbumper2->cov[i * nps + j];
		}
	}
	Inverse(Cov, Curv, nps);
}