// bumper.h
// Definition of the bumper class containing information on a single model

#ifndef _bumper
#define _bumper

#include<stdio.h>

// Bumper class is the penalty function to be used to fill minima in chi square
class bumper{
public:
	double *p0;
	double *dp;
	double *curv, *cov;
	double Amp;
	double tanomaly;
	double resanomaly;
	int nps;
	char modelcode[16];
	char *buffer;
	int il;
	bool duplicate;
	bumper(double *,int);
	~bumper();
	void SetCurvature(double *,double);
	void SetCovariance(double*, double);
	void UpdateCurvature(double);
	void signCovariance(int);
	void flipCovariance(int, int);
	double distance(double *);
	void SetBuffer(FILE*, int, int);
	bumper *next;
};

double Determinant(double *,int);
void Inverse(double*, double*, int);
void CombineCovariances(bumper*, bumper*, double *Cov, double * Curv, int);

#endif