// bumper.h
// Definition of the bumper class containing information on a single model

#ifndef _bumper
#define _bumper

// Bumper class is the penalty function to be used to fill minima in chi square
class bumper{
public:
	double *p0;
	double *dp;
	double *curv, *cov;
	double Amp;
	int nps;
	char modelcode[16];
	int il;
	bumper(double *,int);
	~bumper();
	void SetCurvature(double *,double);
	void SetCovariance(double*, double);
	void UpdateCurvature(double);
	void signCovariance(int);
	void flipCovariance(int, int);
	double distance(double *);
	bumper *next;
};

double Determinant(double *,int);
void Inverse(double*, double*, int);
void CombineCovariances(bumper*, bumper*, double *Cov, double * Curv, int);

#endif