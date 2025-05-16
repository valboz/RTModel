// LevMarFit.h
// Definition of the LevMar and bumper classes for Levenberg-Marquardt fitting

#include <cstdio>
#include "bumper.h"
#include <regex>
#include <filesystem>

using namespace std;
using namespace std::filesystem;

#ifndef _LevMarFit
#define _LevMarFit
#define __unmanaged

#include <VBMicrolensingLibrary.h>

class LevMar {

	VBMicrolensing* VBM;

	char eventname[512], filename[30], filnum[10], outdir[30], satellitedir[256];
	char modelcode[16];
	int error;
	int flagblending;
	path exedir;
	bool astrometric;
	int nlinpar;

	double tim0, tm;

	int nps;
	double* sigmapr, * leftlim, * rightlim;

	int* filter, * satel, nfil, np, OGLE;
	double* t, * y, * w, *y1a, *y2a, * delta, * maxdelta, * Curv, * A, * B, * B0, * Cov, * fb, ** Gr, * dFdp, * errs;
	double* cN, * cE, * wcN, * wcE, *c1s, *c2s, *c1l, *c2l;
	double* pr, * prn, * sumy, * sumy2, * sumsigma, * sumfy, * sumf, * sumf2, * limbdarks;
	double* sumsigmaN, * sumsigmaE, * sumcN, * sumcE, * sumc1, * sumc2;
	int* sizes, * starts;

	int consnumber, *consindex;
	double* constraints, * consleft, * consright, *consvars;
	int modnumber;

	double Tol;

//	void (LevMar::* PrintOut)(double*);
//	void (LevMar::* PrintFile)(FILE*, double, bool);
	void PrintOut(double*);
	void PrintFile(char *filename, double, bool);

	bumper* stepchain, * bumperlist, * laststep;


public:
	LevMar(int, char**);
	~LevMar();

	void ReadFiles(int, char**);
	int InitCond(double* presigmapr, double* preleftlim, double* prerightlim);
	void ReadCurve();
	void ReadAncillary();
	void ReadOptions();
	int Run();
	double ChiSquared(double*);
	void Grad();
	void EvaluateModel(double *pr, int filter, int ips);
	void Covariance();
	double ComputeConstraint(double *pr, int i);

	void PrintOutPS(double*);
	void PrintOutPX(double*);
	void PrintOutBS(double*);
	void PrintOutBO(double*);
	void PrintOutLS(double*);
	void PrintOutLX(double*);
	void PrintOutLO(double*);
	void PrintOutLK(double*);
	void PrintOutTS(double*);
	void PrintOutTX(double*);

	void PrintFilePS(FILE*, double, bool);
	void PrintFilePX(FILE*, double, bool);
	void PrintFileBS(FILE*, double, bool);
	void PrintFileBO(FILE*, double, bool);
	void PrintFileLS(FILE*, double, bool);
	void PrintFileLX(FILE*, double, bool);
	void PrintFileLO(FILE*, double, bool);
	void PrintFileLK(FILE*, double, bool);
	void PrintFileTS(FILE*, double, bool);
	void PrintFileTX(FILE*, double, bool);

};

double Determinant(double*, int);
void Inverse(double*, double*, int);

#endif