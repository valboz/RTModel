// LevMarFit.h
// Definition of the LevMar and bumper classes for Levenberg-Marquardt fitting

#include <stdio.h>
#include "bumper.h"
#include <regex>
#include <filesystem>

using namespace std;
using namespace std::filesystem;

#ifndef _LevMarFit
#define _LevMarFit
#define __unmanaged

#include <VBBinaryLensingLibrary.h>

class LevMar{

		VBBinaryLensing *VBBL;

		char eventname[512],filename[30],filnum[10],outdir[30],satellitedir[256];
		char modelcode[16];
		int error;
		int flagblending;
		path exedir;

		double tim0,tm;

		double (VBBinaryLensing::*model)(double *,double);
		int nps;
		double *sigmapr,*leftlim,*rightlim;

		int *filter,*satel,nfil,np,OGLE;
		double *t,*y,*w,*delta,*maxdelta,*Curv,*A,*B,*B0,*Cov,*fb,**Gr,*dFdp,*errs;
		double *pr,*prn,*sumy, *sumy2, *sumsigma,*sumfy,*sumf,*sumf2,*limbdarks;

		double Tol;

		void (LevMar::*PrintOut)(double *);
		void (LevMar::*PrintFile)(FILE *,double, bool);
		
		bumper *stepchain,*bumperlist,*laststep;


	public:
		LevMar(int,char **);
		~LevMar();

		void ReadFiles(int,char **);
		int InitCond(double *presigmapr, double *preleftlim,double *prerightlim);
		void ReadCurve();
		void Run();
		double ChiSquared(double *);
		void Grad();
		void Covariance();

		void PrintOutPS(double *);
		void PrintOutPX(double *);
		void PrintOutBS(double *);
		void PrintOutBO(double *);
		void PrintOutLS(double *);
		void PrintOutLX(double *);
		void PrintOutLO(double *);
		void PrintOutLK(double*);

		void PrintFilePS(FILE *,double, bool);
		void PrintFilePX(FILE *,double, bool);
		void PrintFileBS(FILE *,double, bool);
		void PrintFileBO(FILE *,double, bool);
		void PrintFileLS(FILE *,double, bool);
		void PrintFileLX(FILE *,double, bool);
		void PrintFileLO(FILE *,double, bool);
		void PrintFileLK(FILE*, double, bool);

};

double Determinant(double *,int);
void Inverse(double*, double*, int);

#endif