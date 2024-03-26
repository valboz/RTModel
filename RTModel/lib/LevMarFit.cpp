// LevMarFit.cpp 
// Implementation of all methods in LevMarFit.h for Levenberg-Marquardt fitting

#define _CRT_SECURE_NO_WARNINGS
#define _USE_MATH_DEFINES
#include "LevMarFit.h"
#include "bumper.h"
#include <VBBinaryLensingLibrary.h>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <malloc.h>
#include <regex>
#include <filesystem>

using namespace std;
using namespace std::filesystem;



int nlc = 5; // Number of models to be calculated from the same initial condition using the bumper method
int maxsteps = 50; // Maximum number of steps in each fit
double maxtime=600.0; // Maximum time in seconds for total execution
double bumperpower = 2.0; // Repulsion factor of bumpers
double maxbumpcount = 25;


const double epsilon = 1.e-100;

LevMar::LevMar(int argc,char *argv[]){
 	printf("******************************************\n");
	printf("*************      LevMar     **********\n");
	printf("******************************************\n\n\n");
	printf("This program fits an event from specific initial conditions\n\n");
	error=0;
// Setting default values
	Tol=1.e-2;
	VBBL =new VBBinaryLensing;
	VBBL->Tol = Tol;
	VBBL->RelTol = 0.001;
	VBBL->parallaxsystem = 1;

	ReadFiles(argc,argv);

}

LevMar::~LevMar(){
	delete VBBL;
	if(error<9){
		free(t);
		free(y);
		free(w);
		free(satel);
		free(filter);
		free(delta);
		free(maxdelta);
		free(Curv);
		free(A);
		free(B);
		free(B0);
		free(Cov);
		free(fb);
		free(pr);
		free(prn);
		free(errs);
		free(sumy);
		free(sumy2);
		free(sumsigma);
		free(sumfy);
		free(sumf);
		free(sumf2);
		free(limbdarks);
		free(dFdp);
		for(int i=0;i<nps;i++){
			free(Gr[i]);
		}
		free(Gr);


		bumper *scanbumper=bumperlist,*scanbumper2;
		while(scanbumper){
			scanbumper2=scanbumper->next;
			delete scanbumper;
			scanbumper=scanbumper2;
		}


		scanbumper=stepchain;
		while(scanbumper){
			scanbumper2=scanbumper->next;
			delete scanbumper;
			scanbumper=scanbumper2;
		}

	}
	if(error<10){
		free(sigmapr);
		free(leftlim);
		free(rightlim);
	}
}

void LevMar::ReadFiles(int argc,char *argv[]){
	FILE *f;
	char buffer[3200], initcondfile[256];
	char command[256];
	double value;

	try{

	// Reading directory and event information

	exedir = current_path();
	strcpy(satellitedir, exedir.string().c_str());

	if (argc > 2) {
		strcpy(eventname, argv[1]);
		strcpy(outdir, argv[2]);
		if (argc > 3) {
			strcpy(satellitedir, argv[3]);
		}
	}
	else {
		printf("\n\nEvent name? ");
		scanf("%s", eventname);
		//sprintf(eventname, "WDC10193");
	}

	printf("\n\n- Event: %s\n", eventname);

	current_path(eventname);

	//if(argc>2){
	//	strcpy(eventname,argv[1]);
	//	strcpy(outdir,argv[2]);
	//}else{
	//	*eventname = 0;
	//	f = fopen("Model.ini", "r");
	//	fscanf(f, "%s", outdir);
	//	fclose(f);
	//	nlc = 1;
	//}
	

/* Reading coordinates of the model */

	if (exists("ini")) {
		current_path("ini");
		f = fopen("LevMar.ini", "r");
		if (f != 0) {
			printf("\n\n- Reading options in LevMar.ini");
			while (!feof(f)) {
				int red = fscanf(f, "%s %s %lf", command, buffer, &value);
				if (red < 1) {
					command[0] = 0;
					//if (red != 0) {
					//	printf("\n\n!!! Bad command in Reader.ini");
					//	return -1;
					//};
				}
				if (strcmp(command, "nfits") == 0) {
					nlc = (int) value;
				}
				if (strcmp(command, "maxsteps") == 0) {
					maxsteps = (int) value;
				}
				if (strcmp(command, "timelimit") == 0) {
					maxtime = value;
				}
				if (strcmp(command, "bumperpower") == 0) {
					bumperpower = value;
				}
			}
			fclose(f);
		}
		else {
			printf("\n\n- Default options:");
		}
	}

	current_path(eventname);
	current_path("Data");

	auto searchstring = regex(".*\\.coordinates");
	for (auto const& itr : directory_iterator(".")) {
		string curfile = (itr).path().filename().string();
		if (regex_match(curfile, searchstring)) {
			VBBL->SetObjectCoordinates((char*) curfile.c_str(), satellitedir);
			printf("\n- Coordinates set.");
			break;
		}
	}

	current_path(eventname);

	// Establishing model to be fit
	strcpy(modelcode,outdir);
	switch(modelcode[0]){
		case 'P':
			if(modelcode[1]=='X'){
				model=&(VBBinaryLensing::ESPLLightCurveParallax);
				nps=6;
				double presigmapr[] ={.5,.5,5.,4.6,1,1};
				double preleftlim[]={-13.,-6.9,-10.e100,-11.5,-10.,-10.};
				double prerightlim[]={.7,6.9,10.e100,0.0,10.,10.};		
				strcpy(initcondfile, "InitCondPX.txt");
				InitCond(presigmapr,preleftlim,prerightlim);
				pr[1]=log(pr[1]);
				pr[3]=log(pr[3]);
				PrintOut=&LevMar::PrintOutPX;
				PrintFile=&LevMar::PrintFilePX;
			}else{
				model=&(VBBinaryLensing::ESPLLightCurve);
				nps=4;
				double presigmapr[] ={.5,.5,5.,4.6};
				double preleftlim[]={-13.,-6.9,-10.e100,-11.5};
				double prerightlim[]={.7,6.9,10.e100,0.0};					
				InitCond(presigmapr, preleftlim, prerightlim);
				pr[0]=log(pr[0]);
				pr[1]=log(pr[1]);
				pr[3]=log(pr[3]);
				PrintOut=&LevMar::PrintOutPS;
				PrintFile=&LevMar::PrintFilePS;
			}
			current_path(exedir);
			current_path("..");
			current_path("data");
			VBBL->LoadESPLTable("ESPL.tbl");
			current_path(eventname);
			break;
		case 'B':
			if(modelcode[1]=='O'){
				model=&(VBBinaryLensing::BinSourceSingleLensXallarap);
				nps=10;
				double presigmapr[] ={1,1,1,15,3,3,1,3,6,3};
				double preleftlim[]={-3.,-1.e100,-6.9,-11.5,-3.,-3.,0,-3,-6,-4.6};
				double prerightlim[]={3.,1.e100,6.9,0.,3.,3.,1,3,6,1,+4.6};					
				InitCond(presigmapr, preleftlim, prerightlim);
				pr[2]=log(pr[2]);
				pr[3]=log(pr[3]);
				pr[9] = log(pr[9]);
				PrintOut=&LevMar::PrintOutBO;
				PrintFile=&LevMar::PrintFileBO;
			}else{
				model=&(VBBinaryLensing::BinSourceExtLightCurve);
				nps=7;
				double presigmapr[] ={.1,.4,1,1,1,1,4.6};
				double preleftlim[]={-6.9,-11.5,0.,0.,-10.e100,-10.e100,-11.5};
				double prerightlim[]={6.9,11.5,3.,3.,10.e100,10.e100,0.0};					
				InitCond(presigmapr, preleftlim, prerightlim);
				pr[0]=log(pr[0]);
				pr[1]=log(pr[1]);
				pr[6] = log(pr[6]);
				PrintOut=&LevMar::PrintOutBS;
				PrintFile=&LevMar::PrintFileBS;
			}			
			current_path(exedir);
			current_path("..");
			current_path("data");
			VBBL->LoadESPLTable("ESPL.tbl");
			current_path(eventname);
			break;
		case 'L':
			if(modelcode[1]=='X'){
				model=&(VBBinaryLensing::BinaryLightCurveParallax);
				nps=9;
				double presigmapr[] ={.1,.4,.1,.1,4.6,.1,1.,1.,1.};
				double preleftlim[]={-4.0,-11.5,-3.,-12.56,-11.5,-6.9,-10.e100,-3.,-3.};
				double prerightlim[]={3.0,11.5,3.,12.56,-2.5,7.6,10.e100,3.,3.};					
				InitCond(presigmapr, preleftlim, prerightlim);
				pr[0]=log(pr[0]);
				pr[1]=log(pr[1]);
				pr[4]=log(pr[4]);
				pr[5]=log(pr[5]);
				PrintOut=&LevMar::PrintOutLX;
				PrintFile=&LevMar::PrintFileLX;
			}else{ 
				if(modelcode[1]=='O'){
					model=&(VBBinaryLensing::BinaryLightCurveOrbital);
					nps=12;
					double presigmapr[]={1.,2.,1.,5.,15.6,2.,10.,3.,3.,1.,1.,3.};
					double preleftlim[]={-4.0,-11.5,-3.,-12.56,-11.5,-6.9,-10.e100,-3.,-3.,-1,-1,1.e-7};
					double prerightlim[]={3.0,11.5,3.,12.56,-2.5,7.6,10.e100,3.,3.,1,1,1};
					InitCond(presigmapr, preleftlim, prerightlim);
					pr[0]=log(pr[0]);
					pr[1]=log(pr[1]);
					pr[4]=log(pr[4]);
					pr[5]=log(pr[5]);
					PrintOut=&LevMar::PrintOutLO;
					PrintFile=&LevMar::PrintFileLO;
				}else{
					if (modelcode[1] == 'P') {
						model = &(VBBinaryLensing::BinaryLightCurveOrbital);
						nps = 12;
						double presigmapr[] = { 1.,2.,1.,5.,15.6,2.,10.,.0001,.0001,1.,1.,3. };
						double preleftlim[] = { -4.0,-11.5,-3.,-12.56,-11.5,-6.9,-10.e100,-.0001,-.0001,-1,-1,1.e-7 };
						double prerightlim[] = { 3.0,11.5,3.,12.56,-2.5,7.6,10.e100,.0001,.0001,1,1,1 };
						InitCond(presigmapr, preleftlim, prerightlim);
						pr[0] = log(pr[0]);
						pr[1] = log(pr[1]);
						pr[4] = log(pr[4]);
						pr[5] = log(pr[5]);
						pr[7] = 0.;
						pr[8] = 0.;
						PrintOut = &LevMar::PrintOutLO;
						PrintFile = &LevMar::PrintFileLO;
					}
					else {
						if (modelcode[1] == 'K') {
							model = &(VBBinaryLensing::BinaryLightCurveKepler);
							nps = 14;
							double presigmapr[] = { 1.,2.,1.,5.,15.6,2.,10.,3.,3.,1.,1.,3., 3., 3. };
							double preleftlim[] = { -4.0,-11.5,-3.,-12.56,-11.5,-6.9,-10.e100,-3.,-3.,-1,-1,1.e-7, -10,0.5001 };
							double prerightlim[] = { 3.0,11.5,3.,12.56,-2.5,7.6,10.e100,3.,3.,1,1,1,10,10 };
							InitCond(presigmapr, preleftlim, prerightlim);
							pr[0] = log(pr[0]);
							pr[1] = log(pr[1]);
							pr[4] = log(pr[4]);
							pr[5] = log(pr[5]);
							PrintOut = &LevMar::PrintOutLK;
							PrintFile = &LevMar::PrintFileLK;
						}
						else {
							model = &(VBBinaryLensing::BinaryLightCurve);
							nps = 7;
							double presigmapr[] = { .1,.4,.1,.1,4.6,.1,1. };
							double preleftlim[] = { -4.0,-11.5,-3.,-12.56,-11.5,-6.9,-10.e100 };
							double prerightlim[] = { 3.0,11.5,3.,12.56,-2.5,7.6,10.e100};
							InitCond(presigmapr, preleftlim, prerightlim);
							pr[0] = log(pr[0]);
							pr[1] = log(pr[1]);
							pr[4] = log(pr[4]);
							pr[5] = log(pr[5]);
							PrintOut = &LevMar::PrintOutLS;
							PrintFile = &LevMar::PrintFileLS;
						}
					}
				}
			}
			break;
	}
	printf("\n- Model: %s",modelcode);

	}catch(...){
		error=10;
	}

	if(!error){
		// Reading Light Curve
		ReadCurve();
	}

}

int LevMar::InitCond(double* presigmapr, double* preleftlim, double* prerightlim) {
	char buffer[3200],initcondfile[256];
	sigmapr = (double*)malloc(sizeof(double) * nps); 
	leftlim = (double*)malloc(sizeof(double) * nps); 
	rightlim = (double*)malloc(sizeof(double) * nps); 
	pr = (double*)malloc(sizeof(double) * nps); 
	sprintf(initcondfile, "InitCond%c%c.txt", modelcode[0],modelcode[1]);
	current_path("InitCond");
	FILE *f = fopen(initcondfile, "r"); 
	int npeaks, ninit, incond; 
	fscanf(f, "%d %d", &npeaks, &ninit); 
	incond = atoi(outdir + 2); 
	if (incond >= ninit) nps=nps/0; 
	fgets(buffer, 3200, f);
	for (int i = 0; i < npeaks + incond; i++) {
			fgets(buffer,3200,f); 
	}
	for (int i = 0; i < nps; i++) {
		sigmapr[i] = presigmapr[i]; 
		leftlim[i] = preleftlim[i]; 
		rightlim[i] = prerightlim[i]; 
		fscanf(f, "%lg", &(pr[i])); 
		if (rightlim[i] > 1.e90) {
			rightlim[i] = pr[i] + 300.0;
			leftlim[i] = pr[i] - 300.0;
		}
	}
	fclose(f); 
	current_path(eventname);
	return 0;
}

void LevMar::ReadCurve(){
	FILE *f;
	int k;
	double *normfacs;

	try{
	printf("\n\n- Reading data \n");

	f=fopen("LCToFit.txt","r");
    fscanf(f,"%d",&np);
	filter=(int *) malloc(sizeof(int)*np);
	t=(double *) malloc(sizeof(double)*np);
	y=(double *) malloc(sizeof(double)*np);
	w=(double *) malloc(sizeof(double)*np);
	satel=(int *) malloc(sizeof(int)*np);
	delta=(double *) malloc(sizeof(double)*nps);
	maxdelta=(double *) malloc(sizeof(double)*nps);
	B=(double *) malloc(sizeof(double)*nps);
	B0=(double *) malloc(sizeof(double)*(nps));
	A=(double *) malloc(sizeof(double)*nps*nps);
	Curv=(double *) malloc(sizeof(double)*nps*nps);
	Cov=(double *) malloc(sizeof(double )*nps*nps);
	fb=(double *) malloc(sizeof(double)*np*(nps+1));


	nfil=1;
	for(int i=0;i<np;i++){
		fscanf(f,"%d %lf %lf %lf %d",&(filter[i]),&(t[i]),&(y[i]),&(w[i]),&(satel[i]));
		if((i!=0)&&(filter[i]>filter[i-1])){
			nfil++;
		}
		w[i]=1/(w[i]);
	}
	fclose(f);

	sumf=(double *) malloc(sizeof(double)*nfil);
	sumf2=(double *) malloc(sizeof(double)*nfil);
	sumfy=(double *) malloc(sizeof(double)*nfil);
	sumy=(double *) malloc(sizeof(double)*nfil);
	sumy2=(double *) malloc(sizeof(double)*nfil);
	sumsigma=(double *) malloc(sizeof(double)*nfil);
	pr=(double *) realloc(pr,sizeof(double)*(nps+2*nfil));
	prn=(double *) malloc(sizeof(double)*(nps+2*nfil));
	errs=(double *) malloc(sizeof(double)*(nps+2*nfil));
	dFdp=(double *) malloc(sizeof(double )*2*nfil*nps);

	Gr=(double **) malloc(sizeof(double *)*nps);
	for(int i=0;i<nps;i++){
		Gr[i]=(double *) malloc(sizeof(double)*np);
	}

	current_path("Data");

	// If Normalization.txt is present, use numbers therein to normalize datasets.

	normfacs = (double *)malloc(sizeof(double)*nfil);
	for (int i = 0; i < nfil; i++) {
		normfacs[i] = 1.;
	}
	if (exists("Normalizations.txt")) {
		f = fopen("Normalizations.txt", "r");
		for (int i = 0; i < nfil; i++) {
			fscanf(f, "%lf", &normfacs[i]);
			printf("\nNormalization %d: %lf", i, normfacs[i]);
		}
		fclose(f);
	}

	k = 0;
	sumsigma[k] = sumy[k] = sumy2[k] = 0;
	for (int i = 0; i<np; i++) {
		w[i] /= normfacs[filter[i]];
		if ((i != 0) && (filter[i]>filter[i - 1])) {
			k++;
			sumsigma[k] = sumy[k] = sumy2[k] = 0;
		}
		sumsigma[k] += w[i] * w[i];
		sumy[k] += w[i] * w[i] * y[i];
		sumy2[k] += w[i] * w[i] * y[i] * y[i];
	}

	for(int i=0;i<nps;i++){
		maxdelta[i]=sigmapr[i];
	}

	// If LimbDarkening.txt is present, use number therein for linear limb darkening coefficient
	limbdarks = (double *)malloc(sizeof(double)*nfil);
	for (int i = 0; i < nfil; i++) {
		limbdarks[i] = 0;
	}
	if(exists("LimbDarkening.txt")){
		f = fopen("LimbDarkening.txt", "r");
		for (int i = 0; i < nfil; i++) {
			fscanf(f, "%lf", &limbdarks[i]);
			printf("\nLimbDarkening %d: %lf", i, limbdarks[i]);
		}
		fclose(f);
	}
	current_path("..");

	//jfile=_findfirst("ZeroBlending.txt",&str2file); // Option for fitting with zero blending, to be elaborated
	//if(jfile!=-1){

	//	f=fopen("FilterToData.txt","r");
	//	k=0;
	//	while(k>=0){
	//		k++;
	//		fscanf(f,"%s",filnum);
	//		if(filnum[0]=='O'){
	//			OGLE=k;
	//			k=-1;
	//		}
	//		if(feof(f)) k=-1;
	//	}
	//	fclose(f);
	//	OGLE--;
	//}else{
	//	OGLE=-1;
	//}

	}catch(...){
		error=9;
	}
	free(normfacs);

}

void LevMar::Run(){
	FILE *f;
	int il,k,ichi,flag,ilam,bumpnum,bumpcounter;
	double minchi,bestchi,c1,c0,oldlambda,lambda,inclambda,fac,fac2;
	bumper *scanbumper,*scanbumper2;

	if(!error){
		try{
		// Initializing step chain and bumper list

		printf("\n\n- Initializing step chain \n");

		stepchain=new bumper(pr,nps);
		laststep=stepchain;
		bumperlist=0;

		// Going to right directory

		if(!exists("PreModels")) create_directory("PreModels");
		current_path("PreModels");
		if (!exists(outdir)) create_directory(outdir);
		current_path(outdir);
		il = 0;
		if (f = fopen("nlc.dat", "r")) {
			int i1, i2;
			fscanf(f, "%d %d", &i1, &i2);
			if (i2 == i1 + 1) {
				il = nlc;
			}
			fclose(f);
		}
		if(il==0){
			// Removing all existing files
			auto searchstring = regex(".*");
			for (auto const& itr : directory_iterator(".")) {
				string curfile = (itr).path().filename().string();
				if (is_regular_file(curfile)) {
					remove(curfile);
				}
			}
			// nlc.dat contains the number of models calculated so far and the total number to be computed
			f = fopen("nlc.dat", "w");
			fprintf(f, "%d %d", -1, nlc);
			fclose(f);
			// minchi.dat contains the best chi square found so far in all kinds of models for this event
			minchi = bestchi = 1.e100;
		}


		time_t ltime,rtime;
		time(&ltime);

		// Calculate nlc models starting from the initial condition
		// il is the number of the current model being calculated
		while(il<nlc){
			printf("\n*************************************");
			printf("\n************ Curve %d *************",il);
			printf("\n*************************************\n");

			for(int i=0;i<nps;i++){
				pr[i]=(laststep->p0)[i];
			}
			bumpcounter = 0;

			(this->*PrintOut)(pr);
			/* Chi Squared at initial conditions */

			c1=ChiSquared(pr);
			printf("\nStarting chi2 = %lf",c1);

			// Saving to file
			sprintf(filename, "%s-stepchain%d.dat", modelcode, il);
			f = fopen(filename, "a");
			(this->*PrintFile)(f, c1, false);
			fprintf(f, "\n");
			fclose(f);

			// Initializing the main cycle

			lambda=3.;
			inclambda=3.;
			bumpnum=1;
			k=1;
			ichi=0;
			// Main fit
			while((ichi<3)&&(lambda<1.e10)&&(k<=maxsteps)){
				c0=c1;
				printf("\nStep %d\n",k++);

				/* Calculation of the gradient */
			
				Grad();

				// Preparing the matrices A=(Curv + lambda diag(Curv))
				// with Curv = Gr^T Gr.
				// B=Gr (y-f)*w

				oldlambda=lambda;
				lambda/=inclambda;

				// Levenberg-Marquardt with parameter lambda
				ilam=0;
				while((c1>=c0)&&(ilam<20)){
					for(int i=0;i<nps;i++){
						for(int j=0;j<nps;j++){
							A[i*nps+j]=Curv[i*nps+j];
							if(i==j){
								A[i*nps+j]+=lambda*Curv[i*nps+i];
							}
						}
						B[i]=B0[i];
					}

					// Triangularizing the equations A.delta = B

					for(int i=0;i<nps-1;i++){
						for(int j=i+1;j<nps;j++){
							fac=-A[j*nps+i]/A[i*nps+i];
							for(int m=i+1;m<nps;m++){
								A[j*nps+m]+=A[i*nps+m]*fac;
							}
							B[j]+=B[i]*fac;
						}
					}

					// Solving for delta

					for(int i=nps-1;i>=0;i--){
						fac=B[i];
						for(int j=i+1;j<nps;j++){
							fac-=delta[j]*A[i*nps+j];
						}
						delta[i]=fac/A[i*nps+i];
					}
					// If we end up out of bounds, the point is taken at the bound.
					for(int i=0;i<nps;i++){
	//					printf("%lf ",delta[i]);
						if(!((delta[i]>0)||(delta[i]<0))){
							delta[i]=0.;
						}
						if(fabs(delta[i])>maxdelta[i]){
							delta[i]*=maxdelta[i]/fabs(delta[i]);
						}
						prn[i]=pr[i]+delta[i];
						if(prn[i]>rightlim[i]){
							prn[i]=(0.99*rightlim[i]+(0.01+ilam)*pr[i])/(1+ilam);
						}else if(prn[i]<leftlim[i]){
							prn[i]= (0.99*leftlim[i] + (0.01 + ilam)*pr[i]) / (1 + ilam);
						}
					}
					// Current parameters, current chi square and current lambda are displayed
					(this->*PrintOut)(prn);
					c1=ChiSquared(prn);
					printf("\nilam= %d lambda= %lf\nc1 = %lf prec=%le",ilam,lambda,c1,Tol);

					lambda*=inclambda;
					ilam++;
				}
				lambda/=inclambda;

	
				// if new point is better than previous, pr is updated.
					if(ilam<20){
						for(int i=0;i<nps+2*nfil;i++){
							pr[i]=prn[i];
						}
					}


					// Bumping mechanism: if point is within the covariance ellipsoid it gets kicked off
					fac=0;
					flag=0;
					while(fac<1 && bumperlist){
						for(scanbumper=bumperlist;scanbumper; scanbumper=scanbumper->next){
							fac=scanbumper->distance(pr);
							if(fac<1){
								printf("\nBumped!");
								bumpcounter++;
								for(int i=0;i<nps;i++){
									fac = 2.0*bumperpower / sqrt(fac);
									prn[i]=pr[i]-fac*scanbumper->dp[i];
									if (prn[i]>rightlim[i]) {
										prn[i] = 0.99*rightlim[i] + 0.01*pr[i];
									}
									else if (prn[i]<leftlim[i]) {
										prn[i] = 0.99*leftlim[i] + 0.01*pr[i];
									}
									pr[i] = prn[i];
								}
								scanbumper->UpdateCurvature(bumperpower);
								flag=1;
								(this->*PrintOut)(pr);
							}
						}
						if (bumpcounter >= maxbumpcount) break;
					}
					if(flag){
						c1=ChiSquared(pr);
						printf("\nc1 = %lf",c1);
					}
					// Add new point to stepchain 
					laststep->next=new bumper(pr,nps);
					laststep=laststep->next;
					// Saving to file
					sprintf(filename, "%s-stepchain%d.dat", modelcode ,il);
					f = fopen(filename, "a");
					(this->*PrintFile)(f, c0, false);
					fprintf(f, "\n");
					fclose(f);

					// If new point's chi square is very close to previous one, ichi is increased.
					// When ichi reaches 3, we declare convergence achieved.
					if(fabs(1-c1/c0)<1.0e-3){
						ichi++;
					}else{
						ichi=0;
					}
					// Time limit is set by maxtime
					// Also close if bumpcounter has reached maxbumpcount
					time(&rtime);
					if(difftime(rtime,ltime)>maxtime || bumpcounter>= maxbumpcount) {
						ichi=3;
						nlc=il;
						printf("\n---- Time limit has been reached. Exiting...");
					}
			

			}

			if(ichi<4){
				// Final chi square of this model
				c0=c1;//ChiSquared(pr);
				printf("\nFinal chi square = %lf\n",c0);

				//Check if this is the best model with this initial condition
				if(c0<bestchi) bestchi=c0;

				//Grad(); // Are these two lines really needed? They seem a repetition
				//Covariance();

				// Add a new bumper at the position of this model
				scanbumper=bumperlist;
				bumperlist=new bumper(pr,nps);
				bumperlist->next=scanbumper;
				// Calculate covariance matrix
				Grad();
				Covariance();
				bumperlist->SetCurvature(Curv,np/c0);


				// Saving the steps done to this model
	//			sprintf(filename, "stepchain%d.dat", il);
				//f = fopen(filename, "w");
				//for (scanbumper = stepchain; scanbumper; scanbumper = scanbumper->next) {
				//	fprintf(f, "%.16le", scanbumper->p0[0]);
				//	for (int j = 1; j < nps; j++) {
				//		fprintf(f, " %.16le", scanbumper->p0[j]);
				//	}
				//	fprintf(f, "\n");
				//}
				//fclose(f);

				// Finding the last step outside any bumper
				k=1;
				for(scanbumper=stepchain;scanbumper&&(scanbumper->next);scanbumper=scanbumper->next){
					k++;
					fac=1.e100;
					for(scanbumper2=bumperlist;scanbumper2;scanbumper2=scanbumper2->next){
						fac2=scanbumper2->distance(scanbumper->next->p0);
						if(fac2<fac) fac=fac2;
					}
					printf("\n%d %lf",k,fac);
					if(fac<1.){
						laststep=scanbumper;
						k--;
						scanbumper2=scanbumper->next;
						while(scanbumper2){
							scanbumper=scanbumper2->next;
							delete scanbumper2;
							scanbumper2=scanbumper;
						}
						laststep->next=0;
						scanbumper=laststep;
					}
				}


				// Check if the model makes sense, otherwise set negative chi square for 
				// processing by subsequent programs.
				if((c0>2)&&(c0<1.e100) && flagblending <np/2){
					c0=c0;
				}else{
					c0=c1=-1.;
				}

				// Check if the parameters make sense, otherwise set negative chi square for 
				// processing by subsequent programs.
				for(int i=0;i<nps;i++){
					if(!((pr[i]>-1.e300)&&(pr[i]<1.e300))){
						pr[i]=1.5;
						c0=c1=-1.;
					}
				}



				// Store model in a file named "(il).txt"
				sprintf(filename, "%d.txt",il);
				f=fopen(filename,"w");

				// Print parameters to file

				(this->*PrintFile)(f,c0, true);

				// Print covariance matrix to file
				fprintf(f,"\n");
				for(int i=0;i<nps;i++){
					if(i>0) fprintf(f,"\n");
					fprintf(f,"%le",Cov[i]);
					for(int j=1;j<nps;j++){
						fprintf(f," %le",Cov[i+nps*j]);
					}
				}
				fclose(f);


				// Updating time count
				//printf("\npartial time=%lf secs\n",(Environment::TickCount-tm)/1000.0);
				//printf("total time = %lf mins\n",(Environment::TickCount-tim0)/60000.0);
				//tm=Environment::TickCount;

				// Updating the nlc.dat file
				f=fopen("nlc.dat","w");
				fprintf(f,"%d %d",il,nlc);
				fclose(f);			
			}

			il++;
		}
		
		// When terminating, update the Termination file, 
		// so that subsequent programs know everything has been concluded with this fit
		sprintf(filename,"t%s.dat",outdir);
		f=fopen(filename,"w");
		fprintf(f,"%d %d",il,nlc);
		fclose(f);

		// If this fit has found the best minimum so far, update minchi.dat
		current_path("..");
		strcpy(filename,"minchi.dat");
		if(exists(filename)){
			f=fopen(filename,"r");
			fscanf(f,"%lf",&minchi);
			fclose(f);
		}
		if(bestchi<minchi){
			f=fopen(filename,"w");
			fprintf(f,"%lf %s",bestchi,outdir);
			fclose(f);
		}

		}catch(...){
			error=8;
		}
	}
}

double LevMar::ChiSquared(double *pr){
	double chi2=0,chi0,chi1,p1;
	int fl;

	fl=0;
	sumf[0]=sumfy[0]=sumf2[0]=0;
	for(int i=0;i<np; i++){
		if (filter[i] != fl) {
			fl++;
			sumf[fl] = sumfy[fl] = sumf2[fl] = 0;
		}
		VBBL->satellite=satel[i];
		VBBL->a1 = limbdarks[fl];
		fb[i]=(VBBL->*model)(pr,t[i]);
		sumf[fl]+=w[i]*w[i]*fb[i];
		sumf2[fl]+=w[i]*w[i]*fb[i]*fb[i];
		sumfy[fl]+=w[i]*w[i]*fb[i]*y[i];
	}
	for(int i=0; i<=fl;i++){
#ifndef NOgOGLE
		p1=sumf[i]*sumf[i]-sumf2[i]*sumsigma[i] + epsilon;
		pr[nps+i*2]=(sumf[i]*sumfy[i]-sumf2[i]*sumy[i])/p1;
        pr[nps+1+i*2]=(sumf[i]*sumy[i]-sumsigma[i]*sumfy[i])/p1;
		//chi2+=fabs(sumy2[i]+(sumfy[i]*sumfy[i]*sumsigma[i]+sumf2[i]*sumy[i]*sumy[i]-2*sumf[i]*sumy[i]*sumfy[i])/p1);
#else
		if(i==OGLE){
			pr[nps+i*2]=0.;
		    pr[nps+1+i*2]=sumfy[i]/sumf2[i];
			chi2+=fabs(sumy2[i]-sumfy[i]*pr[nps+1+i*2]);
		}else{
			p1=sumf[i]*sumf[i]-sumf2[i]*sumsigma[i] + epsilon;
			pr[nps+i*2]=(sumf[i]*sumfy[i]-sumf2[i]*sumy[i])/p1;
			pr[nps+1+i*2]=(sumf[i]*sumy[i]-sumsigma[i]*sumfy[i])/p1;
			chi2+=fabs(sumy2[i]+(sumfy[i]*sumfy[i]*sumsigma[i]+sumf2[i]*sumy[i]*sumy[i]-2*sumf[i]*sumy[i]*sumfy[i])/p1);
		}
#endif
		if(pr[nps+1+i*2]<0){
			pr[nps + 1 + i * 2] = 0;
		}

	}
	chi2 = 0.;
	chi0 = 0;
	flagblending = 0;
	for(int i=0;i<np; i++){
		p1=(y[i]-pr[nps+filter[i]*2]-pr[nps+1+filter[i]*2]*fb[i])*w[i]*w[i]*pr[nps+1+filter[i]*2]*Tol;
		chi0+=p1*p1;
		p1 = (y[i] - pr[nps + filter[i] * 2] - pr[nps + 1 + filter[i] * 2] * fb[i]) * w[i];
		chi2 += p1 * p1;
		if (pr[nps + 1 + filter[i] * 2] > 2 * y[i]) {
			flagblending++;
			chi2 += (pr[nps + 1 + filter[i] * 2] - 2 * y[i]) * (pr[nps + 1 + filter[i] * 2] - 2 * y[i]) * w[i] * w[i];
		}
	}
	chi0=sqrt(2*chi0);
	if(chi0/chi2>0.1) Tol*=0.5;
	if(chi0/chi2<0.01 && Tol<.99e-2) Tol*=2;

	return chi2;
}

void LevMar::Grad(){
	double inc,p1;
	int fl;

	for(int j=0;j<nps;j++){
		printf("%d ",j);
		//if((j<2)||(j>3)){
		//	inc=pr[j]*1.0e-3;
		//}else{
			inc=1.0e-3;
		//}
		pr[j]+=inc;	
		fl=0;
		sumf[0]=sumfy[0]=sumf2[0]=0;
		for(int i=0;i<np;i++){
			if (filter[i] != fl) {
				fl++;
				sumf[fl] = sumfy[fl] = sumf2[fl] = 0;
			}
			VBBL->satellite=satel[i];
			VBBL->a1 = limbdarks[fl];
			fb[i+np*(j+1)]=(VBBL->*model)(pr,t[i]);
			sumf[fl]+=w[i]*w[i]*fb[i+np*(j+1)];
			sumf2[fl]+=w[i]*w[i]*fb[i+np*(j+1)]*fb[i+np*(j+1)];
			sumfy[fl]+=w[i]*w[i]*fb[i+np*(j+1)]*y[i];
		}
		pr[j]-=inc;
		for(int i=0; i<=fl;i++){
#ifndef NOgOGLE
			p1=sumf[i]*sumf[i]-sumf2[i]*sumsigma[i] + epsilon;
			prn[i*2]=(sumf[i]*sumfy[i]-sumf2[i]*sumy[i])/p1;
			prn[1+i*2]=(sumf[i]*sumy[i]-sumsigma[i]*sumfy[i])/p1;
#else
			if(i==OGLE){
				prn[i*2]=0.;
				prn[1+i*2]=sumfy[i]/sumf2[i];
			}else{
				p1=sumf[i]*sumf[i]-sumf2[i]*sumsigma[i] + epsilon;
				prn[i*2]=(sumf[i]*sumfy[i]-sumf2[i]*sumy[i])/p1;
				prn[1+i*2]=(sumf[i]*sumy[i]-sumsigma[i]*sumfy[i])/p1;
			}
#endif
			if(prn[1+i*2]<0)  prn[1+i*2]=0;
			dFdp[(1+2*i)*nps+j]=(prn[i*2]+prn[1+i*2]-pr[nps+i*2]-pr[nps+1+i*2])/inc;  // error on baseline FB+FS
			dFdp[(2*i)*nps+j]=(prn[i*2]/prn[1+i*2]-pr[nps+i*2]/pr[nps+1+i*2])/inc;    // error on blending FB/FS
		}
		for(int i=0;i<np;i++){
			Gr[j][i]=w[i]*(prn[filter[i]*2]+prn[1+filter[i]*2]*fb[i+np*(j+1)]-pr[nps+filter[i]*2]-pr[nps+1+filter[i]*2]*fb[i])/inc;
		}
	}
	printf("OK\n");

	for(int j=0;j<nps;j++){
		for(int i=0;i<nps;i++){
			Curv[i*nps+j]=0;
			for(int k=0;k<np;k++){
				Curv[i*nps+j]+=Gr[i][k]*Gr[j][k];
			}
		}
	}
	for(int i=0;i<nps;i++){
		inc=0;
		for(int k=0;k<np;k++){
			inc+=w[k]*Gr[i][k]*(y[k]-pr[nps+filter[k]*2]-pr[nps+1+filter[k]*2]*fb[k]);
		}
		B0[i]=inc;
	}
}

void LevMar::Covariance(){
	//double p1;
	//int j1,k1;
	//p1=Determinant(Curv,nps);
	//for(int i=0;i<nps;i++){
	//	for(int i2=0;i2<nps;i2++){
	//		for(int j=0;j<nps-1;j++){
	//			for(int k=0;k<nps-1;k++){
	//				j1=(j<i2)? j : j+1;
	//				k1=(k<i)? k : k+1;
	//				A[j*(nps-1)+k]=Curv[j1*nps+k1];
	//			}
	//		}
	//		Cov[i*nps+i2]=Determinant(A,nps-1)/p1*(((i+i2)%2)? -1 : 1);
	//	}
	//}
	Inverse(Curv, Cov, nps);
	for (int i = 0; i < nps; i++) {
		errs[i] = sqrt(fabs(Cov[i * nps + i]));
	}
	for(int i=0;i<2*nfil;i++){
		errs[i+nps]=0;
		for(int j1=0;j1<nps;j1++){
			for(int j2=0;j2<nps;j2++){
				errs[i+nps]+=Cov[j1*nps+j2]*dFdp[i*nps+j1]*dFdp[i*nps+j2];
			}
		}
		errs[i+nps]=sqrt(fabs(errs[i+nps]));
		if (!(errs[i + nps] > 0)) errs[i + nps] = 1.e100;
	}
}

void LevMar::PrintOutPS(double *pr){
	printf("\nu0=%lf tE=%lf t0=%lf RS=%lf",exp(pr[0]),exp(pr[1]),pr[2],exp(pr[3]));
}

void LevMar::PrintOutPX(double *pr){
	printf("\nu0=%lf tE=%lf t0=%lf RS=%lf pai1=%lf pai2=%lf",exp(pr[0]),exp(pr[1]),pr[2],exp(pr[3]),pr[4],pr[5]);
}

void LevMar::PrintOutBS(double *pr){
	printf("\ntE=%lf FR=%lf u1=%lf u2=%lf t1=%lf t2=%lf rho=%lf",exp(pr[0]),exp(pr[1]),pr[2],pr[3],pr[4],pr[5],exp(pr[6]));
}

void LevMar::PrintOutBO(double *pr){
	printf("\nu0=%lf t0=%lf tE=%lf rho=%lf xi1=%lf xi2=%lf\nom=%lf inc=%lf phi=%lf qs=%lf",pr[0],pr[1],exp(pr[2]),exp(pr[3]),pr[4],pr[5],pr[6],pr[7],pr[8],exp(pr[9]));
}

void LevMar::PrintOutLS(double *pr){
	printf("\na=%lf q=%lf uc=%lf th=%lf RS=%lf\ntE=%lf tc=%lf",exp(pr[0]),exp(pr[1]),pr[2],pr[3],exp(pr[4]),exp(pr[5]),pr[6]);
}

void LevMar::PrintOutLX(double *pr){
	printf("\na=%lf q=%lf u0=%lf th=%lf RS=%lf\ntE=%lf t0=%lf pai1=%lf pai2=%lf",exp(pr[0]),exp(pr[1]),pr[2],pr[3],exp(pr[4]),exp(pr[5]),pr[6],pr[7],pr[8]);
}

void LevMar::PrintOutLO(double *pr){
	printf("\na=%lf q=%lf u0=%lf th=%lf RS=%lf\ntE=%lf t0=%lf pai1=%lf pai2=%lf\ndsdt=%lf dthdt=%lf w3=%lf",exp(pr[0]),exp(pr[1]),pr[2],pr[3],exp(pr[4]),exp(pr[5]),pr[6],pr[7],pr[8],pr[9],pr[10],pr[11]);
}

void LevMar::PrintOutLK(double* pr) {
	printf("\na=%lf q=%lf u0=%lf th=%lf RS=%lf\ntE=%lf t0=%lf pai1=%lf pai2=%lf\ndsdt=%lf dthdt=%lf w3=%lf szs=%lf ar=%lf", exp(pr[0]), exp(pr[1]), pr[2], pr[3], exp(pr[4]), exp(pr[5]), pr[6], pr[7], pr[8], pr[9], pr[10], pr[11], pr[12], pr[13]);
}


void LevMar::PrintFilePS(FILE *f,double c0, bool printerrors){
	fprintf(f,"%.16le %.16le %.16le %.16le",exp(pr[0]),exp(pr[1]),pr[2],exp(pr[3]));
	for(int i=nps;i<nps+2*nfil;i++){
		pr[i]=((pr[i]>-1.e300)&&(pr[i]<1.e300))? pr[i] : -1.e300;
		fprintf(f," %le",pr[i]);
	}
	fprintf(f," %.16le",c0);
	if (printerrors) {
		fprintf(f, "\n%le %le %le %le", errs[0] * exp(pr[0]), errs[1] * exp(pr[1]), errs[2], errs[3] * exp(pr[3]));
		for (int i = nps; i < nps + 2 * nfil; i++) {
			fprintf(f, " %le", errs[i]);
		}
	}
}

void LevMar::PrintFilePX(FILE *f,double c0, bool printerrors){
	fprintf(f,"%.16le %.16le %.16le %.16le %.16le %.16le",pr[0],exp(pr[1]),pr[2],exp(pr[3]),pr[4],pr[5]);
	for(int i=nps;i<nps+2*nfil;i++){
		pr[i]=((pr[i]>-1.e300)&&(pr[i]<1.e300))? pr[i] : -1.e300;
		fprintf(f," %le",pr[i]);
	}
	fprintf(f," %.16le",c0);
	if (printerrors) {
		fprintf(f, "\n%le %le %le %le %le %le", errs[0], errs[1] * exp(pr[1]), errs[2], errs[3] * exp(pr[3]), errs[4], errs[5]);
		for (int i = nps; i < nps + 2 * nfil; i++) {
			fprintf(f, " %le", errs[i]);
		}
	}
}

void LevMar::PrintFileBS(FILE *f,double c0, bool printerrors){
	if (pr[1] > 0) { // Invert sources
		double sc;
		pr[1] = -pr[1];
		sc = pr[3];
		pr[3] = pr[2];
		pr[2] = sc;
		sc = pr[5];
		pr[5] = pr[4];
		pr[4] = sc;
		sc = errs[3];
		errs[3] = errs[2];
		errs[2] = sc;
		sc = errs[5];
		errs[5] = errs[4];
		errs[4] = sc;
		for (int k = 0; k < nps; k++) {
			Cov[1 + nps * k] = - Cov[1 + nps * k];
			Cov[k + nps * 1] = - Cov[k + nps * 1];
		}
		for (int k = 0; k < nps; k++) {
			sc = Cov[2 + nps * k];
			Cov[2 + nps * k] = Cov[3 + nps * k];
			Cov[3 + nps * k] = sc;
			sc = Cov[5 + nps * k];
			Cov[5 + nps * k] = Cov[4 + nps * k];
			Cov[4 + nps * k] = sc;
		}
		for (int k = 0; k < nps; k++) {
			sc = Cov[k + nps * 2];
			Cov[k + nps * 2] = Cov[k + nps * 3];
			Cov[k + nps * 3] = sc;
			sc = Cov[k + nps * 5];
			Cov[k + nps * 5] = Cov[k + nps * 4];
			Cov[k + nps * 4] = sc;
		}
	}

	fprintf(f,"%.16le %.16le %.16le %.16le %.16le %.16le %.16le",exp(pr[0]),exp(pr[1]),pr[2],pr[3],pr[4],pr[5],exp(pr[6]));
	for(int i=nps;i<nps+2*nfil;i++){
		pr[i]=((pr[i]>-1.e300)&&(pr[i]<1.e300))? pr[i] : -1.e300;
		fprintf(f," %le",pr[i]);
	}
	fprintf(f," %.16le",c0);
	if (printerrors) {
		fprintf(f, "\n%le %le %le %le %le %le %le", errs[0] * exp(pr[0]), errs[1] * exp(pr[1]), errs[2], errs[3], errs[4], errs[5], errs[6] * exp(pr[6]));
		for (int i = nps; i < nps + 2 * nfil; i++) {
			fprintf(f, " %le", errs[i]);
		}
	}
}

void LevMar::PrintFileBO(FILE *f,double c0, bool printerrors){
	fprintf(f,"%.16lf %.16lf %.16lf %.16lf %.16lf %.16lf %.16lf %.16lf %.16lf %.16lf", pr[0], pr[1], exp(pr[2]), exp(pr[3]), pr[4], pr[5], pr[6], pr[7], pr[8], exp(pr[9]));
	for(int i=nps;i<nps+2*nfil;i++){
		pr[i]=((pr[i]>-1.e300)&&(pr[i]<1.e300))? pr[i] : -1.e300;
		fprintf(f," %le",pr[i]);
	}
	fprintf(f," %.16le",c0);
	if (printerrors) {
		fprintf(f, "\n%le %le %le %le %le %le %le %le %le %le", errs[0], errs[1], errs[2] * exp(pr[2]), errs[3] * exp(pr[3]), errs[4], errs[5], errs[6], errs[7], errs[8], errs[9] * exp(pr[9]));
		for (int i = nps; i < nps + 2 * nfil; i++) {
			fprintf(f, " %le", errs[i]);
		}
	}
}

void LevMar::PrintFileLS(FILE *f,double c0, bool printerrors){
	if (pr[1] > 0) {
		pr[3] = pr[3] - M_PI;
		pr[1] = -pr[1];
		for (int k = 0; k < nps; k++) {
			Cov[1 + nps * k] = -Cov[1 + nps * k];
			Cov[k + nps * 1] = -Cov[k + nps * 1];
		}
	}
	while (pr[3] > 2*M_PI) pr[3] -= 2*M_PI;
	while (pr[3] < 0) pr[3] += 2*M_PI;
	if (pr[2] < 0) {
		pr[3] = 2 * M_PI - pr[3];
		pr[2] = -pr[2];
		for (int k = 0; k < nps; k++) {
			Cov[2 + nps * k] = -Cov[2 + nps * k];
			Cov[k + nps * 2] = -Cov[k + nps * 2];
			Cov[3 + nps * k] = -Cov[3 + nps * k];
			Cov[k + nps * 3] = -Cov[k + nps * 3];
		}
	}

	fprintf(f,"%.16le %.16le %.16le %.16le %.16le %.16le %.16le",exp(pr[0]),exp(pr[1]),pr[2],pr[3],exp(pr[4]),exp(pr[5]),pr[6]);
	for(int i=nps;i<nps+2*nfil;i++){
		pr[i]=((pr[i]>-1.e300)&&(pr[i]<1.e300))? pr[i] : -1.e300;
		fprintf(f," %le",pr[i]);
	}
	fprintf(f," %.16le",c0);
	if (printerrors) {
		fprintf(f, "\n%le %le %le %le %le %le %le", errs[0] * exp(pr[0]), errs[1] * exp(pr[1]), errs[2], errs[3], errs[4] * exp(pr[4]), errs[5] * exp(pr[5]), errs[6]);
		for (int i = nps; i < nps + 2 * nfil; i++) {
			fprintf(f, " %le", errs[i]);
		}
	}
}

void LevMar::PrintFileLX(FILE *f,double c0, bool printerrors){
	if (pr[1] > 0) {
		pr[3] = pr[3] - M_PI;
		pr[1] = -pr[1];
		for (int k = 0; k < nps; k++) {
			Cov[1 + nps * k] = -Cov[1 + nps * k];
			Cov[k + nps * 1] = -Cov[k + nps * 1];
		}
	}
	while (pr[3] > 2 * M_PI) pr[3] -= 2 * M_PI;
	while (pr[3] < 0) pr[3] += 2 * M_PI;

	fprintf(f,"%.16le %.16le %.16le %.16le %.16le %.16le %.16le %.16le %.16le",exp(pr[0]),exp(pr[1]),pr[2],pr[3],exp(pr[4]),exp(pr[5]),pr[6],pr[7],pr[8]);
	for(int i=nps;i<nps+2*nfil;i++){
		pr[i]=((pr[i]>-1.e300)&&(pr[i]<1.e300))? pr[i] : -1.e300;
		fprintf(f," %le",pr[i]);
	}
	fprintf(f," %.16le",c0);
	if (printerrors) {
		fprintf(f, "\n%le %le %le %le %le %le %le %le %le", errs[0] * exp(pr[0]), errs[1] * exp(pr[1]), errs[2], errs[3], errs[4] * exp(pr[4]), errs[5] * exp(pr[5]), errs[6], errs[7], errs[8]);
		for (int i = nps; i < nps + 2 * nfil; i++) {
			fprintf(f, " %le", errs[i]);
		}
	}
}

void LevMar::PrintFileLO(FILE *f,double c0, bool printerrors){
	if (pr[1] > 0) {
		pr[3] = pr[3] - M_PI;
		pr[1] = -pr[1];
		for (int k = 0; k < nps; k++) {
			Cov[1 + nps * k] = -Cov[1 + nps * k];
			Cov[k + nps * 1] = -Cov[k + nps * 1];
		}
	}
	while (pr[3] > 2 * M_PI) pr[3] -= 2 * M_PI;
	while (pr[3] < 0) pr[3] += 2 * M_PI;

	fprintf(f,"%.16le %.16le %.16le %.16le %.16le %.16le %.16le %.16le %.16le %.16le %.16le %.16le",exp(pr[0]),exp(pr[1]),pr[2],pr[3],exp(pr[4]),exp(pr[5]),pr[6],pr[7],pr[8],pr[9],pr[10],pr[11]);
	for(int i=nps;i<nps+2*nfil;i++){
		pr[i]=((pr[i]>-1.e300)&&(pr[i]<1.e300))? pr[i] : -1.e300;
		fprintf(f," %le",pr[i]);
	}
	fprintf(f," %.16le",c0);
	if (printerrors) {
		fprintf(f, "\n%le %le %le %le %le %le %le %le %le %le %le %le", errs[0] * exp(pr[0]), errs[1] * exp(pr[1]), errs[2], errs[3], errs[4] * exp(pr[4]), errs[5] * exp(pr[5]), errs[6], errs[7], errs[8], errs[9], errs[10], errs[11]);
		for (int i = nps; i < nps + 2 * nfil; i++) {
			fprintf(f, " %le", errs[i]);
		}
	}
}


void LevMar::PrintFileLK(FILE* f, double c0, bool printerrors) {
	fprintf(f, "%.16le %.16le %.16le %.16le %.16le %.16le %.16le %.16le %.16le %.16le %.16le %.16le %.16le %.16le", exp(pr[0]), exp(pr[1]), pr[2], pr[3], exp(pr[4]), exp(pr[5]), pr[6], pr[7], pr[8], pr[9], pr[10], pr[11],pr[12],pr[13]);
	for (int i = nps; i < nps + 2 * nfil; i++) {
		pr[i] = ((pr[i] > -1.e300) && (pr[i] < 1.e300)) ? pr[i] : -1.e300;
		fprintf(f, " %le", pr[i]);
	}
	fprintf(f, " %.16le", c0);
	if (printerrors) {
		fprintf(f, "\n%le %le %le %le %le %le %le %le %le %le %le %le %le %le", errs[0] * exp(pr[0]), errs[1] * exp(pr[1]), errs[2], errs[3], errs[4] * exp(pr[4]), errs[5] * exp(pr[5]), errs[6], errs[7], errs[8], errs[9], errs[10], errs[11], errs[12], errs[13]);
		for (int i = nps; i < nps + 2 * nfil; i++) {
			fprintf(f, " %le", errs[i]);
		}
	}
}
