// LevMarFit.cpp 
// Implementation of all methods in LevMarFit.h for Levenberg-Marquardt fitting

#define _CRT_SECURE_NO_WARNINGS
#define _USE_MATH_DEFINES
#include "LevMarFit.h"
#include "bumper.h"
#include <VBMicrolensingLibrary.h>
#include <cstdio>
#include <ctime>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <regex>
#include <filesystem>

//using namespace std;
using std::regex, std::string, std::regex_match;
using namespace std::filesystem;

int nlc = 6; // Number of models to be calculated from the same initial condition using the bumper method
int offsetdegeneracy = 3; // Number of models to be calculated with offset degeneracy
int maxsteps = 50; // Maximum number of steps in each fit
double maxtime = 1.e100; // 600.0; // Maximum time in seconds for total execution (no longer controlled within LevMar)
double bumperpower = 2.0; // Repulsion factor of bumpers
double maxbumpcount = 25;
char parametersfile[256] = ""; // File from which user parameters are read, if any
double mass_luminosity_exponent = 4.0; // Mass-luminosity exponent for the sources
double mass_radius_exponent = 0.9; // Mass-radius exponent for the sources
double lens_mass_luminosity_exponent = 4.0; // Mass-luminosity exponent for the lenses
bool stepchainsave = false;

std::vector<std::vector<std::string>> parnames = { {"u0","tE","t0","rho"},
					{"u0","tE","t0","rho","piN","piE"},
					{"tE","FR","u01","u02","t0","t02","rho"},
					{"tE","FR","u01","u02","t0","t02","rho","piN","piE","gamma1","gamma2","gamma3"},
					{"s","q","u0","alpha","rho","tE","t0"},
					{"s","q","u0","alpha","rho","tE","t0","piN","piE"},
					{"s","q","u0","alpha","rho","tE","t0","piN","piE","gamma1","gamma2","gammaz"},
					{"s","q","u0","alpha","rho","tE","t0","piN","piE","gamma1","gamma2","gammaz","sz_s","a_s3d"},
					{"s","q","u0","alpha","rho","tE","t0","s2","q2","beta"},
					{"s","q","u0","alpha","rho","tE","t0","s2","q2","beta","piN","piE"} };
std::vector<std::string> astroparnames = { "muSDec", "muSRA", "piS", "thetaE" };

std::vector<std::vector<int>> logposs = { {0, 1, 3},
										{1, 3},
										{0, 1, 6},
										{0, 1, 6},
										{0, 1, 4, 5},
										{0, 1, 4, 5},
										{0, 1, 4, 5},
										{0, 1, 4, 5},
										{0, 1, 4, 5, 7, 8},
										{0, 1, 4, 5, 7, 8} };

const double epsilon = 1.e-100;

LevMar::LevMar(int argc, char* argv[]) {
	setbuf(stdout, nullptr);
	printf("******************************************\n");
	printf("*************      LevMar     **********\n");
	printf("******************************************\n\n\n");
	printf("This program fits an event from specific initial conditions\n\n");
	error = 0;
	it0 = it02 = -1;
	astrometric = false;
	nlinpar = 2;
	// Setting default values
	Tol = 1.e-2;
	VBM = new VBMicrolensing;
	VBM->Tol = Tol;
	VBM->RelTol = 0.001;
	VBM->parallaxsystem = 1;
	VBM->SetMethod(VBMicrolensing::Method::Multipoly);

	ReadFiles(argc, argv);

}

LevMar::~LevMar() {
	delete VBM;
	if (error < 9) {
		free(t);
		free(y);
		free(w);
		free(y1a);
		free(y2a);
		free(seps);
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
		free(starts);
		free(sizes);
		free(dFdp);
		free(cN);
		free(cE);
		free(wcN);
		free(wcE);
		free(c1l);
		free(c2l);
		free(c1s);
		free(c2s);
		if (astrometric) {

			free(sumsigmaN);
			free(sumsigmaE);
			free(sumcN);
			free(sumcE);
			free(sumc1);
			free(sumc2);
		}

		for (int i = 0; i < nps; i++) {
			free(Gr[i]);
		}
		free(Gr);
		if (consnumber > 0) {
			free(consindex);
			free(constraints);
			free(consleft);
			free(consright);
			free(consvars);
		}


		bumper* scanbumper = bumperlist, * scanbumper2;
		while (scanbumper) {
			scanbumper2 = scanbumper->next;
			delete scanbumper;
			scanbumper = scanbumper2;
		}


		scanbumper = stepchain;
		while (scanbumper) {
			scanbumper2 = scanbumper->next;
			delete scanbumper;
			scanbumper = scanbumper2;
		}

	}
	if (error < 9) {
		free(sigmapr);
		free(leftlim);
		free(rightlim);
	}
}

void LevMar::ReadFiles(int argc, char* argv[]) {

	try {

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

		current_path(eventname);
		current_path("Data");

		/* Reading coordinates */

		auto searchstring = regex(".*\\.coordinates");
		for (auto const& itr : directory_iterator(".")) {
			string curfile = (itr).path().filename().string();
			if (regex_match(curfile, searchstring)) {
				VBM->SetObjectCoordinates((char*)curfile.c_str(), satellitedir);
				printf("\n- Coordinates set.");
				break;
			}
		}

		current_path(eventname);

		// Read light curve
		ReadCurve();
		if (error) return;

		// Establishing model to be fit
		nps = (astrometric) ? 4 : 0;
		strcpy(modelcode, outdir);
		switch (modelcode[0]) {
		case 'P':
			if (modelcode[1] == 'X') {
				modnumber = 1;
				nps += 6;
				double presigmapr[] = { .5,1.0,5.,2.3,0.1,0.1,         1.0, 1.0, 0.1, 0.2 };
				double preleftlim[] = { -3.0,-6.9,-10.e100,-11.5,-10.,-10.,       -30.0, -30.0, 0.05, 0.001 };
				double prerightlim[] = { 3.0,6.9,10.e100,0.0,10.,10.,          30.0, 30.0, 1.0, 30.0 };
				ReadOptions(preleftlim, prerightlim, presigmapr);
				it0 = 2;
				error = InitCond(presigmapr, preleftlim, prerightlim);
				pr[1] = log(pr[1]);
				pr[3] = log(pr[3]);
				current_path(exedir);
				//current_path("..");
				current_path("data");
				VBM->LoadSunTable("SunEphemeris.txt");
				current_path(eventname);
			}

			else {
				modnumber = 0;
				nps = 4;
				double presigmapr[] = { 1.15,1.0,5.,2.3 };
				double preleftlim[] = { -11.0,-6.9,-10.e100,-11.5 };
				double prerightlim[] = { 1.0,6.9,10.e100,0.0 };
				ReadOptions(preleftlim, prerightlim, presigmapr);
				it0 = 2;
				error = InitCond(presigmapr, preleftlim, prerightlim);
				pr[0] = log(pr[0]);
				pr[1] = log(pr[1]);
				pr[3] = log(pr[3]);
			}
			current_path(exedir);
			//current_path("..");
			current_path("data");
			VBM->LoadESPLTable("ESPL.tbl");
			current_path(eventname);
			break;

		case 'B':
			if (modelcode[1] == 'O') {
				modnumber = 3;
				nps += 12;
				double presigmapr[] = { 1.0,.5,0.5,0.5,1,1,0.5,0.03,0.03,0.01,0.01,0.01,          1.0, 1.0, 0.1, 0.2 };
				double preleftlim[] = { -6.9,-11.5,0.,0.,-10.e100,-10.e100,-11.5,-3.,-3.,-1,-1,1.e-7,       -30.0, -30.0, 0.05, 0.001 };
				double prerightlim[] = { 6.9,0.0,3.,3.,10.e100,10.e100,0.0,3.,3.,1.,1.,1.,          30.0, 30.0, 1.0, 30.0 };
				ReadOptions(preleftlim, prerightlim, presigmapr);
				it0 = 4;
				it02 = 5;
				error = InitCond(presigmapr, preleftlim, prerightlim);
				pr[0] = log(pr[0]);
				pr[1] = log(pr[1]);
				pr[6] = log(pr[6]);

				current_path(exedir);
				//current_path("..");
				current_path("data");
				VBM->LoadSunTable("SunEphemeris.txt");
				current_path(eventname);
			}
			else {
				modnumber = 2;
				nps = 7;
				double presigmapr[] = { 1.0,.5,0.5,0.5,1,1,0.5 };
				double preleftlim[] = { -6.9,-11.5,0.,0.,-10.e100,-10.e100,-11.5 };
				double prerightlim[] = { 6.9,0.0,3.,3.,10.e100,10.e100,0.0 };
				ReadOptions(preleftlim, prerightlim, presigmapr);
				it0 = 4;
				it02 = 5;
				error = InitCond(presigmapr, preleftlim, prerightlim);
				pr[0] = log(pr[0]);
				pr[1] = log(pr[1]);
				pr[6] = log(pr[6]);
			}
			current_path(exedir);
			//current_path("..");
			current_path("data");
			VBM->LoadESPLTable("ESPL.tbl");
			current_path(eventname);
			break;
		case 'L':
			if (modelcode[1] == 'X') {
				modnumber = 5;
				nps += 9;
				double presigmapr[] = { .1,0.5,.1,.1,0.3,.6,5.,0.03,0.03,          1.0, 1.0, 0.1, 0.2 };
				double preleftlim[] = { -4.0,-16.1,-3.,-12.56,-11.5,-6.9,-10.e100,-3.,-3.,       -30.0, -30.0, 0.05, 0.001 };
				double prerightlim[] = { 3.0,16.1,3.,12.56,-2.5,7.6,10.e100,3.,3.,          30.0, 30.0, 1.0, 30.0 };
				ReadOptions(preleftlim, prerightlim, presigmapr);
				it0 = 6;
				error = InitCond(presigmapr, preleftlim, prerightlim);
				pr[0] = log(pr[0]);
				pr[1] = log(pr[1]);
				pr[4] = log(pr[4]);
				pr[5] = log(pr[5]);
				current_path(exedir);
				current_path("..");
				current_path("data");
				VBM->LoadSunTable("SunEphemeris.txt");
				current_path(eventname);
			}
			else {
				if (modelcode[1] == 'O') {
					modnumber = 6;
					nps += 12;
					double presigmapr[] = { .1,0.5,.1,.1,0.3,.6,5.,0.03,0.03,0.01,0.01,0.01,          1.0, 1.0, 0.1, 0.2 };
					double preleftlim[] = { -4.0,-16.1,-3.,-12.56,-11.5,-6.9,-10.e100,-3.,-3.,-1,-1,1.e-7,       -30.0, -30.0, 0.05, 0.001 };
					double prerightlim[] = { 3.0,16.1,3.,12.56,-2.5,7.6,10.e100,3.,3.,1,1,1,          30.0, 30.0, 1.0, 30.0 };
					ReadOptions(preleftlim, prerightlim, presigmapr);
					it0 = 6;
					error = InitCond(presigmapr, preleftlim, prerightlim);
					pr[0] = log(pr[0]);
					pr[1] = log(pr[1]);
					pr[4] = log(pr[4]);
					pr[5] = log(pr[5]);
					current_path(exedir);
					//current_path("..");
					current_path("data");
					VBM->LoadSunTable("SunEphemeris.txt");
					current_path(eventname);
				}
				else {
					if (modelcode[1] == 'K') {
						modnumber = 7;
						nps += 14;
						double presigmapr[] = { .1,0.5,.1,.1,0.3,.6,5.,0.03,0.03,0.01,0.01,0.01,0.1, 0.1,          1.0, 1.0, 0.1, 0.2 };
						double preleftlim[] = { -4.0,-16.1,-3.,-12.56,-11.5,-6.9,-10.e100,-3.,-3.,-1,-1,1.e-7, -10,0.5001,       -30.0, -30.0, 0.05, 0.001 };
						double prerightlim[] = { 3.0,16.1,3.,12.56,-2.5,7.6,10.e100,3.,3.,1,1,1,10,10,          30.0, 30.0, 1.0, 30.0 };
						ReadOptions(preleftlim, prerightlim, presigmapr);
						it0 = 6;
						error = InitCond(presigmapr, preleftlim, prerightlim);
						pr[0] = log(pr[0]);
						pr[1] = log(pr[1]);
						pr[4] = log(pr[4]);
						pr[5] = log(pr[5]);
						current_path(exedir);
						//current_path("..");
						current_path("data");
						VBM->LoadSunTable("SunEphemeris.txt");
						current_path(eventname);
					}
					else {
						modnumber = 4;
						nps = 7;
						double presigmapr[] = { .1,0.5,.1,.1,0.3,.6,5. };
						double preleftlim[] = { -4.0,-16.1,-3.,-12.56,-11.5,-6.9,-10.e100 };
						double prerightlim[] = { 3.0,16.1,3.,12.56,-2.5,7.6,10.e100 };
						ReadOptions(preleftlim, prerightlim, presigmapr);
						it0 = 6;
						error = InitCond(presigmapr, preleftlim, prerightlim);
						pr[0] = log(pr[0]);
						pr[1] = log(pr[1]);
						pr[4] = log(pr[4]);
						pr[5] = log(pr[5]);
					}
				}
			}
			break;
		case 'T':
			if (modelcode[1] == 'X') {
				modnumber = 9;
				nps += 12;
				double presigmapr[] = { .1,.4,.1,.1,4.6,.1,1.,.1,.4,.1, 1.,1.,          1.0, 1.0, 0.1, 0.2 };
				double preleftlim[] = { -4.0,-11.5,-3.,-12.56,-11.5,-6.9,-10.e100,-4.0,-11.5,-12.56, -3.,-3.,       -30.0, -30.0, 0.05, 0.001 };
				double prerightlim[] = { 3.0,11.5,3.,12.56,-2.5,7.6,10.e100,3.0, 11.5, 12.56, 3.,3.,          30.0, 30.0, 1.0, 30.0 };
				ReadOptions(preleftlim, prerightlim, presigmapr);
				it0 = 6;
				error = InitCond(presigmapr, preleftlim, prerightlim);
				pr[0] = log(pr[0]);
				pr[1] = log(pr[1]);
				pr[4] = log(pr[4]);
				pr[5] = log(pr[5]);
				pr[7] = log(pr[7]);
				pr[8] = log(pr[8]);
				current_path(exedir);
				//current_path("..");
				current_path("data");
				VBM->LoadSunTable("SunEphemeris.txt");
				current_path(eventname);
			}
			else {
				modnumber = 8;
				nps = 10;
				double presigmapr[] = { .1,0.5,.1,.1,0.3,.6,5., 0.3, 0.5, 0.3 };
				double preleftlim[] = { -4.0,-11.5,-3.,-12.56,-11.5,-6.9,-10.e100,-4.0,-11.5,-12.56 };
				double prerightlim[] = { 3.0,11.5,3.,12.56,-2.5,7.6,10.e100,3.0, 11.5, 12.56 };
				ReadOptions(preleftlim, prerightlim, presigmapr);
				it0 = 6;
				error = InitCond(presigmapr, preleftlim, prerightlim);
				pr[0] = log(pr[0]);
				pr[1] = log(pr[1]);
				pr[4] = log(pr[4]);
				pr[5] = log(pr[5]);
				pr[7] = log(pr[7]);
				pr[8] = log(pr[8]);
			}
			break;
		}
		printf("\n- Model: %s", modelcode);

	}
	catch (...) {
		error = 10;
	}
	current_path(eventname);

	if (!error) {
		// Reading Light Curve
		ReadAncillary();
	}

}

void LevMar::ReadCurve() {
	FILE* f;

	try {
		printf("\n\n- Reading data \n");

		f = fopen("LCToFit.txt", "r");
		fscanf(f, "%d", &np);
		filter = (int*)malloc(sizeof(int) * np);
		t = (double*)malloc(sizeof(double) * np);
		y = (double*)malloc(sizeof(double) * np);
		w = (double*)malloc(sizeof(double) * np);
		y1a = (double*)malloc(sizeof(double) * np);
		y2a = (double*)malloc(sizeof(double) * np);
		satel = (int*)malloc(sizeof(int) * np);
		cN = (double*)malloc(sizeof(double) * np);
		cE = (double*)malloc(sizeof(double) * np);
		wcN = (double*)malloc(sizeof(double) * np);
		wcE = (double*)malloc(sizeof(double) * np);


		nfil = 1;
		for (int i = 0; i < np; i++) {
			fscanf(f, "%d %lg %lg %lg %d %lg %lg %lg %lg", &(filter[i]), &(t[i]), &(y[i]), &(w[i]), &(satel[i]), &(cN[i]), &(wcN[i]), &(cE[i]), &(wcE[i]));
			if (wcN[i] > 0) astrometric = true; // Check if astrometric set.
			if (satel[i] > VBM->nsat) {
				printf("\n! Satellite table %d not found in satellite dir \"%s\"", satel[i], satellitedir);
				error = 9;
				break;
			}
			if ((i != 0) && (filter[i] > filter[i - 1])) {
				nfil++;
			}
			w[i] = 1 / (w[i]);
			wcN[i] = 1 / (wcN[i]);
			wcE[i] = 1 / (wcE[i]);
		}
		fclose(f);
		if (astrometric) {
			nlinpar = 4;
			VBM->astrometry = true;
		}

	}
	catch (...) {
		error = 9;
	}

}

void LevMar::ReadOptions(double* preleftlim, double* prerightlim, double* presigmapr) {
	char buffer[3200];
	char command[256];
	char value[256], value2[256], value3[256];

	consnumber = 0;
	FILE* f;
	/* Reading ini files */

	if (exists("ini")) {
		current_path("ini");

		/* Reading parameters ranges */

		f = fopen("Parameters_Ranges.ini", "r");
		if (f != 0) {
			strcpy(buffer, modelcode);
			buffer[2] = 0;
			while (!feof(f)) {
				fscanf(f, "%s", command);
				if (strcmp(buffer, command) == 0) {
					int npp = (astrometric) ? nps - 4 : nps;
					for (int i = 0; i < npp; i++) {
						fscanf(f, "%lg %lg %lg", &preleftlim[i], &prerightlim[i], &presigmapr[i]);
					}
				}
				if (astrometric && strcmp("astrometry", command) == 0) {
					for (int i = nps - 4; i < nps; i++) {
						fscanf(f, "%lg %lg %lg", &preleftlim[i], &prerightlim[i], &presigmapr[i]);
					}
				}
			}
			fclose(f);
		}
		/* Reading constraints */

		f = fopen("Constraints.ini", "r");
		if (f != 0) {
			while (!feof(f)) {
				fscanf(f, "%[^\n]s", command);
				fscanf(f, "%s", command);
				consnumber++;
			}
			fclose(f);
			consindex = (int*)malloc(sizeof(int) * consnumber);
			constraints = (double*)malloc(sizeof(double) * consnumber);
			consleft = (double*)malloc(sizeof(double) * consnumber);
			consright = (double*)malloc(sizeof(double) * consnumber);
			consvars = (double*)malloc(sizeof(double) * consnumber * (nps + 1));
			printf("\n\n- Reading Constraints.ini");
			f = fopen("Constraints.ini", "r");
			int conscurrent = 0;
			while (!feof(f)) {
				int red = fscanf(f, "%s %s %s %s %s", command, buffer, value, value2, value3);
				printf("\n%s %s %s %s %s", command, buffer, value, value2, value3);
				if (red < 5) {
					command[0] = 0;
					//if (red != 0) {
				}
				else {
					consindex[conscurrent] = -1;
					strcpy(buffer, command);
					int flaglog = 0, flaglog2 = 0;
					buffer[4] = 0;
					for (int i = 0; i < (int)parnames[modnumber].size(); i++) {
						if (strcmp(buffer, "log_") == 0) { // Log parameter
							if (strcmp(parnames[modnumber][i].c_str(), &command[4]) == 0) {
								consindex[conscurrent] = i;
								flaglog = 1;
							}
						}
						else {
							if (strcmp(parnames[modnumber][i].c_str(), command) == 0) {
								consindex[conscurrent] = i;
							}
						}
					}
					for (int j = 0; j < (int)logposs[modnumber].size(); j++) {
						if (logposs[modnumber][j] == consindex[conscurrent]) flaglog2 = 1;
					}

					for (int i = 0; i < (int)astroparnames.size(); i++) { // Astrometric parameters
						if (strcmp(buffer, "log_") == 0) { // Log parameter
							if (strcmp(astroparnames[i].c_str(), &command[4]) == 0) {
								consindex[conscurrent] = 100 + i;
								flaglog = 1;
							}
						}
						else {
							if (strcmp(astroparnames[i].c_str(), command) == 0) {
								consindex[conscurrent] = 100 + i;
							}
						}
					}

					if (command[0] == 'g' && command[1] == '_') { // Blending
						sscanf(&command[2], "%d", &consindex[conscurrent]);
						consindex[conscurrent] += 10000;
					}
					// special combinations
					if (strcmp(command, "muangle") == 0) {
						for (int i = 0; i < nps; i++) {
							if (strcmp(parnames[modnumber][i].c_str(), "piN") == 0) consindex[conscurrent] = 30000;
						}
					}
					if (strcmp(command, "t*") == 0) {
						consindex[conscurrent] = 30001;
					}
					if (consindex[conscurrent] >= 0) {
						sscanf(value, "%lg", &constraints[conscurrent]);
						sscanf(value2, "%lg", &consleft[conscurrent]);
						sscanf(value3, "%lg", &consright[conscurrent]);
						consleft[conscurrent] = fabs(consleft[conscurrent]);
						consright[conscurrent] = fabs(consright[conscurrent]);
						if (flaglog == 0 && flaglog2 == 1) {
							consleft[conscurrent] /= constraints[conscurrent];
							consright[conscurrent] /= constraints[conscurrent];
							constraints[conscurrent] = log(constraints[conscurrent]);
						}
						if (flaglog == 1 && flaglog2 == 0) {
							constraints[conscurrent] = exp(constraints[conscurrent] * log(10));
							consleft[conscurrent] *= constraints[conscurrent] * log(10);
							consright[conscurrent] *= constraints[conscurrent] * log(10);
						}
						if (flaglog == 1 && flaglog2 == 1) {
							constraints[conscurrent] *= log(10);
							consleft[conscurrent] *= log(10);
							consright[conscurrent] *= log(10);
						}
						conscurrent++;
					}
				}
			}
			fclose(f);
			consnumber = conscurrent;
			printf("\nConstraints: %d", consnumber);
		}

		/* Reading options */

		f = fopen("LevMar.ini", "r");
		if (f != 0) {
			printf("\n\n- Reading options in LevMar.ini");
			while (!feof(f)) {
				int red = fscanf(f, "%s %s %s", command, buffer, value);
				if (red < 1) {
					command[0] = 0;
					//if (red != 0) {
					//	printf("\n\n!!! Bad command in Reader.ini");
					//	return -1;
					//};
				}
				if (strcmp(command, "nfits") == 0) {
					sscanf(value, "%d", &nlc);
				}
				if (strcmp(command, "offsetdegeneracy") == 0) {
					sscanf(value, "%d", &offsetdegeneracy);
				}
				if (strcmp(command, "maxsteps") == 0) {
					sscanf(value, "%d", &maxsteps);
				}
				if (strcmp(command, "timelimit") == 0) {
					// sscanf(value, "%lf", &maxtime);
					maxtime = maxtime;  // No longer controlled within LevMar
				}
				if (strcmp(command, "bumperpower") == 0) {
					sscanf(value, "%lg", &bumperpower);
				}
				if (strcmp(command, "parametersfile") == 0) {
					strcpy(parametersfile, value);
				}
				if (strcmp(command, "mass_luminosity_exponent") == 0) {
					sscanf(value, "%lg", &mass_luminosity_exponent);
				}
				if (strcmp(command, "mass_radius_exponent") == 0) {
					sscanf(value, "%lg", &mass_radius_exponent);
				}
				if (strcmp(command, "lens_mass_luminosity_exponent") == 0) {
					sscanf(value, "%lg", &lens_mass_luminosity_exponent);
				}
				if (strcmp(command, "turn_off_secondary_lens") == 0 && strcmp(value, "True") == 0) {
					VBM->turn_off_secondary_lens = true;
				}
				if (strcmp(command, "turn_off_secondary_source") == 0 && strcmp(value, "True") == 0) {
					VBM->turn_off_secondary_source = true;
				}
				if (strcmp(command, "stepchainsave") == 0 && strcmp(value, "True") == 0) {
					stepchainsave = true;
				}
			}
			fclose(f);
		}
		else {
			printf("\n\n- Default options:");
		}
	}
	VBM->mass_luminosity_exponent = mass_luminosity_exponent;
	VBM->mass_radius_exponent = mass_radius_exponent;
	VBM->lens_mass_luminosity_exponent = lens_mass_luminosity_exponent;
	current_path(eventname);
}

int LevMar::InitCond(double* presigmapr, double* preleftlim, double* prerightlim) {
	char buffer[3200], initcondfile[256];
	int npeaks, ninit, incond;
	FILE* f;
	sigmapr = (double*)malloc(sizeof(double) * nps);
	leftlim = (double*)malloc(sizeof(double) * nps);
	rightlim = (double*)malloc(sizeof(double) * nps);
	pr = (double*)malloc(sizeof(double) * nps);
	if (parametersfile[0] == 0) {
		sprintf(initcondfile, "InitCond%c%c.txt", modelcode[0], modelcode[1]);
		current_path("InitCond");
		f = fopen(initcondfile, "r");
		fscanf(f, "%d %d", &npeaks, &ninit);
		incond = atoi(outdir + 2);
		if (incond >= ninit) return 10;
		fgets(buffer, 3200, f);
		for (int i = 0; i < npeaks + incond; i++) {
			fgets(buffer, 3200, f);
		}
	}
	else {
		f = fopen(parametersfile, "r");
	}
	for (int i = 0; i < nps; i++) {
		sigmapr[i] = presigmapr[i];
		leftlim[i] = preleftlim[i];
		rightlim[i] = prerightlim[i];
		fscanf(f, "%lg", &(pr[i]));
		if (i == it0 || i == it02) {
			rightlim[i] += pr[i];
			leftlim[i] += pr[i];
		}
	}
	fclose(f);
	current_path(eventname);
	return 0;
}


void LevMar::ReadAncillary() {
	FILE* f;
	int k;
	double* normfacs;

	try {
		sumf = (double*)malloc(sizeof(double) * nfil);
		sumf2 = (double*)malloc(sizeof(double) * nfil);
		sumfy = (double*)malloc(sizeof(double) * nfil);
		sumy = (double*)malloc(sizeof(double) * nfil);
		sumy2 = (double*)malloc(sizeof(double) * nfil);
		sumsigma = (double*)malloc(sizeof(double) * nfil);

		if (astrometric) {
			sumsigmaN = (double*)malloc(sizeof(double) * nfil);
			sumsigmaE = (double*)malloc(sizeof(double) * nfil);
			sumcN = (double*)malloc(sizeof(double) * nfil);
			sumcE = (double*)malloc(sizeof(double) * nfil);
			sumc1 = (double*)malloc(sizeof(double) * nfil);
			sumc2 = (double*)malloc(sizeof(double) * nfil);
		}

		starts = (int*)malloc(sizeof(int) * nfil);
		sizes = (int*)malloc(sizeof(int) * nfil);
		pr = (double*)realloc(pr, sizeof(double) * (nps + nlinpar * nfil));
		prn = (double*)malloc(sizeof(double) * (nps + nlinpar * nfil));
		errs = (double*)malloc(sizeof(double) * (nps + nlinpar * nfil));
		dFdp = (double*)malloc(sizeof(double) * nlinpar * nfil * nps);
		delta = (double*)malloc(sizeof(double) * nps);
		maxdelta = (double*)malloc(sizeof(double) * nps);
		B = (double*)malloc(sizeof(double) * nps);
		B0 = (double*)malloc(sizeof(double) * (nps));
		A = (double*)malloc(sizeof(double) * nps * nps);
		Curv = (double*)malloc(sizeof(double) * nps * nps);
		Cov = (double*)malloc(sizeof(double) * nps * nps);
		fb = (double*)malloc(sizeof(double) * np * (nps + 1));
		c1l = (double*)malloc(sizeof(double) * np * (nps + 1));
		c2l = (double*)malloc(sizeof(double) * np * (nps + 1));
		c1s = (double*)malloc(sizeof(double) * np * (nps + 1));
		c2s = (double*)malloc(sizeof(double) * np * (nps + 1));
		seps = (double*)malloc(sizeof(double) * np);

		Gr = (double**)malloc(sizeof(double*) * nps);
		for (int i = 0; i < nps; i++) {
			Gr[i] = (double*)malloc(sizeof(double) * np * (nlinpar - 1));
		}

		current_path("Data");

		// If Normalization.txt is present, use numbers therein to normalize datasets.

		normfacs = (double*)malloc(sizeof(double) * nfil);
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
		if (astrometric) sumcN[k] = sumcE[k] = sumsigmaN[k] = sumsigmaE[k] = 0;
		starts[k] = 0;
		for (int i = 0; i < np; i++) {
			w[i] /= normfacs[filter[i]];
			if ((i != 0) && (filter[i] > filter[i - 1])) {
				sizes[k] = i - starts[k];
				k++;
				sumsigma[k] = sumy[k] = sumy2[k] = 0;
				if (astrometric) sumcN[k] = sumcE[k] = sumsigmaN[k] = sumsigmaE[k] = 0;
				starts[k] = i;
			}
			sumsigma[k] += w[i] * w[i];
			sumy[k] += w[i] * w[i] * y[i];
			sumy2[k] += w[i] * w[i] * y[i] * y[i];

			if (astrometric) {
				wcN[i] /= normfacs[filter[i]];
				wcE[i] /= normfacs[filter[i]];
				sumsigmaN[k] += wcN[i] * wcN[i];
				sumsigmaE[k] += wcE[i] * wcE[i];
				sumcN[k] += wcN[i] * wcN[i] * cN[i];
				sumcE[k] += wcE[i] * wcE[i] * cE[i];
			}
		}
		sizes[k] = np - starts[k];

		for (int i = 0; i < nps; i++) {
			maxdelta[i] = sigmapr[i];
		}

		// If LimbDarkening.txt is present, use number therein for linear limb darkening coefficient
		limbdarks = (double*)malloc(sizeof(double) * nfil);
		for (int i = 0; i < nfil; i++) {
			limbdarks[i] = 0;
		}
		if (exists("LimbDarkening.txt")) {
			f = fopen("LimbDarkening.txt", "r");
			for (int i = 0; i < nfil; i++) {
				fscanf(f, "%lf", &limbdarks[i]);
				printf("\nLimbDarkening %d: %lf", i, limbdarks[i]);
			}
			fclose(f);
		}
		current_path("..");

	}
	catch (...) {
		error = 9;
	}
	free(normfacs);
}

int LevMar::Run() {
	FILE* f;
	int il, k, ichi, flag, ilam, bumpnum, bumpcounter;
	double minchi, bestchi, c1, c0, oldlambda, lambda, inclambda, fac, fac2;
	bumper* scanbumper, * scanbumper2;

	if (!error) {
		try {
			// Initializing step chain and bumper list

			printf("\n\n- Initializing step chain \n");

			stepchain = new bumper(pr, nps);
			laststep = stepchain;
			bumperlist = 0;

			// Going to right directory

			if (!exists("PreModels")) create_directory("PreModels");
			current_path("PreModels");

			// Removing all existing files
			sprintf(filename, "%s.*", modelcode);
			auto searchstring = regex(filename);
			for (auto const& itr : directory_iterator(".")) {
				string curfile = (itr).path().filename().string();
				if (regex_match(curfile, searchstring)) {
					remove(curfile);
				}
			}

			il = 0;
			minchi = bestchi = 1.e100;

			time_t ltime, rtime;
			time(&ltime);

			// Calculate nlc models starting from the initial condition
			// il is the number of the current model being calculated
			while (il < nlc) {
				printf("\n*************************************");
				printf("\n************ Curve %d *************", il);
				printf("\n*************************************\n");

				for (int i = 0; i < nps; i++) {
					pr[i] = (laststep->p0)[i];
				}
				bumpcounter = 0;

				PrintOut(pr);
				/* Chi Squared at initial conditions */

				c1 = ChiSquared(pr);
				printf("\nStarting chi2 = %lf", c1);

				// Saving to stepchain
				if (stepchainsave) {
					sprintf(filename, "%s-stepchain%d.dat", modelcode, il);
					PrintFile(filename, il, c1, false);
				}
				// Initializing the main cycle

				lambda = 3.;
				inclambda = 3.;
				bumpnum = 1;
				k = 1;
				ichi = 0;
				// Main fit
				while ((ichi < 3) && (lambda < 1.e10) && (k <= maxsteps)) {
					c0 = c1;
					printf("\nStep %d\n", k++);
					/* Calculation of the gradient */
					//getchar();

					Grad();

					// Preparing the matrices A=(Curv + lambda diag(Curv))
					// with Curv = Gr^T Gr.
					// B=Gr (y-f)*w

					oldlambda = lambda;
					lambda /= inclambda;

					// Levenberg-Marquardt with parameter lambda
					ilam = 0;
					while ((c1 >= c0) && (ilam < 20)) {
						for (int i = 0; i < nps; i++) {
							for (int j = 0; j < nps; j++) {
								A[i * nps + j] = Curv[i * nps + j];
								if (i == j) {
									A[i * nps + j] += lambda * Curv[i * nps + i];
								}
							}
							B[i] = B0[i];
						}

						// Triangularizing the equations A.delta = B

						for (int i = 0; i < nps - 1; i++) {
							for (int j = i + 1; j < nps; j++) {
								fac = -A[j * nps + i] / A[i * nps + i];
								for (int m = i + 1; m < nps; m++) {
									A[j * nps + m] += A[i * nps + m] * fac;
								}
								B[j] += B[i] * fac;
							}
						}

						// Solving for delta

						for (int i = nps - 1; i >= 0; i--) {
							fac = B[i];
							for (int j = i + 1; j < nps; j++) {
								fac -= delta[j] * A[i * nps + j];
							}
							delta[i] = fac / A[i * nps + i];
						}
						// If we end up out of bounds, the point is taken at the bound.
						for (int i = 0; i < nps; i++) {
							//					printf("%lf ",delta[i]);
							if (!((delta[i] > 0) || (delta[i] < 0))) {
								delta[i] = 0.;
							}
							if (fabs(delta[i]) > maxdelta[i]) {
								delta[i] *= maxdelta[i] / fabs(delta[i]);
							}
							prn[i] = pr[i] + delta[i];
							if (prn[i] > rightlim[i]) {
								prn[i] = (0.99 * rightlim[i] + (0.01 + ilam) * pr[i]) / (1 + ilam);
							}
							else if (prn[i] < leftlim[i]) {
								prn[i] = (0.99 * leftlim[i] + (0.01 + ilam) * pr[i]) / (1 + ilam);
							}
						}
						// Current parameters, current chi square and current lambda are displayed
						PrintOut(prn);
						c1 = ChiSquared(prn);
						printf("\nilam= %d lambda= %lf\nc1 = %lf prec=%le", ilam, lambda, c1, Tol);

						lambda *= inclambda;
						ilam++;
					}
					lambda /= inclambda;


					// if new point is better than previous, pr is updated.
					if (ilam < 20) {
						for (int i = 0; i < nps + nlinpar * nfil; i++) {
							pr[i] = prn[i];
						}
					}


					// Bumping mechanism: if point is within the covariance ellipsoid it gets kicked off
					fac = 0;
					flag = 0;
					while (fac < 1 && bumperlist) {
						for (scanbumper = bumperlist; scanbumper; scanbumper = scanbumper->next) {
							fac = scanbumper->distance(pr);
							if (fac < 1) {
								printf("\nBumped!");
								bumpcounter++;
								for (int i = 0; i < nps; i++) {
									//fac = 2.0 * bumperpower / sqrt(fac);
									//prn[i] = pr[i] - fac * scanbumper->dp[i];
									prn[i] = scanbumper->p0[i] - 2.0 * bumperpower * (laststep->p0[i] - scanbumper->p0[i]) / sqrt(fac);
									if (prn[i] > rightlim[i]) {
										prn[i] = 0.99 * rightlim[i] + 0.01 * pr[i];
									}
									else if (prn[i] < leftlim[i]) {
										prn[i] = 0.99 * leftlim[i] + 0.01 * pr[i];
									}
									pr[i] = prn[i];
								}
								scanbumper->UpdateCurvature(bumperpower);
								flag = 1;
								PrintOut(pr);
							}
						}
						if (bumpcounter >= maxbumpcount) break;
					}
					if (flag) {
						c1 = ChiSquared(pr);
						printf("\nc1 = %lf", c1);
					}
					// Add new point to stepchain 
					laststep->next = new bumper(pr, nps);
					laststep = laststep->next;
					laststep->Amp = c1;
					// Saving to file
					if (stepchainsave) {
						sprintf(filename, "%s-stepchain%d.dat", modelcode, il);
						PrintFile(filename, il, c1, false);
					}
					// If new point's chi square is very close to previous one, ichi is increased.
					// When ichi reaches 3, we declare convergence achieved.
					//if (fabs(1 - c1 / c0) < 1.0e-3) {
					if (fabs(c0 - c1) * np < c0) {
						ichi++;
					}
					else {
						ichi = 0;
					}
					// Time limit is set by maxtime
					// Also close if bumpcounter has reached maxbumpcount
					time(&rtime);
					if (difftime(rtime, ltime) > maxtime || bumpcounter >= maxbumpcount) {
						ichi = 10;
						nlc = il;
						printf("\n---- Time limit has been reached. Exiting...");
					}


				}

				if (ichi < 4) {
					// Final chi square of this model
					c0 = ChiSquared(pr);
					printf("\nFinal chi square = %lf\n", c0);

					//Check if this is the best model with this initial condition
					if (c0 < bestchi) bestchi = c0;

					//Grad(); // Are these two lines really needed? They seem a repetition
					//Covariance();

					// Add a new bumper at the position of this model
					scanbumper = bumperlist;
					bumperlist = new bumper(pr, nps);
					bumperlist->Amp = c0;
					bumperlist->next = scanbumper;
					// Calculate covariance matrix
					Grad();
					Covariance();
					bumperlist->SetCurvature(Curv, np / c0);

					// Check if the model makes sense, otherwise set negative chi square for 
					// processing by subsequent programs.
					if ((c0 > 2) && (c0 < 1.e100) && flagblending < np / 2) {
						c0 = c0;
					}
					else {
						c0 = c1 = -1.;
					}

					// Check if the parameters make sense, otherwise set negative chi square for 
					// processing by subsequent programs.
					for (int i = 0; i < nps; i++) {
						if (!((pr[i] > -1.e300) && (pr[i] < 1.e300))) {
							pr[i] = 1.5;
							c0 = c1 = -1.;
						}
					}



					// Store model in file
					if (stepchainsave) {
						sprintf(filename, "%s-%d.txt", modelcode, il);
					}
					else {
						sprintf(filename, "%s.txt", modelcode);
					}
					PrintFile(filename, il, c0, true);

					if (modelcode[0] == 'L' && il == nlc - offsetdegeneracy - 1) {
						// Find best model found so far and fit offset degeneracy model
						scanbumper2 = bumperlist;
						for (scanbumper = bumperlist; scanbumper; scanbumper = scanbumper->next) {
							if (scanbumper->Amp < scanbumper2->Amp) scanbumper2 = scanbumper;
						}
						for (int i = 0; i < nps; i++) {
							pr[i] = scanbumper2->p0[i];
						}

						if (sin(pr[3]) != 0) {
							double s = exp(pr[0]);
							double q = exp(pr[1]);
							q /= 1 + q;
							double xc = pr[2] / sin(pr[3]) + s * q;
							double ac = 2 * s;
							double b = 1 - s * s + ac * xc;
							double snew = (b + sqrt(ac * ac + b * b)) / ac;
							pr[0] = log(snew);
							pr[2] = (xc - snew * q) * sin(pr[3]);
						}
						scanbumper2 = stepchain;
						while (scanbumper2) {
							scanbumper = scanbumper2->next;
							delete scanbumper2;
							scanbumper2 = scanbumper;
						}

						stepchain = laststep = new bumper(pr, nps);
					}
					else {

						// Finding the last step outside any bumper
						k = 1;
						bumper* closestbumper = stepchain, * tentbumper = laststep;
						for (scanbumper = stepchain; scanbumper && (scanbumper->next); scanbumper = scanbumper->next) {
							k++;
							fac = 1.e100;
							for (scanbumper2 = bumperlist; scanbumper2; scanbumper2 = scanbumper2->next) {
								fac2 = scanbumper2->distance(scanbumper->next->p0);
								if (fac2 < fac) {
									fac = fac2;
									tentbumper = scanbumper2;
								}
							}
							printf("\n%d %lf", k, fac);
							if (fac < 1.) {
								laststep = scanbumper;
								k--;
								scanbumper2 = scanbumper->next;
								while (scanbumper2) {
									scanbumper = scanbumper2->next;
									delete scanbumper2;
									scanbumper2 = scanbumper;
								}
								laststep->next = 0;
								scanbumper = laststep;
							}
							else closestbumper = tentbumper;
						}
						//						if (closestbumper != stepchain) {
													// Start next chain on the other side of closest bumper
						for (int i = 0; i < nps; i++) {
							prn[i] = closestbumper->p0[i] - bumperpower * (laststep->p0[i] - closestbumper->p0[i]);
							if (prn[i] > rightlim[i]) {
								prn[i] = 0.99 * rightlim[i] + 0.01 * pr[i];
							}
							else if (prn[i] < leftlim[i]) {
								prn[i] = 0.99 * leftlim[i] + 0.01 * pr[i];
							}
							pr[i] = prn[i];

						}
						laststep->next = new bumper(pr, nps);
						laststep = laststep->next;
						//						}
					}
					// Updating time count
					//printf("\npartial time=%lf secs\n",(Environment::TickCount-tm)/1000.0);
					//printf("total time = %lf mins\n",(Environment::TickCount-tim0)/60000.0);
					//tm=Environment::TickCount;

					// Updating the nlc.dat file
					//f = fopen("nlc.dat", "w");
					//fprintf(f, "%d %d", il, nlc);
					//fclose(f);
				}

				il++;
			}

			//// When terminating, update the Termination file, 
			//// so that subsequent programs know everything has been concluded with this fit
			//sprintf(filename, "t%s.dat", outdir);
			//f = fopen(filename, "w");
			//fprintf(f, "%d %d", il, nlc);
			//fclose(f);

			// If this fit has found the best minimum so far, update minchi.dat
//			current_path("..");
			strcpy(filename, "minchi.dat");
			if (exists(filename)) {
				f = fopen(filename, "r");
				fscanf(f, "%lf", &minchi);
				fclose(f);
			}
			if (bestchi < minchi) {
				f = fopen(filename, "w");
				fprintf(f, "%lf %s", bestchi, outdir);
				fclose(f);
			}

		}
		catch (...) {
			error = 8;
		}
	}
	return error;
}

void LevMar::EvaluateModel(double* pr, int fl, int ips) {
	double* tfl, * fbfl, * c1sfl, * c2sfl, * c1lfl, * c2lfl;
	tfl = &(t[starts[fl]]);
	fbfl = &(fb[starts[fl] + np * ips]);
	c1sfl = &(c1s[starts[fl] + np * ips]);
	c2sfl = &(c2s[starts[fl] + np * ips]);
	c1lfl = &(c1l[starts[fl] + np * ips]);
	c2lfl = &(c2l[starts[fl] + np * ips]);
	switch (modnumber) {
	case 0:
		VBM->ESPLLightCurve(pr, tfl, fbfl, y1a, y2a, sizes[fl]);
		break;
	case 1:
		if (astrometric) {
			VBM->ESPLAstroLightCurve(pr, tfl, fbfl, c1sfl, c2sfl, c1lfl, c2lfl, y1a, y2a, sizes[fl]);
		}
		else {
			VBM->ESPLLightCurveParallax(pr, tfl, fbfl, y1a, y2a, sizes[fl]);
		}
		break;
	case 2:
		VBM->BinSourceExtLightCurve(pr, tfl, fbfl, y1a, y2a, sizes[fl]);
		break;
	case 3:
		if (astrometric) {
			VBM->BinSourceAstroLightCurveXallarap(pr, tfl, fbfl, c1sfl, c2sfl, c1lfl, c2lfl, y1a, y2a, y1a, y2a, sizes[fl]);
		}
		else {
			VBM->BinSourceExtLightCurveXallarap(pr, tfl, fbfl, y1a, y2a, y1a, y2a, sizes[fl]);
		}
		break;
	case 4:
		VBM->BinaryLightCurve(pr, tfl, fbfl, y1a, y2a, sizes[fl]);
		break;
	case 5:
		if (astrometric) {
			VBM->BinaryAstroLightCurve(pr, tfl, fbfl, c1sfl, c2sfl, c1lfl, c2lfl, y1a, y2a, sizes[fl]);
		}
		else {
			VBM->BinaryLightCurveParallax(pr, tfl, fbfl, y1a, y2a, sizes[fl]);
		}
		break;
	case 6:
		if (astrometric) {
			VBM->BinaryAstroLightCurveOrbital(pr, tfl, fbfl, c1sfl, c2sfl, c1lfl, c2lfl, y1a, y2a, seps, sizes[fl]);
		}
		else {
			VBM->BinaryLightCurveOrbital(pr, tfl, fbfl, y1a, y2a, seps, sizes[fl]);
		}
		break;
	case 7:
		if (astrometric) {
			VBM->BinaryAstroLightCurveKepler(pr, tfl, fbfl, c1sfl, c2sfl, c1lfl, c2lfl, y1a, y2a, seps, sizes[fl]);
		}
		else {
			VBM->BinaryLightCurveKepler(pr, tfl, fbfl, y1a, y2a, seps, sizes[fl]);
		}
		break;
	case 8:
		VBM->TripleLightCurve(pr, tfl, fbfl, y1a, y2a, sizes[fl]);
		break;
	case 9:
		if (astrometric) {
			VBM->TripleAstroLightCurve(pr, tfl, fbfl, c1sfl, c2sfl, c1lfl, c2lfl, y1a, y2a, sizes[fl]);
		}
		else {
			VBM->TripleLightCurveParallax(pr, tfl, fbfl, y1a, y2a, sizes[fl]);
		}
		break;
	}
}

double LevMar::ChiSquared(double* pr) {
	double chi2 = 0, chi0, chia, p1;
	double p1max = 0, maxsum = 0;
	double tmax = 0;

	maxmaxsum = 0;

	for (int fl = 0; fl < nfil; fl++) {
		VBM->satellite = satel[starts[fl]];
		VBM->a1 = limbdarks[fl];
		EvaluateModel(pr, fl, 0); // Populates fb[starts[fl]] with magnifications for given parameters
		// For astrometric datasets, populates also source and lens centroids
		sumf[fl] = sumfy[fl] = sumf2[fl] = 0;
		for (int i = starts[fl]; i < starts[fl] + sizes[fl]; i++) {
			sumf[fl] += w[i] * w[i] * fb[i];
			sumf2[fl] += w[i] * w[i] * fb[i] * fb[i];
			sumfy[fl] += w[i] * w[i] * fb[i] * y[i];
		}
	}
	for (int i = 0; i < nfil; i++) {
		if (w[starts[i]] > 0) {
			p1 = sumf[i] * sumf[i] - sumf2[i] * sumsigma[i] + epsilon;
			pr[nps + i * nlinpar] = (sumf[i] * sumfy[i] - sumf2[i] * sumy[i]) / p1;
			pr[nps + 1 + i * nlinpar] = (sumf[i] * sumy[i] - sumsigma[i] * sumfy[i]) / p1;
			//chi2+=fabs(sumy2[i]+(sumfy[i]*sumfy[i]*sumsigma[i]+sumf2[i]*sumy[i]*sumy[i]-2*sumf[i]*sumy[i]*sumfy[i])/p1);

			if (pr[nps + 1 + i * nlinpar] < 0) {
				pr[nps + 1 + i * nlinpar] = 0;
			}
		}
		else {
			if (w[starts[i]] > -0.75) { // Lens astrometry
				pr[nps + i * nlinpar] = 1;
				pr[nps + 1 + i * nlinpar] = 0;
			}
			else {					// Source astrometry
				pr[nps + i * nlinpar] = 0;
				pr[nps + 1 + i * nlinpar] = 1;
			}
		}

	}
	chi2 = chi0 = chia = 0;

	if (astrometric) { // Astrometric chi square
		double g;
		for (int fl = 0; fl < nfil; fl++) {
			if (wcN[starts[fl]] > 1.e-90) { // Only calculate for datasets with astrometric data
				g = pr[nps + fl * nlinpar] / (pr[nps + fl * nlinpar + 1] + pr[nps + fl * nlinpar] * 1.e-8); // Blending for this dataset
				sumc1[fl] = sumc2[fl] = 0;
				for (int i = starts[fl]; i < starts[fl] + sizes[fl]; i++) {
					c1s[i] = (c1s[i] * fb[i] + c1l[i] * g) / (fb[i] + g); // Weighted centroid
					c2s[i] = (c2s[i] * fb[i] + c2l[i] * g) / (fb[i] + g);
					sumc1[fl] += wcN[i] * wcN[i] * c1s[i];
					sumc2[fl] += wcE[i] * wcE[i] * c2s[i];
				}
				pr[nps + fl * nlinpar + 2] = (sumcN[fl] - sumc1[fl]) / sumsigmaN[fl]; // Origin shift
				pr[nps + fl * nlinpar + 3] = (sumcE[fl] - sumc2[fl]) / sumsigmaE[fl];
				for (int i = starts[fl]; i < starts[fl] + sizes[fl]; i++) {
					p1 = (cN[i] - pr[nps + fl * nlinpar + 2] - c1s[i]) * wcN[i];
					chia += p1 * p1;
					p1 = (cE[i] - pr[nps + fl * nlinpar + 3] - c2s[i]) * wcE[i];
					chia += p1 * p1;

				}
			}
			else {
				pr[nps + fl * nlinpar + 2] = pr[nps + fl * nlinpar + 3] = 0;
			}
		}
	}

	// Photometric chi square
	flagblending = 0;
	for (int i = 0; i < np; i++) {
		if (w[i] > 0) {
			p1 = (y[i] - pr[nps + filter[i] * nlinpar] - pr[nps + 1 + filter[i] * nlinpar] * fb[i]) * w[i] * w[i] * pr[nps + 1 + filter[i] * nlinpar] * Tol;
			chi0 += p1 * p1;
			p1 = (y[i] - pr[nps + filter[i] * nlinpar] - pr[nps + 1 + filter[i] * nlinpar] * fb[i]) * w[i];

			if (p1 > 0) {
				maxsum += p1; // somma i residui 

				if (p1 > p1max) {
					p1max = p1; // massimo residuo positivo
					//std::cout << "\nIl massimo residuo positivo è: " << p1max << " al punto " << i << std::endl;
					tmax = t[i]; // tempo in cui si ha il massimo residuo positivo
				}

			}
			else {
				//printf("\nResiduo negativo o nullo %d: %lf", i, p1); // stampa il residuo negativo
				// Reset in caso di interruzione della sequenza positiva
				maxsum = 0;                  // azzera la somma dei residui

			}
			//Alla fine del ciclo, se la somma dei residui positivi consecutivi è maggiore della somma massima trovata finora allora aggiorna la somma massima e il tempo corrispondente
			if (maxsum > maxmaxsum) {
				maxmaxsum = maxsum; // somma massima dei residui positivi consecutivi	
				//std::cout << "\nLa somma massima dei residui positivi consecutivi è: " << maxmaxsum << " al tempo " << tmax << std::endl;
				tmaxmax = tmax; // tempo in cui si ha la somma massima dei residui positivi consecutivi
			}
			//std::cout << std::endl;
		}
		//else {
			//printf("\nNessun residuo positivo trovato.");
		//}


		//robustezza del fit 
		//bisogna considerare il residuo di una seguenza di punti consecutivi che hanno lo stesso segno
		//**** il segno è dato dal confronto tra il modello ed il flusso misurato
		//**** se il modello è minore del flusso misurato allora il flusso ha un picco che il modello non ha
		//**** viceversa si ha un deep 
		//quindi bisogna considerare la somma dei residui dei punti consecutivi che hanno lo stesso segno
		//in questo modo bisogna tener conto del segno di p1 ad ogni passo e confrontarlo con il segno del residuo precedente
		//se ha lo stesso segno allora si sommano i residui
		// in caso di più sequenze si considera qualla in cui la somma dei residui è max
		//infine si calcola la massima deviazione standard dei residui all'interno di una sequenza di punti consecutivi che hanno lo stesso segno

		chi2 += p1 * p1;
		if (pr[nps + 1 + filter[i] * nlinpar] > 2 * y[i]) {
			flagblending++;
			chi2 += (pr[nps + 1 + filter[i] * nlinpar] - 2 * y[i]) * (pr[nps + 1 + filter[i] * nlinpar] - 2 * y[i]) * w[i] * w[i];
		}
	}
	if (maxsum > maxmaxsum) {
		maxmaxsum = maxsum; // somma massima dei residui positivi consecutivi	
		//std::cout << "\nLa somma massima dei residui positivi consecutivi è: " << maxmaxsum << " al tempo " << tmax << std::endl;
		tmaxmax = tmax; // tempo in cui si ha la somma massima dei residui positivi consecutivi
	}
	chi0 = sqrt(2 * chi0); // Error in chi square
	if (chi0 / chi2 > 0.1) Tol *= 0.5;
	if (chi0 / chi2 < 0.01 && Tol < .99e-2) Tol *= 2;
	chi2 += chia; // Add astrometric chi square

	// Constraints
	for (int icons = 0; icons < consnumber; icons++) {
		consvars[icons] = ComputeConstraint(pr, icons);
		p1 = consvars[icons] - constraints[icons];
		p1 /= (p1 > 0) ? consright[icons] : consleft[icons];
		chi2 += p1 * p1;
	}

	return chi2;
}

void LevMar::Grad() {
	static double inc = 1.0e-3, p1;
	for (int i = 0; i < nps + nfil * nlinpar; i++) {
		prn[i] = pr[i];
	}
	for (int j = 0; j < nps; j++) {
		printf("%d ", j);
		prn[j] += inc;
		for (int fl = 0; fl < nfil; fl++) {
			VBM->satellite = satel[starts[fl]];
			VBM->a1 = limbdarks[fl];
			EvaluateModel(prn, fl, j + 1); // Populates fb[starts[fl]] with magnifications for given parameters
			// For astrometric datasets, populates also source and lens centroids
			sumf[fl] = sumfy[fl] = sumf2[fl] = 0;
			for (int i = starts[fl]; i < starts[fl] + sizes[fl]; i++) {
				sumf[fl] += w[i] * w[i] * fb[i + np * (j + 1)];
				sumf2[fl] += w[i] * w[i] * fb[i + np * (j + 1)] * fb[i + np * (j + 1)];
				sumfy[fl] += w[i] * w[i] * fb[i + np * (j + 1)] * y[i];
			}
		}
		for (int i = 0; i < nfil; i++) {
			if (w[starts[i]] > 0) {
				p1 = sumf[i] * sumf[i] - sumf2[i] * sumsigma[i] + epsilon;
				prn[nps + i * nlinpar] = (sumf[i] * sumfy[i] - sumf2[i] * sumy[i]) / p1;
				prn[nps + 1 + i * nlinpar] = (sumf[i] * sumy[i] - sumsigma[i] * sumfy[i]) / p1;

				if (prn[nps + 1 + i * nlinpar] < 0)  prn[nps + 1 + i * nlinpar] = 0;
			}
			else {
				if (w[starts[i]] > -0.75) { // Lens astrometry
					prn[nps + i * nlinpar] = 1;
					prn[nps + 1 + i * nlinpar] = 0;
				}
				else {					// Source astrometry
					prn[nps + i * nlinpar] = 0;
					prn[nps + 1 + i * nlinpar] = 1;
				}
			}
			dFdp[(1 + nlinpar * i) * nps + j] = -1.08574 * ((prn[nps + i * nlinpar] + prn[nps + 1 + i * nlinpar]) / (pr[nps + i * nlinpar] + pr[nps + 1 + i * nlinpar]) - 1) / inc;  // error on baseline -2.5log(FB+FS)/log(10) expanded to first order
			dFdp[(nlinpar * i) * nps + j] = (prn[nps + i * nlinpar] / prn[nps + 1 + i * nlinpar] - pr[nps + i * nlinpar] / pr[nps + 1 + i * nlinpar]) / inc;    // error on blending FB/FS
		}

		if (astrometric) { // Astrometric calculations
			double g;
			for (int fl = 0; fl < nfil; fl++) {
				if (wcN[starts[fl]] > 1.e-90) { // Only calculate for datasets with astrometric data
					g = prn[nps + fl * nlinpar] / (prn[nps + fl * nlinpar + 1] + 1.e-12 * prn[nps + fl * nlinpar]); // Blending for this dataset
					sumc1[fl] = sumc2[fl] = 0;
					for (int i = starts[fl]; i < starts[fl] + sizes[fl]; i++) {
						c1s[i + np * (j + 1)] = (c1s[i + np * (j + 1)] * fb[i + np * (j + 1)] + c1l[i + np * (j + 1)] * g) / (fb[i + np * (j + 1)] + g); // Weighted centroid
						c2s[i + np * (j + 1)] = (c2s[i + np * (j + 1)] * fb[i + np * (j + 1)] + c2l[i + np * (j + 1)] * g) / (fb[i + np * (j + 1)] + g);
						sumc1[fl] += wcN[i] * wcN[i] * c1s[i + np * (j + 1)];
						sumc2[fl] += wcE[i] * wcE[i] * c2s[i + np * (j + 1)];
					}
					prn[nps + fl * nlinpar + 2] = (sumcN[fl] - sumc1[fl]) / sumsigmaN[fl]; // Origin shift
					prn[nps + fl * nlinpar + 3] = (sumcE[fl] - sumc2[fl]) / sumsigmaE[fl];
					dFdp[(2 + nlinpar * fl) * nps + j] = (prn[nps + fl * nlinpar + 2] - pr[nps + fl * nlinpar + 2]) / inc;
					dFdp[(3 + nlinpar * fl) * nps + j] = (prn[nps + fl * nlinpar + 3] - pr[nps + fl * nlinpar + 3]) / inc;
				}
				else {
					prn[nps + fl * nlinpar + 2] = prn[nps + fl * nlinpar + 3] = 0;
				}
			}
		}

		for (int icons = 0; icons < consnumber; icons++) { // Gradient of constraints
			consvars[icons + (j + 1) * consnumber] = ComputeConstraint(prn, icons);
			p1 = (consvars[icons] - constraints[icons]);
			consvars[icons + (j + 1) * consnumber] = (consvars[icons + (j + 1) * consnumber] - consvars[icons]) / (((p1 > 0) ? consright[icons] : consleft[icons]) * inc);
		}
		prn[j] -= inc;
		//		double errgrad = 0,errterm, grad=0;
		for (int i = 0; i < np; i++) {
			if (w[i] > 0) {
				Gr[j][i] = w[i] * (prn[nps + filter[i] * nlinpar] + prn[nps + 1 + filter[i] * nlinpar] * fb[i + np * (j + 1)] - pr[nps + filter[i] * nlinpar] - pr[nps + 1 + filter[i] * nlinpar] * fb[i]) / inc;
				//errterm = w[i] / inc * Tol * pr[nps + 1 + filter[i] * nlinpar];
				//errgrad += errterm * errterm;
				//grad += Gr[j][i] * Gr[j][i];
			}
			if (wcN[i] > 0) {
				Gr[j][i + np] = wcN[i] * (prn[nps + filter[i] * nlinpar + 2] + c1s[i + np * (j + 1)] - pr[nps + filter[i] * nlinpar + 2] - c1s[i]) / inc;
				Gr[j][i + np * 2] = wcE[i] * (prn[nps + filter[i] * nlinpar + 3] + c2s[i + np * (j + 1)] - pr[nps + filter[i] * nlinpar + 3] - c2s[i]) / inc;
				//errterm= (50* exp(pr[3])*Tol*pr[9]) / inc;
				//errgrad += errterm * errterm * (wcN[i] * wcN[i] + wcE[i] * wcE[i]);
				//grad += Gr[j][i+np] * Gr[j][i+np] + Gr[j][i+2*np] * Gr[j][i + 2 * np];
			}
		}
		//printf("%d ", j);
	}
	printf("OK\n");
	// Curvature matrix
	for (int j = 0; j < nps; j++) {
		for (int i = 0; i < nps; i++) {
			Curv[i * nps + j] = 0;
			for (int k = 0; k < np; k++) {
				if (w[k] > 0) {
					Curv[i * nps + j] += Gr[i][k] * Gr[j][k];
				}
				if (wcN[k] > 0) {
					Curv[i * nps + j] += Gr[i][k + np] * Gr[j][k + np];
					Curv[i * nps + j] += Gr[i][k + np * 2] * Gr[j][k + np * 2];
				}
			}
			// Constraints in curvature
			for (int icons = 0; icons < consnumber; icons++) {
				Curv[i * nps + j] += consvars[icons + (i + 1) * consnumber] * consvars[icons + (j + 1) * consnumber];
			}
		}
	}

	// Offset
	for (int i = 0; i < nps; i++) {
		p1 = 0;
		for (int k = 0; k < np; k++) {
			if (w[k] > 0) {
				p1 += w[k] * Gr[i][k] * (y[k] - pr[nps + filter[k] * nlinpar] - pr[nps + 1 + filter[k] * nlinpar] * fb[k]);
			}
			if (wcN[k] > 0) {
				p1 += wcN[k] * Gr[i][k + np] * (cN[k] - pr[nps + filter[k] * nlinpar + 2] - c1s[k]);
				p1 += wcE[k] * Gr[i][k + np * 2] * (cE[k] - pr[nps + filter[k] * nlinpar + 3] - c2s[k]);
			}
		}
		B0[i] = p1;
		for (int icons = 0; icons < consnumber; icons++) {
			p1 = (consvars[icons] - constraints[icons]);
			p1 /= (p1 > 0) ? consright[icons] : consleft[icons];
			B0[i] -= p1 * consvars[icons + (i + 1) * consnumber];
		}
	}
}

inline double LevMar::ComputeConstraint(double* pr, int ic) {
	int i = consindex[ic];
	if (i < 100) {
		return pr[i];
	}
	if (i < 10000) {
		return pr[(i - 100) + nps - 4];
	}
	if (i < 20000) {
		return pr[nps + (i - 10000) * 2] / pr[nps + (i - 10000) * 2 + 1];
	}
	if (i == 30000) {
		int posN = -1, posE = -1;
		for (int i = 0; i < nps; i++) {
			if (strcmp(parnames[modnumber][i].c_str(), "piN") == 0) posN = i;
			if (strcmp(parnames[modnumber][i].c_str(), "piE") == 0) posE = i;
		}
		return atan2(pr[posE], pr[posN]);
	}
	if (i == 30001) {
		int postE = -1, posrho = -1;
		for (int i = 0; i < nps; i++) {
			if (strcmp(parnames[modnumber][i].c_str(), "tE") == 0) postE = i;
			if (strcmp(parnames[modnumber][i].c_str(), "rho") == 0) posrho = i;
		}
		return exp(pr[postE] + pr[posrho]);
	}

	return 0;
}

void LevMar::Covariance() {
	Inverse(Curv, Cov, nps);
	for (int i = 0; i < nps; i++) {
		errs[i] = sqrt(fabs(Cov[i * nps + i]));
	}
	for (int i = 0; i < nlinpar * nfil; i++) {
		errs[i + nps] = 0;
		for (int j1 = 0; j1 < nps; j1++) {
			for (int j2 = 0; j2 < nps; j2++) {
				errs[i + nps] += Cov[j1 * nps + j2] * dFdp[i * nps + j1] * dFdp[i * nps + j2];
			}
		}
		errs[i + nps] = sqrt(fabs(errs[i + nps]));
		if (!(errs[i + nps] > 0)) errs[i + nps] = 1.e100;
	}
}



void LevMar::PrintOut(double* pr) {
	int npp = nps, ilog = 0, logsize = logposs[modnumber].size();
	double fl;
	if (astrometric) npp -= 4;

	for (int i = 0; i < npp; i++) {
		if (i % 4 == 0) {
			printf("\n");
		}
		if (ilog < logsize && i == logposs[modnumber][ilog]) {
			fl = exp(pr[i]);
			ilog++;
		}
		else {
			fl = pr[i];
		}
		printf("%s=%lf ", parnames[modnumber][i].c_str(), fl);
	}
	if (astrometric) {
		printf("\n");
		for (int i = 0; i < 4; i++) {
			printf("%s=%lf ", astroparnames[i].c_str(), pr[npp + i]);
		}
	}
}

void LevMar::PrintFile(char* filename, int il, double c0, bool printerrors) {
	int npp = nps, ilog, logsize = logposs[modnumber].size();
	double fl;
	FILE* f;
	f = fopen(filename, "a");
	switch (modnumber) {
		//case 2:
		//	if (pr[1] > 0 && !VBM->turn_off_secondary_source) { // Invert sources
		//		double sc;
		//		pr[1] = -pr[1];
		//		sc = pr[3];
		//		pr[3] = pr[2];
		//		pr[2] = sc;
		//		sc = pr[5];
		//		pr[5] = pr[4];
		//		pr[4] = sc;
		//		sc = errs[3];
		//		errs[3] = errs[2];
		//		errs[2] = sc;
		//		sc = errs[5];
		//		errs[5] = errs[4];
		//		errs[4] = sc;
		//		for (int k = 0; k < nps; k++) {
		//			Cov[1 + nps * k] = -Cov[1 + nps * k];
		//			Cov[k + nps * 1] = -Cov[k + nps * 1];
		//		}
		//		for (int k = 0; k < nps; k++) {
		//			sc = Cov[2 + nps * k];
		//			Cov[2 + nps * k] = Cov[3 + nps * k];
		//			Cov[3 + nps * k] = sc;
		//			sc = Cov[5 + nps * k];
		//			Cov[5 + nps * k] = Cov[4 + nps * k];
		//			Cov[4 + nps * k] = sc;
		//		}
		//		for (int k = 0; k < nps; k++) {
		//			sc = Cov[k + nps * 2];
		//			Cov[k + nps * 2] = Cov[k + nps * 3];
		//			Cov[k + nps * 3] = sc;
		//			sc = Cov[k + nps * 5];
		//			Cov[k + nps * 5] = Cov[k + nps * 4];
		//			Cov[k + nps * 4] = sc;
		//		}
		//	}
		//	break;
	case 4:
		if (pr[1] > 0 && !VBM->turn_off_secondary_lens) {
			pr[3] = pr[3] - M_PI;
			pr[1] = -pr[1];
			for (int k = 0; k < nps; k++) {
				Cov[1 + nps * k] = -Cov[1 + nps * k];
				Cov[k + nps * 1] = -Cov[k + nps * 1];
			}
		}
		while (pr[3] > 2 * M_PI) pr[3] -= 2 * M_PI;
		while (pr[3] < 0) pr[3] += 2 * M_PI;
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
		break;
	case 5:
		if (pr[1] > 0 && !VBM->turn_off_secondary_lens) {
			pr[3] = pr[3] - M_PI;
			pr[1] = -pr[1];
			for (int k = 0; k < nps; k++) {
				Cov[1 + nps * k] = -Cov[1 + nps * k];
				Cov[k + nps * 1] = -Cov[k + nps * 1];
			}
		}
		while (pr[3] > 2 * M_PI) pr[3] -= 2 * M_PI;
		while (pr[3] < 0) pr[3] += 2 * M_PI;
		break;
	case 6:
		if (pr[1] > 0 && !VBM->turn_off_secondary_lens) {
			pr[3] = pr[3] - M_PI;
			pr[1] = -pr[1];
			for (int k = 0; k < nps; k++) {
				Cov[1 + nps * k] = -Cov[1 + nps * k];
				Cov[k + nps * 1] = -Cov[k + nps * 1];
			}
		}
		while (pr[3] > 2 * M_PI) pr[3] -= 2 * M_PI;
		while (pr[3] < 0) pr[3] += 2 * M_PI;
		break;
	case 8:
		while (pr[3] > 2 * M_PI) pr[3] -= 2 * M_PI;
		while (pr[3] < 0) pr[3] += 2 * M_PI;
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

		while (pr[9] > 2 * M_PI) pr[9] -= 2 * M_PI;
		while (pr[9] < 0) pr[9] += 2 * M_PI;
		break;
	case 9:
		while (pr[3] > 2 * M_PI) pr[3] -= 2 * M_PI;
		while (pr[3] < 0) pr[3] += 2 * M_PI;

		while (pr[9] > 2 * M_PI) pr[9] -= 2 * M_PI;
		while (pr[9] < 0) pr[9] += 2 * M_PI;
		break;
	}

	if (astrometric) npp -= 4;

	if (printerrors) {
		//Write model name
		filename[strlen(filename) - 4] = 0;
		fprintf(f, "%s-%d\n", filename, il);
	}

	// Write values of the parameters (without log)
	ilog = 0;
	for (int i = 0; i < npp; i++) {
		if (ilog < logsize && i == logposs[modnumber][ilog]) {
			fl = exp(pr[i]);
			ilog++;
		}
		else {
			fl = pr[i];
		}
		fprintf(f, "%.16le ", fl);
	}
	// Write astrometric parameters
	if (astrometric) {
		for (int i = 0; i < 4; i++) {
			fprintf(f, "%.16le ", pr[npp + i]);
		}
	}
	for (int i = nps; i < nps + nlinpar * nfil; i++) {
		pr[i] = ((pr[i] > -1.e300) && (pr[i] < 1.e300)) ? pr[i] : -1.e300;
		fprintf(f, "%le ", pr[i]);
	}
	//stampa tmaxmax 
	fprintf(f, "%.16le ", tmaxmax);

	//stampa maxmaxsum
	fprintf(f, "%.16le ", maxmaxsum);
	// Write chi square
	fprintf(f, "%.16le\n", c0);

	if (printerrors) {
		// Write errors
		ilog = 0;
		for (int i = 0; i < npp; i++) {
			if (ilog < logsize && i == logposs[modnumber][ilog]) {
				fl = exp(pr[i]) * errs[i];
				ilog++;
			}
			else {
				fl = errs[i];
			}
			fprintf(f, "%.16le ", fl);
		}
		// Write astrometric errors
		if (astrometric) {
			for (int i = 0; i < 4; i++) {
				fprintf(f, "%.16le ", errs[npp + i]);
			}
		}
		// Write errors for fluxes
		fprintf(f, "%le", errs[nps]);
		for (int i = nps + 1; i < nps + nlinpar * nfil; i++) {
			fprintf(f, " %le", errs[i]);
		}
		fprintf(f, "\n");

		// Print covariance matrix to file
		for (int i = 0; i < nps; i++) {
			if (i > 0) fprintf(f, "\n");
			fprintf(f, "%le", Cov[i]);
			for (int j = 1; j < nps; j++) {
				fprintf(f, " %le", Cov[i + nps * j]);
			}
		}
		fprintf(f, "\n");
	}
	fclose(f);
}