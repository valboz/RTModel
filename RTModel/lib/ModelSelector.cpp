// ModelSelector.cpp : main project file.
// This program select the best model and viable alternatives from PreModels of a given class

#define _USE_MATH_DEFINES
#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <malloc.h>
#include <filesystem>
#include <regex>
#include "bumper.h"

using namespace std;
using namespace std::filesystem;

double supfac = 3.0; // factor multiplying the inverse covariance in search for superpositions (models are incompatible if farther than supfac*sigma)
double chifac = 1.; // number of sigmas in the chi square distribution for accepting alternative models after the best one
int maxmodels = 10; // maximum number of models returned

int main(int argc, char* argv[]) {

	char eventname[256] = "";
	char filename[256] = "", filnum[256] = "", filnum2[10] = "";
	char modelcode[256];
	char command[256], buffer[256];
	double value;
	double t,  y,  w, * pr, * sigmapr, * Cov, *Curv, fac, facr, c1, c0, chithr, bestplan = 1.e100, bestbin = 1.e100;
	double renorm;
	int nfil, il, nlc, nmod, np, k;
	int nps=4,dof;
	FILE* f, * g;

	bumper* bumperlist=0, *scanbumper, *scanbumper2;


	printf("******************************************\n");
	printf("**********     Model selector    *********\n");
	printf("******************************************\n\n\n");
	printf("This program finds the best independent models of a given class\n\n");

	// Directory preliminaries. Reads event name from arguments.

	auto exedir = current_path();

	if (argc > 2) {
		strcpy(eventname, argv[1]);
		strcpy(modelcode, argv[2]);
	}
	else {
		printf("\n\nEvent name? ");
		scanf("%s", eventname);
		printf("\n\nModel code? ");
		scanf("%s", modelcode);
	}

	printf("\n\n- Event: %s\n", eventname);

	current_path(eventname);

	if (exists("Models")) {
		current_path("Models");
		auto searchstring = regex(string(modelcode) + ".*");
		for (auto const& itr : directory_iterator(".")) {
			string curfile = (itr).path().filename().string();
			if (regex_match(curfile, searchstring)) {
				printf("\n\nModels already selected");
				return 0;
			}
			itr;
		}
	}
	current_path(eventname);

	if (exists("ini")) {
		current_path("ini");
		f = fopen("ModelSelector.ini", "r");
		if (f != 0) {
			printf("\n\n- Reading options in ModelSelector.ini");
			while (!feof(f)) {
				int red = fscanf(f, "%s %s %lf", command, buffer, &value);
				if (red < 1) {
					command[0] = 0;
					//if (red != 0) {
					//	printf("\n\n!!! Bad command in Reader.ini");
					//	return -1;
					//};
				}
				if (strcmp(command, "maxmodels") == 0) {
					maxmodels = (int) value;
				}
				if (strcmp(command, "sigmasoverlap") == 0) {
					supfac = value;
				}
				if (strcmp(command, "sigmaschisquare") == 0) {
					chifac = value;
				}
			}
			fclose(f);
		}
		else {
			printf("\n\n- Default options:");
		}
	}
	printf("\nNumber of sigmas used to declare overlap between two models: %lf", supfac);
	printf("\nNumber of sigmas in chi square distribution for model acceptance: %lf", chifac);
	printf("\nMaximum number of models reported: %d", maxmodels);
	supfac *= supfac;

	switch (modelcode[0]) {
	case 'P':
		nps = 4;
		break;
	case 'B':
		nps = 7;
		if (modelcode[1] == 'O') {
			nps += 3;
		}
		break;
	case 'L':
		nps = 7;
		if (modelcode[1] == 'O') {
			nps += 5;
		}
		break;
	default:
		printf("\n\n - Invalid model code!!!");
		return -1;
	}
	if (modelcode[1] == 'X') {
		nps += 2;
	}

	supfac *= nps;

	printf("\n- Model code: %s", modelcode);

	// Read curve to fit

	printf("\n\nReading data\n");

	
	current_path(eventname);
	
	f = fopen("LCToFit.txt", "r");
	fscanf(f, "%d", &np);
	nfil = 1;
	for (int i = 0; i < np; i++) {
		fscanf(f, "%d %lf %lf %lf %d", &il, &t, &y, &w, &k);
		if ((i != 0) && il!=nfil-1) {
			nfil++;
		}
	}
	fclose(f);

	pr = (double*)malloc(sizeof(double) * (nps + 2 * nfil));
	sigmapr = (double*)malloc(sizeof(double) * (nps + 2 * nfil));
	Cov = (double*)malloc(sizeof(double) * (nps * nps));
	Curv = (double*)malloc(sizeof(double) * (nps * nps));

	dof = np - nps; // degrees of freedom

	// Read models

	printf("\n\nReading models\n");

	if (!exists("PreModels")) {
		printf("\n! No models for this class");
		return 0;
	}
	current_path("PreModels");
	nmod = 0;

	auto searchstring = regex(string(modelcode) + "[0-9]*");
	for (auto const& itr : directory_iterator(".")) {
		if (regex_match((itr).path().filename().string(), searchstring) && (itr).is_directory()) {
			current_path(itr);

			f = fopen("nlc.dat", "r");
			fscanf(f, "%d %d", &il, &nlc);
			fclose(f);

			for (int i1 = 0; i1 <= il; i1++) {
				sprintf(filename, "%d.txt", i1);
				//					if(_findfirst(filename,&strfile)==-1) continue;
				f = fopen(filename, "r");
				for (int j = 0; j < nps + 2 * nfil; j++) {
					fscanf(f, "%le", &(pr[j]));
				}
				fscanf(f, "%le", &(c0));

				renorm = sqrt(c0 / dof);

				for (int j = 0; j < nps + 2 * nfil; j++) {
					fscanf(f, "%le", &(sigmapr[j]));
					if (!(sigmapr[j] >= 0)) c0 = -1;
				}

				for (int i = 0; i < nps; i++) {
					for (int j = 0; j < nps; j++) {
						fscanf(f, "%le", &(Cov[i + j * nps]));
					}
				}

				switch (modelcode[0]) {
				case 'P':
					sigmapr[1] = sigmapr[1] / pr[1];
					sigmapr[3] = sigmapr[3] / pr[3];
					pr[1] = log(pr[1]);
					pr[3] = log(pr[3]);
					if (modelcode[1] == 'S') {
						sigmapr[0] = sigmapr[0] / pr[0];
						pr[0] = log(pr[0]);
					}
					break;
				case 'B':
					if (modelcode[1] == 'O') {
						sigmapr[2] = sigmapr[2] / pr[2];							
						sigmapr[3] = sigmapr[3] / pr[3];
						sigmapr[9] = sigmapr[9] / pr[9];
						pr[2] = log(pr[2]);
						pr[3] = log(pr[3]);
						pr[9] = log(pr[9]);
					}
					else {
						sigmapr[0] = sigmapr[0] / pr[0];
						sigmapr[1] = sigmapr[1] / pr[1];
						sigmapr[6] = sigmapr[6] / pr[6];
						pr[0] = log(pr[0]);
						pr[1] = log(pr[1]);
						pr[6] = log(pr[6]);
					}
					break;
				case 'L':
					sigmapr[0] = sigmapr[0] / pr[0];
					sigmapr[1] = sigmapr[1] / pr[1];
					sigmapr[4] = sigmapr[4] / pr[4];
					sigmapr[5] = sigmapr[5] / pr[5];
					pr[0] = log(pr[0]);
					pr[1] = log(pr[1]);
					pr[4] = log(pr[4]);
					pr[5] = log(pr[5]);
					//// If q>1, exchange the two masses
					//if (pr[1] > 0) {
					//	pr[1] = -pr[1];
					//	pr[3] += M_PI;
					//}
					// Put alpha in [0,2\pi]
					pr[3] -= floor(pr[3] / (2 * M_PI)) * 2 * M_PI;
					if (modelcode[1] == 'S') {
						//// If u0<0, use parametrization with u0>0
						//if (pr[2] < 0) {
						//	pr[3] = 2 * M_PI - pr[3];
						//	pr[2] = -pr[2];
						//}
					}
					break;
				}

				fclose(f);

				//					printf("\n%s %d",str2file.name,i1);
				// Check for extremely negative blending
				//for (int j = nps; j < nps + 2 * nfil; j += 2) {
				//	if (pr[j] < -pr[j + 1] * 0.5 /*|| pr[j + 1]<1.e-4*pr[j]*/) c0 = -1;
				//}

				if (c0 > 0) {
					if (nmod) {
						if (c0 < bumperlist->Amp) {
							scanbumper = bumperlist;
							bumperlist = new bumper(pr, nps);
							bumperlist->SetCovariance(Cov, dof / c0);
							strcpy(bumperlist->modelcode, (char *) (itr).path().filename().string().c_str());
							bumperlist->il = i1;
							bumperlist->Amp = c0;
							bumperlist->next = scanbumper;
						}
						else {
							scanbumper = bumperlist;
							while ((scanbumper->next) && (scanbumper->next->Amp < c0)) scanbumper = scanbumper->next;
							scanbumper2 = new bumper(pr, nps);
							scanbumper2->SetCovariance(Cov, dof / c0);
							strcpy(scanbumper2->modelcode, (char*)(itr).path().filename().string().c_str());
							scanbumper2->il = i1;
							scanbumper2->Amp = c0;
							scanbumper2->next = scanbumper->next;
							scanbumper->next = scanbumper2;
						}
					}
					else {
						bumperlist = new bumper(pr, nps);
						bumperlist->SetCovariance(Cov, dof / c0);
						strcpy(bumperlist->modelcode, (char*)(itr).path().filename().string().c_str());
						bumperlist->il = i1;
						bumperlist->Amp = c0;
					}
					nmod++;
				}
			}

			current_path("..");
		}
	}
	if (nmod == 0) {
		printf("\n! No models for this class");
		return 0;
	}
	printf("\nnmod: %d", nmod);
	c1 = (bumperlist) ? bumperlist->Amp : 1.e100;
	printf("\nNumber of initial models = %d\nMinimum chi square = %lf\n", nmod, c1);
	chithr = c1 + c1 / dof * chifac * sqrt(2 * dof);
	printf("\nThreshold for acceptable models = %lf\n", chithr);

	// Apply reflections and check for redundancy

	nmod = 0;
	for (scanbumper = bumperlist; scanbumper; scanbumper = scanbumper->next) {
		// Use these parameters as reference parameters for comparison with other models
		for (int i = 0; i < nps; i++) {
			pr[i] = scanbumper->p0[i];
		}
		// Check for overlap with all other models
		for (scanbumper2 = bumperlist; scanbumper2 != scanbumper; scanbumper2 = scanbumper2->next) {
			// Build inverse covariance as the inverse of the sum of covariances
			CombineCovariances(scanbumper, scanbumper2, Cov, Curv, nps);

			switch (modelcode[0]) {
				// Single lens
				case 'P':
					facr = 0.;
					for (int i = 0; i < nps; i++) {
						for (int j = 0; j < nps; j++) {
							facr += (pr[i] - (scanbumper2->p0)[i]) * (pr[j] - (scanbumper2->p0)[j]) * Curv[i * nps + j];
						}
					}
					fac = facr;
					break;
				// Binary lens
				case 'L':
					// Check with no reflections
					pr[1] = scanbumper->p0[1];
					pr[2] = scanbumper->p0[2];
					pr[3] = scanbumper->p0[3];
					while (pr[3] - scanbumper2->p0[3] > M_PI) pr[3] -= 2*M_PI;
					while (pr[3] - scanbumper2->p0[3] < -M_PI) pr[3] += 2*M_PI;
					facr = 0.;
					for (int i = 0; i < nps; i++) {
						for (int j = 0; j < nps; j++) {
							facr += (pr[i] - (scanbumper2->p0)[i]) * (pr[j] - (scanbumper2->p0)[j]) * Curv[i * nps + j];
						}
					}
					fac = facr;
					if(modelcode[1]=='S') {
						// Repeat with reflection along x1 axis
						pr[3] = 2 * M_PI - scanbumper->p0[3];
						pr[2] = -scanbumper->p0[2];
						while (pr[3] - scanbumper2->p0[3] > M_PI) pr[3] -= 2*M_PI;
						while (pr[3] - scanbumper2->p0[3] < -M_PI) pr[3] += 2*M_PI;

						scanbumper->signCovariance(2);
						scanbumper->signCovariance(3);
						CombineCovariances(scanbumper, scanbumper2, Cov, Curv, nps);

						facr = 0.;
						for (int i = 0; i < nps; i++) {
							for (int j = 0; j < nps; j++) {
								facr += (pr[i] - (scanbumper2->p0)[i]) * (pr[j] - (scanbumper2->p0)[j]) * Curv[i * nps + j];
							}
						}
						if (facr < fac) fac = facr;
						scanbumper->signCovariance(2);
						scanbumper->signCovariance(3);

						// Repeat with reflection between the two masses
						pr[1] = -scanbumper->p0[1];
						pr[3] = M_PI - scanbumper->p0[3];
						pr[2] = -scanbumper->p0[2];
						while (pr[3] - scanbumper2->p0[3] > M_PI) pr[3] -= 2*M_PI;
						while (pr[3] - scanbumper2->p0[3] < -M_PI) pr[3] += 2*M_PI;
						scanbumper->signCovariance(1);
						scanbumper->signCovariance(2);
						scanbumper->signCovariance(3);
						CombineCovariances(scanbumper, scanbumper2, Cov, Curv, nps);
						facr = 0.;
						for (int i = 0; i < nps; i++) {
							for (int j = 0; j < nps; j++) {
								facr += (pr[i] - (scanbumper2->p0)[i]) * (pr[j] - (scanbumper2->p0)[j]) * Curv[i * nps + j];
							}
						}
						if (facr < fac) fac = facr;
						scanbumper->signCovariance(1);
						scanbumper->signCovariance(2);
						scanbumper->signCovariance(3);
					}

					// Repeat with both reflections
					pr[3] = scanbumper->p0[3] - M_PI;
					pr[2] = scanbumper->p0[2];
					pr[1] = -scanbumper->p0[1];
					while (pr[3] - scanbumper2->p0[3] > M_PI) pr[3] -= 2*M_PI;
					while (pr[3] - scanbumper2->p0[3] < -M_PI) pr[3] += 2*M_PI;
					scanbumper->signCovariance(1);
					CombineCovariances(scanbumper, scanbumper2, Cov, Curv, nps);
					facr = 0.;
					for (int i = 0; i < nps; i++) {
						for (int j = 0; j < nps; j++) {
							facr += (pr[i] - (scanbumper2->p0)[i]) * (pr[j] - (scanbumper2->p0)[j]) * Curv[i * nps + j];
						}
					}
					if (facr < fac) fac = facr;
					scanbumper->signCovariance(1);
					break;
					// Binary source
				case 'B':
					// Check with no reflections
					pr[1] = scanbumper->p0[1];
					pr[2] = scanbumper->p0[2];
					pr[3] = scanbumper->p0[3];
					pr[4] = scanbumper->p0[4];
					pr[5] = scanbumper->p0[5];
					facr = 0.;
					for (int i = 0; i < nps; i++) {
						for (int j = 0; j < nps; j++) {
							facr += (pr[i] - (scanbumper2->p0)[i]) * (pr[j] - (scanbumper2->p0)[j]) * Curv[i * nps + j];
						}
					}
					fac = facr;
					if (modelcode[1] == 'S') {
						// Exchange two sources
						pr[1] = -scanbumper->p0[1];
						pr[2] = scanbumper->p0[3];
						pr[3] = scanbumper->p0[2];
						pr[4] = scanbumper->p0[5];
						pr[5] = scanbumper->p0[4];
						scanbumper->signCovariance(1);
						scanbumper->flipCovariance(2, 3);
						scanbumper->flipCovariance(4, 5);
						CombineCovariances(scanbumper, scanbumper2, Cov, Curv, nps);
						facr = 0.;
						for (int i = 0; i < nps; i++) {
							for (int j = 0; j < nps; j++) {
								facr += (pr[i] - (scanbumper2->p0)[i]) * (pr[j] - (scanbumper2->p0)[j]) * Curv[i * nps + j];
							}
						}
						if (facr < fac) fac = facr;
						scanbumper->signCovariance(1);
						scanbumper->flipCovariance(2, 3);
						scanbumper->flipCovariance(4, 5);
					}
			}
			// If models are closer than threshold, the higher chi square model is set beyond acceptance threshold
			if (fac < supfac) {
				if (scanbumper->Amp < scanbumper2->Amp) {
					scanbumper2->Amp = chithr + 1;
				}
				else {
					scanbumper->Amp = chithr + 1;
				}
			}
		}
		nmod++;
		if (nmod - floor(nmod / 50.) * 50. == 0) printf("\nn: %d", nmod);
	}

	// Remove models beyond threshold and all duplicates which have been sent beyond threshold
	if (bumperlist) {
		while (bumperlist->Amp>chithr) {
			nmod--;
			scanbumper = bumperlist->next;
			delete bumperlist;
			bumperlist = scanbumper;
		}
		scanbumper = bumperlist;
		fac = bumperlist->Amp;
		while (scanbumper->next) {
			if (scanbumper->next->Amp>chithr) {
				nmod--;
				scanbumper2 = scanbumper->next->next;
				delete scanbumper->next;
				scanbumper->next = scanbumper2;
			}
			else {
				scanbumper = scanbumper->next;
			}
		}
	}

	printf("\nNumber of independent models with chi square less than \n%.1lf times the minimum = %d\n", chithr, nmod);

	if (nmod > maxmodels) nmod = maxmodels;
	printf("\nModels to be saved = %d\n", nmod);


	// Store best models passing the selection in directory "Models"

	current_path("..");
	create_directory("Models");

	scanbumper = bumperlist;
	for (il = 1; il <= nmod; il++) {
		// Copy corresponding file to directory Models
		auto frompath = path("PreModels") / path(scanbumper->modelcode) / path(to_string(scanbumper->il) + ".txt");
		auto topath = path("Models") / path(string(scanbumper->modelcode) + "-" + to_string(scanbumper->il) + ".txt");
		copy_file(frompath, topath,copy_options::overwrite_existing);
		scanbumper = scanbumper->next;
	}


	current_path("InitCond");

	switch(modelcode[0]) {
	case 'L':
		if (modelcode[1] == 'S') {
			// Preparation of initial conditions for parallax
			printf("\n- Preparing initial conditions for parallax");
			strcpy(filename, "PreInitCondLX.txt");
			if (f = fopen("PreInitCondLX.txt", "r")) {
				g = fopen("InitCondLX.txt", "w");
				fscanf(f, "%d", &np);
				fprintf(g, "0 %d\n", np + 2 * nmod);
				printf("\nNumber of initial conditions: %d", np + 2 * nmod);
				for (int i = 0; i < np; i++) {
					for (int j = 0; j < 9; j++) {
						fscanf(f, "%lg", &pr[j]);
						fprintf(g, "%.10le ", pr[j]);
					}
					fprintf(g, "\n");
				}
				fclose(f);

				scanbumper = bumperlist;
				for (il = 1; il <= nmod; il++) {
					for (int i = 0; i < nps; i++) {
						pr[i] = scanbumper->p0[i];
					}
					fprintf(g, "%.10le %.10le %.10le %.10le %.10le %.10le %.10le 0.0e-003 0.0e-003\n", exp(pr[0]), exp(pr[1]), pr[2], pr[3], exp(pr[4]), exp(pr[5]), pr[6]);
					pr[2] = -pr[2];  // Reflected solution
					pr[3] = 2 * M_PI - pr[3];
					fprintf(g, "%.10le %.10le %.10le %.10le %.10le %.10le %.10le 0.0e-003 0.0e-003\n", exp(pr[0]), exp(pr[1]), pr[2], pr[3], exp(pr[4]), exp(pr[5]), pr[6]);
					scanbumper = scanbumper->next;
				}
				fclose(g);
			}
		}	
		if (modelcode[1] == 'X') {
			printf("\n- Preparing initial conditions for orbital motion");
			if (f = fopen("PreInitCondLO.txt", "r")) {
				g = fopen("InitCondLO.txt", "w");
				fscanf(f, "%d", &np);
				fprintf(g, "0 %d\n", np + nmod);
				printf("\nNumber of initial conditions: %d", np + nmod);
				for (int i = 0; i < np; i++) {
					for (int j = 0; j < 12; j++) {
						fscanf(f, "%lg", &pr[j]);
						fprintf(g, "%.10le ", pr[j]);
					}
					fprintf(g, "\n");
				}
				fclose(f);

				scanbumper = bumperlist;
				for (il = 1; il <= nmod; il++) {
					for (int i = 0; i < nps; i++) {
						pr[i] = scanbumper->p0[i];
					}
					fprintf(g, "%.10le %.10le %.10le %.10le %.10le %.10le %.10le %.10le %.10le 0.0e-006 0.0e-006 1.0e-006\n", exp(pr[0]), exp(pr[1]), pr[2], pr[3], exp(pr[4]), exp(pr[5]), pr[6], pr[7], pr[8]);
					scanbumper = scanbumper->next;
				}
				fclose(g);
			}
		}
		break;
	case 'P':
		if (modelcode[1] == 'S') {
			// Preparation of initial conditions for parallax
			printf("\n- Preparing initial conditions for parallax");
			if (f = fopen("PreInitCondPX.txt", "r")) {				
				g = fopen("InitCondPX.txt", "w");
				fscanf(f, "%d", &np);
				fprintf(g, "0 %d\n", np + 2 * nmod);
				printf("\nNumber of initial conditions: %d", np + 2 * nmod);
				for (int i = 0; i < np; i++) {
					for (int j = 0; j < 6; j++) {
						fscanf(f, "%lg", &pr[j]);
						fprintf(g, "%.10le ", pr[j]);
					}
					fprintf(g, "\n");
				}
				fclose(f);

				scanbumper = bumperlist;
				for (il = 1; il <= nmod; il++) {
					for (int i = 0; i < nps; i++) {
						pr[i] = scanbumper->p0[i];
					}
					fprintf(g, "%.10le %.10le %.10le %.10le 0.0e-003 0.0e-003\n", exp(pr[0]), exp(pr[1]), pr[2], exp(pr[3]));
					// Reflected solution
					fprintf(g, "%.10le %.10le %.10le %.10le 0.0e-003 0.0e-003\n", -exp(pr[0]), exp(pr[1]), pr[2], exp(pr[3]));
					scanbumper = scanbumper->next;
				}
				fclose(g);
			}
			printf("\n- Adding initial conditions for planets");
			if (f = fopen("PreInitCondLS.txt", "r")) {
				int npeaks=0;
				double* peaks;
				g = fopen("InitCondLS.txt", "w");
				fscanf(f, "%d %d", &npeaks, &np);
				fprintf(g, "%d %d\n",npeaks, np + 4 * nmod);
				peaks = (double*)malloc(sizeof(double) * npeaks);

				printf("\nNumber of initial conditions: %d", np + 4 * nmod);
				for (int i = 0; i < npeaks; i++) {
					fscanf(f, "%lg", &peaks[i]);
					fprintf(g, "%le", peaks[i]);
					fscanf(f, "%[^\n]s", &buffer);
					fprintf(g, "%s\n", buffer);
				}
				for (int i = 0; i < np; i++) {
					for (int j = 0; j < 7; j++) {
						fscanf(f, "%lg", &pr[j]);
						fprintf(g, "%.10le ", pr[j]);
					}
					fprintf(g, "\n");
				}
				fclose(f);


				scanbumper = bumperlist;
				for (il = 1; il <= nmod; il++) {
					for (int i = 0; i < nps; i++) {
						pr[i] = scanbumper->p0[i];
					}

					double mindpeak = 1.e100;
					int idpeak = 0;
					for (int ipeak = 0; ipeak < npeaks; ipeak++) {
						double ddpeak = fabs(peaks[ipeak] - pr[2]);
						if (ddpeak < mindpeak) {
							idpeak = ipeak;
							mindpeak = ddpeak;
						}
					}

					double u0 = exp(pr[0]);
					double s0,s;
					for (int ipeak = 0; ipeak < npeaks; ipeak++) {
						if (ipeak == idpeak) continue;
						double q = 0.001;
						double dt = (peaks[ipeak] - pr[2]) / exp(pr[1]);
						double xc = sqrt(u0 * u0 + dt * dt);
						double alpha0 = atan2(u0, -dt), alpha;
						double rho = 0.001; // exp(pr[3] - sqrt(scanbumper->cov[3 * nps + 3]));

						s0 = 0.5 * (sqrt(4 + xc * xc) + xc);
						while (xc < 4 * sqrt(q) / (s0 * s0)) q *= 0.1;
						s = s0 + 4 * sqrt(q) / (s0 * s0);
						alpha = alpha0;
						fprintf(g, "%.10le %.10le %.10le %.10le %.10le %.10le %.10le\n", s, q, u0, alpha, rho,exp(pr[1]),pr[2]);
						s = s0 - 4 * sqrt(q) / (s0 * s0);
						fprintf(g, "%.10le %.10le %.10le %.10le %.10le %.10le %.10le\n", s, q, u0, alpha, rho, exp(pr[1]), pr[2]);

						s = 0.5 * (sqrt(4 + xc * xc) - xc);

						while (xc < 2 * sqrt(q) / s) q *= 0.1;
						alpha = alpha0 + M_PI + asin(fabs(2*sqrt(q)/s)/xc);
						fprintf(g, "%.10le %.10le %.10le %.10le %.10le %.10le %.10le\n", s, q, u0, alpha, rho, exp(pr[1]), pr[2]);
						alpha = alpha0 + M_PI - asin(fabs(2 * sqrt(q) / s) / xc);
						fprintf(g, "%.10le %.10le %.10le %.10le %.10le %.10le %.10le\n", s, q, u0, alpha, rho, exp(pr[1]), pr[2]);
					}
					scanbumper = scanbumper->next;
				}
				fclose(g);
				free(peaks);
			}
		}
		else {
			printf("\n- Adding initial conditions for planets");
			if (f = fopen("PreInitCondLX.txt", "r")) {
				int npeaks = 0;
				double* peaks;
				g = fopen("InitCondLX-2.txt", "w");
				fscanf(f, "%d %d", &npeaks, &np);
				fprintf(g, "%d %d\n", npeaks, np + 4 * nmod);
				peaks = (double*)malloc(sizeof(double) * npeaks);

				printf("\nNumber of initial conditions: %d", np + 4 * nmod);
				for (int i = 0; i < npeaks; i++) {
					fscanf(f, "%lg", &peaks[i]);
					fprintf(g, "%le", peaks[i]);
					fscanf(f, "%[^\n]s", &buffer);
					fprintf(g, "%s\n", buffer);
				}
				for (int i = 0; i < np; i++) {
					for (int j = 0; j < 9; j++) {
						fscanf(f, "%lg", &pr[j]);
						fprintf(g, "%.10le ", pr[j]);
					}
					fprintf(g, "\n");
				}
				fclose(f);

				double mindpeak = 1.e100;
				int idpeak = 0;
				for (int ipeak = 0; ipeak < npeaks; ipeak++) {
					double ddpeak = fabs(peaks[ipeak] - pr[2]);
					if (ddpeak < mindpeak) {
						idpeak = ipeak;
						mindpeak = ddpeak;
					}
				}

				scanbumper = bumperlist;
				for (il = 1; il <= nmod; il++) {
					for (int i = 0; i < nps; i++) {
						pr[i] = scanbumper->p0[i];
					}

					double mindpeak = 1.e100;
					int idpeak = 0;
					for (int ipeak = 0; ipeak < npeaks; ipeak++) {
						double ddpeak = fabs(peaks[ipeak] - pr[2]);
						if (ddpeak < mindpeak) {
							idpeak = ipeak;
							mindpeak = ddpeak;
						}
					}

					double u0 = exp(pr[0]);
					double s0, s;
					for (int ipeak = 0; ipeak < npeaks; ipeak++) {
						if (ipeak == idpeak) continue;
						double q = 0.001;
						double dt = (peaks[ipeak] - pr[2]) / exp(pr[1]);
						double xc = sqrt(u0 * u0 + dt * dt);
						double alpha0 = atan2(u0, -dt), alpha;
						double rho = 0.001; // exp(pr[3] - sqrt(scanbumper->cov[3 * nps + 3]));

						s0 = 0.5 * (sqrt(4 + xc * xc) + xc);
						while (xc < 4 * sqrt(q) / (s0 * s0)) q *= 0.1;
						s = s0 + 4 * sqrt(q) / (s0 * s0);
						alpha = alpha0;
						fprintf(g, "%.10le %.10le %.10le %.10le %.10le %.10le %.10le %.10le %.10le\n", s, q, u0, alpha, rho, exp(pr[1]), pr[2], pr[4], pr[5]);
						s = s0 - 4 * sqrt(q) / (s0 * s0);
						fprintf(g, "%.10le %.10le %.10le %.10le %.10le %.10le %.10le %.10le %.10le\n", s, q, u0, alpha, rho, exp(pr[1]), pr[2], pr[4], pr[5]);

						s = 0.5 * (sqrt(4 + xc * xc) - xc);

						while (xc < 2 * sqrt(q) / s) q *= 0.1;
						alpha = alpha0 + M_PI + asin(fabs(2 * sqrt(q) / s) / xc);
						fprintf(g, "%.10le %.10le %.10le %.10le %.10le %.10le %.10le %.10le %.10le\n", s, q, u0, alpha, rho, exp(pr[1]), pr[2], pr[4], pr[5]);
						alpha = alpha0 + M_PI - asin(fabs(2 * sqrt(q) / s) / xc);
						fprintf(g, "%.10le %.10le %.10le %.10le %.10le %.10le %.10le %.10le %.10le\n", s, q, u0, alpha, rho, exp(pr[1]), pr[2], pr[4], pr[5]);
					}
					scanbumper = scanbumper->next;
				}

				fclose(g);
				free(peaks);

				if (g = fopen("InitCondLS.txt", "r")) {
					fclose(g);
					remove("PreInitCondLX.txt");
					rename("InitCondLX-2.txt", "PreInitCondLX.txt");
				}
				else {
					rename("InitCondLX-2.txt", "InitCondLX.txt");
				}
			}
			else {
				if (f = fopen("PreInitCondLO.txt", "r")) {
					int npeaks = 0;
					double* peaks;
					g = fopen("InitCondLO-2.txt", "w");
					fscanf(f, "%d %d", &npeaks, &np);
					fprintf(g, "%d %d\n", npeaks, np + 4 * nmod);
					peaks = (double*)malloc(sizeof(double) * npeaks);

					printf("\nNumber of initial conditions: %d", np + 4 * nmod);
					for (int i = 0; i < npeaks; i++) {
						fscanf(f, "%lg", &peaks[i]);
						fprintf(g, "%le", peaks[i]);
						fscanf(f, "%[^\n]s", &buffer);
						fprintf(g, "%s\n", buffer);
					}
					for (int i = 0; i < np; i++) {
						for (int j = 0; j < 12; j++) {
							fscanf(f, "%lg", &pr[j]);
							fprintf(g, "%.10le ", pr[j]);
						}
						fprintf(g, "\n");
					}
					fclose(f);

					double mindpeak = 1.e100;
					int idpeak = 0;
					for (int ipeak = 0; ipeak < npeaks; ipeak++) {
						double ddpeak = fabs(peaks[ipeak] - pr[2]);
						if (ddpeak < mindpeak) {
							idpeak = ipeak;
							mindpeak = ddpeak;
						}
					}

					scanbumper = bumperlist;
					for (il = 1; il <= nmod; il++) {
						for (int i = 0; i < nps; i++) {
							pr[i] = scanbumper->p0[i];
						}

						double mindpeak = 1.e100;
						int idpeak = 0;
						for (int ipeak = 0; ipeak < npeaks; ipeak++) {
							double ddpeak = fabs(peaks[ipeak] - pr[2]);
							if (ddpeak < mindpeak) {
								idpeak = ipeak;
								mindpeak = ddpeak;
							}
						}

						double u0 = exp(pr[0]);
						double s0, s;
						for (int ipeak = 0; ipeak < npeaks; ipeak++) {
							if (ipeak == idpeak) continue;
							double q = 0.001;
							double dt = (peaks[ipeak] - pr[2]) / exp(pr[1]);
							double xc = sqrt(u0 * u0 + dt * dt);
							double alpha0 = atan2(u0, -dt), alpha;
							double rho = 0.001; // exp(pr[3] - sqrt(scanbumper->cov[3 * nps + 3]));

							s0 = 0.5 * (sqrt(4 + xc * xc) + xc);
							while (xc < 4 * sqrt(q) / (s0 * s0)) q *= 0.1;
							s = s0 + 4 * sqrt(q) / (s0 * s0);
							alpha = alpha0;
							fprintf(g, "%.10le %.10le %.10le %.10le %.10le %.10le %.10le %.10le %.10le 0.0001 0.0001 0.0001\n", s, q, u0, alpha, rho, exp(pr[1]), pr[2], pr[4], pr[5]);
							s = s0 - 4 * sqrt(q) / (s0 * s0);
							fprintf(g, "%.10le %.10le %.10le %.10le %.10le %.10le %.10le %.10le %.10le 0.0001 0.0001 0.0001\n", s, q, u0, alpha, rho, exp(pr[1]), pr[2], pr[4], pr[5]);

							s = 0.5 * (sqrt(4 + xc * xc) - xc);

							while (xc < 2 * sqrt(q) / s) q *= 0.1;
							alpha = alpha0 + M_PI + asin(fabs(2 * sqrt(q) / s) / xc);
							fprintf(g, "%.10le %.10le %.10le %.10le %.10le %.10le %.10le %.10le %.10le 0.0001 0.0001 0.0001\n", s, q, u0, alpha, rho, exp(pr[1]), pr[2], pr[4], pr[5]);
							alpha = alpha0 + M_PI - asin(fabs(2 * sqrt(q) / s) / xc);
							fprintf(g, "%.10le %.10le %.10le %.10le %.10le %.10le %.10le %.10le %.10le 0.0001 0.0001 0.0001\n", s, q, u0, alpha, rho, exp(pr[1]), pr[2], pr[4], pr[5]);
						}
						scanbumper = scanbumper->next;
					}

					fclose(g);
					free(peaks);

					rename("InitCondLO-2.txt", "InitCondLO.txt");					
				}
			}
		}
		break;
	case 'B':
		if (modelcode[1] == 'S') {
			// Preparation of initial conditions for xallarap
			printf("\n- Preparing initial conditions for xallarap");
			if (f = fopen("PreInitCondBO.txt", "r")) {
				g = fopen("InitCondBO.txt", "w");
				fscanf(f, "%d", &np);
				fprintf(g, "0 %d\n", np + 2 * nmod);
				printf("\nNumber of initial conditions: %d", np + 2 * nmod);
				for (int i = 0; i < np; i++) {
					for (int j = 0; j < 10; j++) {
						fscanf(f, "%lg", &pr[j]);
						fprintf(g, "%.10le ", pr[j]);
					}
					fprintf(g, "\n");
				}
				fclose(f);

				scanbumper = bumperlist;
				for (il = 1; il <= nmod; il++) {
					for (int i = 0; i < nps; i++) {
						pr[i] = scanbumper->p0[i];
					}
					double u01 = pr[2];
					double u02 = pr[3];
					double tE = exp(pr[0]);
					double qs = exp(pr[1]/4);
					fprintf(g, "%.10le %.10le %.10le %.10le %.10le %.10le %.10le %.10le %.10le %.10le\n", u01, pr[4], tE, exp(pr[6]), (pr[5] - pr[4]) / tE / (1 + qs) * qs, (-u02 + u01) / (1 + qs) * qs, 0.000001, 0.0001, 0.00001, qs);
					// Reflected solution
					fprintf(g, "%.10le %.10le %.10le %.10le %.10le %.10le %.10le %.10le %.10le %.10le\n", u01, pr[4], tE, exp(pr[6]), (pr[5] - pr[4]) / tE / (1 + qs) * qs, (u02 + u01) / (1 + qs) * qs, 0.000001, 0.0001, 0.00001, qs);
					scanbumper = scanbumper->next;
				}
				fclose(g);
			}
		}
		break;
	}

	printf("\n\n- Closing...");

	scanbumper = bumperlist;
	while (scanbumper) {
		scanbumper2 = scanbumper->next;
		delete scanbumper;
		scanbumper = scanbumper2;
	}

	free(pr);
	free(sigmapr);
	free(Cov);
	free(Curv);

	return 0;
}
