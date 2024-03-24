// Finalizer.cpp : main project file.

#define _USE_MATH_DEFINES
#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <malloc.h>
#include <regex>
#include <filesystem>
#include "bumper.h"

using namespace std;
using namespace std::filesystem;

const double failthr = 20000000.0; // threshold for chisquare/dof for declaring failure
const int ncategories = 7;

int main(int argc, char *argv[]){

	char eventname[512]="";
	char filename[30]="";
	char command[256], buffer[256];
	double value;
	int nfil, * satel;
	double *t,*y,*w,*pr,*sigmapr,thsigma;
    double thrs[10] = { 0,36.,40.0872,43.4518,46.4625,49.2497,51.878 ,54.3854 ,56.7964 ,59.1282 }; // thresholds at 6 sigma for n more parameters
	string modelcodes[ncategories] = { "PS","PX","BS","BO","LS","LX","LO" };
	string modelnames[ncategories] = { "Single-Lens-Single-Source","Single-Lens-Single-Source with parallax",\
									"Binary Source","Binary Source with xallarap",\
									""," with parallax"," with orbital motion" };
	int npss[ncategories] = { 4,6,7,10,7,9,12 };
	double chis[ncategories];
	double cmin, c0, c1, c2, cflat;
	double chiblp=1.e100, chiblb = 1.e100;
	int mn;
	char final[2000],stringa[1024];
	int flag,nmod,ngoodmod,*filter,np;
	FILE *f,*g;
	bumper* bumperlist=0, * scanbumper, * scanbumper2;


	printf("******************************************\n");
	printf("**********     Finalizer    *********\n");
	printf("******************************************\n\n\n");
	printf("This program proposes a final interpretation for the event\n\n");

// Directory preliminaries. Reads event name from arguments.

	auto exedir = current_path();

	if (argc > 1) {
		strcpy(eventname, argv[1]);
	}
	else {
		printf("\n\nEvent name? ");
		scanf("%s", eventname);
	}

	printf("\n\n- Event: %s\n", eventname);

	current_path(eventname);
	
	if (exists("Nature.txt")) {
		printf("\n\nEvent already finalized");
		return 0;
	}

	if (exists("ini")) {
		current_path("ini");
		f = fopen("Finalizer.ini", "r");
		if (f != 0) {
			printf("\n\n- Reading options in Finalizer.ini");
			while (!feof(f)) {
				int red = fscanf(f, "%s %s %lf", command, buffer, &value);
				if (red < 1) {
					command[0] = 0;
					//if (red != 0) {
					//	printf("\n\n!!! Bad command in Reader.ini");
					//	return -1;
					//};
				}
				//if (strcmp(command, "maxmodels") == 0) {
				//	maxmodels = value;
				//}

			}
			fclose(f);
		}
		else {
			printf("\n\n- Default options:");
		}
	}

	/*printf("\nNumber of sigmas used to declare overlap between two models: %lf", supfac);
	printf("\nNumber of sigmas in chi square distribution for model acceptance: %lf", chifac);
	printf("\nMaximum number of models reported: %d", maxmodels);*/


// Read curve to fit

	current_path(eventname);

	printf("\n\nReading data\n");

	f=fopen("LCToFit.txt","r");
    fscanf(f,"%d",&np);
	filter=(int *) malloc(sizeof(int)*np);
	satel=(int *) malloc(sizeof(int)*np);
	t=(double *) malloc(sizeof(double)*np);
	y=(double *) malloc(sizeof(double)*np);
	w=(double *) malloc(sizeof(double)*np);

	nfil=1;
	for(int i=0;i<np;i++){
		fscanf(f,"%d %lf %lf %lf %d",&(filter[i]),&(t[i]),&(y[i]),&(w[i]),&(satel[i]));
		if((i!=0)&&(filter[i]!=filter[i-1])){
			nfil++;
		}
		w[i]=1/(w[i]);
	}
	fclose(f);


	pr=(double *) malloc(sizeof(double)*(12+2*nfil));
	sigmapr = (double *)malloc(sizeof(double)*(12 + 2 * nfil));

// Calculate flat chi square

	c0 = c1=c2=cflat=0;
	nfil = 0;
	for (int i = 0; i < np; i++) {
		if (filter[i] == nfil) {
			c0 += w[i] * w[i];
			c1 += y[i] * w[i] * w[i];
			c2 += y[i] * y[i] * w[i] * w[i];
		}
		else {
			nfil++;
			cflat += c2 - c1 * c1 / c0;
			c0 = c1 = c2 = 0;
		}
	}
	nfil++;
	cflat += c2 - c1 * c1 / c0;

// Compare all models and draw conclusions

	printf("\n\n********************");
	printf("\n   Finalization   ");
	printf("\n********************\n");

	current_path(eventname);

	auto selectmodelsname = path("FinalModels");
	create_directory(selectmodelsname);


	g=fopen("Nature.txt","w");

	current_path("Models");

// Load chi square of best models of each type

	printf("\n- Reading models\n\n");

	nmod = 0;
	for (int icat = 0; icat < ncategories; icat++) {
		chis[icat] = 1.e100;
		auto searchstring = regex(modelcodes[icat] + "[0-9]*-[0-9]*.txt");
		for (auto const& itr : directory_iterator(".")) {
			string curfile = (itr).path().filename().string();
			if (regex_match(curfile, searchstring)) {
				strcpy(filename, curfile.c_str());
				f = fopen(filename, "r");
				for (int j = 0; j < npss[icat] + 2 * nfil; j++) {
					fscanf(f, "%le", &(pr[j]));
				}
				fscanf(f, "%le", &(c0));
				fclose(f);

				if (c0 < chis[icat]) chis[icat] = c0;

				if (modelcodes[icat][0] == 'L') {
					if (pr[1] <= 0.03) {
						if (c0 < chiblp) chiblp = c0;
					}
					else {
						if (c0 < chiblb) chiblb = c0;
					}
				}
				if (nmod) {
					if (c0 < bumperlist->Amp) {
						scanbumper = bumperlist;
						bumperlist = new bumper(pr, npss[icat]);
						strcpy(bumperlist->modelcode, (char*)(itr).path().filename().string().c_str());
						bumperlist->il = icat;
						bumperlist->Amp = c0;
						bumperlist->next = scanbumper;
					}
					else {
						scanbumper = bumperlist;
						while ((scanbumper->next) && (scanbumper->next->Amp < c0)) scanbumper = scanbumper->next;
						scanbumper2 = new bumper(pr, npss[icat]);
						strcpy(scanbumper2->modelcode, (char*)(itr).path().filename().string().c_str());
						scanbumper2->il = icat;
						scanbumper2->Amp = c0;
						scanbumper2->next = scanbumper->next;
						scanbumper->next = scanbumper2;
					}
				}
				else {
					bumperlist = new bumper(pr, npss[icat]);
					strcpy(bumperlist->modelcode, (char*)(itr).path().filename().string().c_str());
					bumperlist->il = icat;
					bumperlist->Amp = c0;
				}
				nmod++;
			}
		}
		printf("%s: ", modelcodes[icat].c_str());
		fprintf(g, "%s: ", modelcodes[icat].c_str());
		if (chis[icat] < 1.e99) {
			printf("%lf\n", chis[icat]);
			fprintf(g,"%lf\n", chis[icat]);
		}
		else {
			printf("N/A\n");
			fprintf(g,"N/A\n");
		}
	}


	//// Check that there are no weird chi squares
	//for (int i = 0; i<9; i++) {
	//	if (chis[i]<np/2) {
	//		chis[i] = 1.e100;
	//	}
	//}

	// Minimum chi square

	mn=0;
	cmin=bumperlist->Amp;

	thsigma = cmin + cmin / np * sqrt(2 * np);
	for (int i = 0; i < 10; i++) {
		thrs[i] *= cmin / np;
	}

	if (chiblp < 1.e99) {
		printf("\nBestPlanetary: %lf", chiblp);
		fprintf(g, "\nBestPlanetary: %lf", chiblp);
	}
	else {
		printf("\nBestPlanetary: N/A");
		fprintf(g, "\nBestPlanetary: N/A");
	}

	if (chiblb < 1.e99) {
		printf("\nBestBinary: %lf\n\n", chiblb);
		fprintf(g, "\nBestBinary: %lf\n\n", chiblb);
	}
	else {
		printf("\nBestBinary: N/A\n\n");
		fprintf(g, "\nBestBinary: N/A\n\n");
	}

	if (chiblb + thsigma - cmin < chiblp) {
		strcpy(stringa, "Binary lens");
	}
	else {
		if (chiblp + thsigma - cmin < chiblb) {
			strcpy(stringa, "Planetary lens");
		}
		else {
			strcpy(stringa, "Binary or Planetary lens");
		}
	}

	// Making assessment. Each category is tested against nested ones

	if (chis[1] > chis[0] - thrs[npss[1] - npss[0]]) {
		chis[1] = 1.e100;
	}
	if(chis[2] > chis[0] - thrs[npss[2] - npss[0]]) {
		chis[2] = 1.e100;
	}
	if (chis[3] > chis[0] - thrs[npss[3] - npss[0]] || chis[3] > chis[1] - thrs[npss[3] - npss[1]] || chis[3] > chis[2] - thrs[npss[3] - npss[2]]) {
		chis[3] = 1.e100;
	}
	if (chis[4] > chis[0] - thrs[npss[4] - npss[0]]) {
		chis[4] = 1.e100;
	}
	if (chis[5] > chis[0] - thrs[npss[5] - npss[0]] || chis[5] > chis[1] - thrs[npss[5] - npss[1]] || chis[5] > chis[4] - thrs[npss[5] - npss[4]]) {
		chis[5] = 1.e100;
	}
	if (chis[6] > chis[0] - thrs[npss[6] - npss[0]] || chis[6] > chis[1] - thrs[npss[6] - npss[1]] || chis[6] > chis[4] - thrs[npss[6] - npss[4]] || chis[6] > chis[5] - thrs[npss[6] - npss[5]]) {
		chis[6] = 1.e100;
	}

	// If more complicated category has survived, all nested categories are removed

	if (chis[6] < 1.e99) {
		chis[5] = chis[4] = chis[1] = chis[0] = 1.e100;
	}
	if (chis[5] < 1.e99) {
		chis[4] = chis[1] = chis[0] = 1.e100;
	}
	if (chis[4] < 1.e99) {
		chis[0] = 1.e100;
	}
	if (chis[3] < 1.e99) {
		chis[2] = chis[1] = chis[0] = 1.e100;
	}
	if (chis[2] < 1.e99) {
		chis[0] = 1.e100;
	}
	if (chis[1] < 1.e99) {
		chis[0] = 1.e100;
	}

	// Models of discarded categories or with chi square higher than threshold are removed

	for (scanbumper = bumperlist; scanbumper; scanbumper = scanbumper->next) {
		if (scanbumper->Amp > thsigma || chis[scanbumper->il] > 1.e99) scanbumper->modelcode[0] = 'N';
	}

	// Counting good models

	ngoodmod = 0;
	for (scanbumper = bumperlist; scanbumper; scanbumper = scanbumper->next) {
		if (scanbumper->modelcode[0] != 'N') ngoodmod++;
	}

	// Decide designation

	if (cflat < thsigma) {
		strcpy(final, "Flat");
	}else{
		if (cmin > np * failthr) {
			strcpy(final, "Uncertain: variable or failed");
		}
		else {
			strcpy(final, "Successful: ");
			flag = 0;
			for (int icat = 0; icat < ncategories; icat++) {
				for (scanbumper = bumperlist; scanbumper; scanbumper = scanbumper->next) {
					if (scanbumper->il == icat && scanbumper->modelcode[0] != 'N') {
						if (flag > 0) {
							strncat(final, " or ", 60);
						}
						else {
							flag = 1;
						}
						if (icat >= 4) strncat(final, stringa, 60);
						strncat(final, modelnames[icat].c_str(), 60);
						break;
					}
				}
			}
		}
	}

	printf("%s\n", final);	
// Writing file Nature.txt with the conclusion drawn before
	fprintf(g,"%s\n",final);

	printf("----\n");
	fprintf(g,"----\n");


	printf("Number of alternative models: %d\n\nchisquare   model\n", ngoodmod);
	fprintf(g, "Number of alternative models: %d\n\nchisquare   model\n", ngoodmod);


// Writing the names of the files containing alternative models in Nature.txt
	current_path("..");
	for (scanbumper = bumperlist; scanbumper;scanbumper=scanbumper->next) {
		if (scanbumper->modelcode[0] != 'N') {
			fprintf(g, "%lf %s\n", scanbumper->Amp, scanbumper->modelcode);
			printf("%lf %s", scanbumper->Amp, scanbumper->modelcode);
			copy_file(path("Models") / path(string(scanbumper->modelcode)), path("FinalModels") / path(string(scanbumper->modelcode)),copy_options::overwrite_existing);
		}
	}
	fclose(g);



	printf("\n\n- Done");

//	Sleep(5000l);
	
	scanbumper = bumperlist;
	while (scanbumper) {
		scanbumper2 = scanbumper->next;
		delete scanbumper;
		scanbumper = scanbumper2;
	}

	free(pr);
	free(sigmapr);
	free(filter);
	free(satel);
	free(t);
	free(y);
	free(w);


    return 0;
}
