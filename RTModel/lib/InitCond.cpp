// InitCond.cpp : main project file.
// This program finds the peaks in the data and sets initial conditions for fitting

#define _CRT_SECURE_NO_WARNINGS
#define _USE_MATH_DEFINES
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <malloc.h>
#include <regex>
#include <filesystem>

using namespace std;
using namespace std::filesystem;

#ifdef _WIN32
char systemslash = '\\';
#else
char systemslash = '/';
#endif

// Main global parameters
int nobspeaks=2; // Number of peaks in the observed light curve to be considered for setting initial conditions.
double sigmathr=5.; // Number of sigmas for spline approximation
double peakthr = 10.; // Number of sigmnas necessary for a deviation to be identified as a maximum or a minimum.
int maxoldmodels = 4; // Maximum number of old models to include in new run as initial conditions
bool override = false; // Override peak identification and manually set peak times
bool nostatic = false; // No static models will be calculated.
bool noparallax = false; // Only orbital motion models will be calculated.
int usesatellite = 0; // Satellite to be used for initial conditions. Ground telescopes by default.

// Structure datapoint: stores time (t), flux (y), error (yerr), significance (sig), time range for the uncertainty (tl-tr)
struct datapoint{
	datapoint *prev,*next;
	int i;
	double t,y,tl,tr,yerr,sig;
};

// Structure dataset: a list of datapoints with obvious functions
struct dataset{
	dataset *prev,*next;
	datapoint *first,*last;
	int length;

	dataset();
	~dataset();
	datapoint *addpoint(int,double,double,double);
	void addpoint(double *,double *,double *,int,int,int);
	void remove(datapoint *);
	void clone(datapoint *);
};

int main(int argc, char *argv[])
{
	char eventname[256] = "";
	regex filebest;
	char fileinit[256]="";
	char command[256], buffer[256];
	double value,value2;
	double tv,yv,errv,t0,tE,t1,t2,tasy,mean,sw;
	double *t,*y,*err,*tt,*yy,*eerr;
	double maxdev, curdev, asydev=0,w1,w2;
	FILE *f;
	int flag,np,nps,npc,*nnp,nfil,ifil,dn,satellite;
	double fint;
	dataset **peaklist, *cpeaks, *newpeaks;
	datapoint *p,*pm,*pmm,*pl,*pr,*pasy=0, *highestpeak, *minimum, *startsection, *endsection, *sectionpeak;


// Directory preliminaries. Reads event name from arguments.

	printf("******************************************\n");
	printf("*************     InitCond      **********\n");
	printf("******************************************\n\n\n");
	printf("***This program finds peaks and sets initial conditions for fitting\n\n");

	auto exedir = current_path();

	if (argc > 1) {
		strcpy(eventname, argv[1]);
	}
	else {
		printf("\n\nEvent name? ");
		scanf("%s", eventname);
		//sprintf(eventname, "WDC10193");
	}

	printf("\n\n- Event: %s\n", eventname);

	current_path(eventname);


	if (exists("InitCond")) {
		current_path("InitCond");
		auto searchstring = regex("Init.*");
		for (auto const& itr : directory_iterator(".")) {
			string curfile = (itr).path().filename().string();
			if (regex_match(curfile, searchstring)) {
				printf("\n\nInitial conditions already selected");
				return 0;
			}
			itr;
		}
	}
	current_path(eventname);

// Look for previous runs to include previously found best models

	auto searchstring = regex("run.*");
	int lastrun = -1;
	string runstring = "";
	for (auto const& itr : directory_iterator(".")) {
		string curfile = (itr).path().filename().string();
		if (regex_match(curfile, searchstring)) {
			int currun = atoi(curfile.c_str() + 4);
			if (currun > lastrun) {
				lastrun = currun;
				runstring = curfile;
			}
		}
	}

// Reading InitCond.ini and set parameters accordingly

	if (exists("ini")) {
		current_path("ini");
		f = fopen("InitCond.ini", "r");
		if (f != 0) {
			printf("\n\n- Reading options in InitCond.ini");
			while (!feof(f)) {
				int red = fscanf(f, "%s %s %lf", command, buffer, &value);
				if (red < 1) {
					command[0] = 0;
					//if (red != 0) {
					//	printf("\n\n!!! Bad command in Reader.ini");
					//	return -1;
					//};
				}
				if (strcmp(command, "override") == 0) {
					fscanf(f, " %lf", &value2);
					override = true;
					t1 = value;
					t2 = value2;
				}
				if (strcmp(command, "npeaks") == 0) {
					nobspeaks = (int) value;
				}
				if (strcmp(command, "peakthreshold") == 0) {
					peakthr = value;
				}
				if (strcmp(command, "oldmodels") == 0) {
					maxoldmodels = (int) value;
				}
				if (strcmp(command, "nostatic") == 0) {
					nostatic = true;
				}
				if (strcmp(command, "noparallax") == 0) {
					noparallax = true;
					nostatic = true;
				}
				if (strcmp(command, "usesatellite") == 0) {
					usesatellite = (int) value;
				}
			}
			fclose(f);
		}
		else {
			printf("\n\n- Default options:");
		}
		current_path(eventname);
	}

	printf("\nNumber of peaks to be considered in light curve %d", nobspeaks);
	printf("\nThreshold for peak identification %lf", peakthr);
	sigmathr = peakthr * 0.5;
	printf("\nOld models to be included in next run %d", maxoldmodels);
	if(noparallax) printf("\nOnly orbital motion models will be initialized");
	if (nostatic && !noparallax) printf("\nStatic models will not be initialized");

	// Check for Override. If you do not like the initial conditions found by this algorithm, 
	// you may optopnally override them by writing the times of two peaks in InitCond.ini.
	
	if(override){
		// Use hand-made initial conditions if desired.
		printf("\n- Overriding peak identification -\n");
		newpeaks = new dataset;
		pl = new datapoint;
		pl->t = t1;
		pr = new datapoint;
		pr->t = t2;
		pl->next = pr;
		pl->prev = 0;
		pr->prev=pl;
		pr->next = 0;
		newpeaks->first=pl;
		newpeaks->length = 2;
		newpeaks->last = pr;
	}
	else { // Normal flow if override is not used

		/* Reading the curve to fit */

		printf("\n\n- Reading data \n");

		f = fopen("LCToFit.txt", "r");
		fscanf(f, "%d", &np);

		nfil = 1;
		dn = 0;
		for (int i = 0; i < np; i++) {
			fscanf(f, "%d %lf %lf %lf %d", &(ifil), &(tv), &(yv), &(errv), &(satellite));
			if (satellite != usesatellite) {
				np--;
				i--;
			}
			else {
				if ((i != 0) && (ifil != dn)) {
					nfil++;
				}
			}
			dn = ifil;
		}
		fclose(f);
		np++;

		tt = (double*)malloc(sizeof(double) * np);
		yy = (double*)malloc(sizeof(double) * np);
		eerr = (double*)malloc(sizeof(double) * np);
		nnp = (int*)calloc(sizeof(int), nfil);
		peaklist = (dataset**)malloc(sizeof(dataset*) * nfil);

		f = fopen("LCToFit.txt", "r");
		fscanf(f, "%d", &np);

		nfil = 0;
		dn = 0;
		for (int i = 0; i < np; i++) {
			fscanf(f, "%d %lf %lf %lf %d", &(ifil), &(tt[i]), &(yy[i]), &(eerr[i]), &(satellite));
			if (satellite > 0) {
				np--;
				i--;
			}
			else {
				if ((i != 0) && (ifil != dn)) {
					nfil++;
				}
				nnp[nfil]++;
			}
			dn = ifil;
		}
		nfil++;
		fclose(f);

		// End of data reading. Here is the relevant information:
		// nfil = number of datasets (one per telescope and filter)
		// nnp[nfil] = number of points in each dataset
		// np = total number of points
		// tt[np] = array of all times of observations. Data are sorted by dataset and then by time.
		// yy[np] = Corresponding fluxes
		// eerr[np] = Corresponding errors
		// peaklist[nfil] will contain the list of peaks for each dataset.

		// tasy = time of more relevant asymmetry (to be calculated later)
		// asydev = relevance of the asymmetry
		fint = -1.e100;
		tasy = 0;
		asydev = 0;
		f = fopen("spline.txt", "w");
		fclose(f);
		///////////////////////////////////
		// Here we start the main algorithm
		printf("\n-- Processing\n");

		// Cycle on each dataset
		for (int ic = 0; ic < nfil; ic++) {
			printf("\n- Dataset: %d", ic);
			// Find starting point in the arrays for this dataset
			npc = 0;
			for (int i = 0; i < ic; i++) {
				npc += nnp[i];
			}
			// t, y and err point to the first point in the dataset and will be used henceforth
			t = &(tt[npc]);
			y = &(yy[npc]);
			err = &(eerr[npc]);
			npc = nnp[ic]; // npc = Number of points in this dataset
			fint = (t[npc - 1] > fint) ? t[npc - 1] : fint; // fint = Last observation time, useful for setting initial conditions on incomplete events

			//// This block is used to skip all data points in the previous years for the analysis.
			//	// It is useful to avoid picking spurious peaks in the baseline of previous years.
			//	// Feel free to remove it or elaborate on this (e.g. by selecting only the peak season)
			//	dn = 0;
			//	while (t[dn]<5197 + (yr - 2010)*365.25) dn++;
			//	t = &(t[dn]);
			//	y = &(y[dn]);
			//	err = &(err[dn]);
			//	npc -= dn;

			//	printf("\n+ Smoothing data\n");
			//// Each data point is replaced by the weighted mean with the previous and following points.
			//	dept1=t[0];
			//	depy1=y[0];
			//	deperr1=err[0];
			//	sw = mean=0;
			//	for(int i=1;i<=npc-2;i++){
			//		dept2=t[i];
			//		depy2=y[i];
			//		deperr2=err[i];
			//		err[i]=pow(1/(deperr1*deperr1)+1/(deperr2*deperr2)+1/(err[i+1]*err[i+1]),-0.5);
			//		y[i]=(depy1/(deperr1*deperr1)+depy2/(deperr2*deperr2)+y[i+1]/(err[i+1]*err[i+1]))*err[i]*err[i];
			//		t[i]=(dept1/(deperr1*deperr1)+dept2/(deperr2*deperr2)+t[i+1]/(err[i+1]*err[i+1]))*err[i]*err[i];
			//		dept1=dept2;
			//		depy1=depy2;
			//		deperr1=deperr2;
			//		sw += 1 / (deperr2*deperr2);
			//		mean += depy2 / (deperr2*deperr2);
			//	}
			//	mean /= sw; // mean = mean flux for this dataset.

				// Calculation of the mean for the dataset

			sw = mean = 0;
			for (int i = 0; i <= npc - 1; i++) {
				sw += 1 / (err[i] * err[i]);
				mean += y[i] / (err[i] * err[i]);
			}
			mean /= sw;

			// Spline approximation 
			printf("\n- Spline approximation");
			peaklist[ic] = new dataset;
			cpeaks = peaklist[ic];
			cpeaks->addpoint(-1, t[0] - 1, mean, 1 / sqrt(sw));
			cpeaks->addpoint(npc, t[npc - 1] + 1, mean, 1 / sqrt(sw));
			// cpeaks = spline approximating the light curve. We start from a straight horizontal line at the mean flux.
			int iter = 0, flag = 0, imax, jmax;
			datapoint* curpeak;
			// Cycle until maxdev becomes less than sigmathr or all points have been added.
			while (iter < npc && (maxdev > sigmathr || flag < 2)) {
				iter++;
				imax = 0;
				maxdev = 0;
				curpeak = cpeaks->first;
				// Cycle on all points in the dataset
				for (int i = 1; i < npc - 1; i++) {
					if (t[i] > curpeak->next->t) curpeak = curpeak->next;
					// curdev = deviation of the current point from the spline
					w1 = (t[i] - curpeak->t) / (curpeak->next->t - curpeak->t);
					w2 = 1 - w1;
					curdev = fabs((y[i] - curpeak->y - (curpeak->next->y - curpeak->y) * w1) / sqrt(err[i] * err[i] + curpeak->next->yerr * curpeak->next->yerr * w1 * w1 + curpeak->yerr * curpeak->yerr * w2 * w2));
					if (curdev > maxdev) {
						maxdev = curdev;
						imax = i;
					}
				}
				// point with maximum deviation is added to the spline
				cpeaks->addpoint(imax, t[imax], y[imax], err[imax]);

				imax = 10000;
				jmax = 0;
				//				maxdev = 1.e100;
							// Significance of each peak in the spline is updated
				for (curpeak = cpeaks->first->next; curpeak->next; curpeak = curpeak->next) {
					// curdev = deviation of current point in the spline from the spline without this point
					w1 = (curpeak->t - curpeak->prev->t) / (curpeak->next->t - curpeak->prev->t);
					w2 = 1 - w1;
					curdev = curpeak->y - curpeak->prev->y - (curpeak->next->y - curpeak->prev->y) * w1;
					// Only for positive deviations (concavities in the spline) we calculate the significance
					if (curdev > 0) {
						flag = (flag > jmax - imax) ? flag : jmax - imax;
						imax = jmax;
						curpeak->sig = curdev / sqrt(curpeak->yerr * curpeak->yerr + curpeak->next->yerr * curpeak->next->yerr * w1 * w1 + curpeak->prev->yerr * curpeak->prev->yerr * w2 * w2);
					}
					else {
						curpeak->sig = 0;
					}
					jmax++;
				}

			}

			// Print all points in the spline with their significance (0 for convex sections, positive for concave sections) 
			f = fopen("spline.txt", "a+");
			for (curpeak = cpeaks->first; curpeak; curpeak = curpeak->next) {
				printf("\n+ t: %lf y: %lg sig: %lf", curpeak->t, curpeak->y, curpeak->sig);
				fprintf(f, "%lf %lg %lf\n", curpeak->t, curpeak->y, curpeak->sig);
			}
			fclose(f);


			// Find highest peak
			highestpeak = cpeaks->first->next;
			for (curpeak = cpeaks->first->next; curpeak->next; curpeak = curpeak->next) {
				if (curpeak->y > highestpeak->y) {
					highestpeak = curpeak;
				}
			}
			tv = highestpeak->t;
			// Find global minimum
			minimum = cpeaks->first->next;
			for (curpeak = cpeaks->first->next; curpeak->next; curpeak = curpeak->next) {
				if (curpeak->y < minimum->y) {
					minimum = curpeak;
				}
			}
			// Calculate prominence of highest peak with respect to global minimum of this dataset
			highestpeak->sig = (highestpeak->y - minimum->y) / sqrt(highestpeak->yerr * highestpeak->yerr + minimum->yerr * minimum->yerr);
			printf("\n- Highest peak\nt: %lf y: %lg sig: %lf", tv, highestpeak->y, highestpeak->sig);


			// Store maximal asymmetry for later use
			curpeak = cpeaks->last;
			int it = 0;
			for (it = 0; it < npc - 1 && 2 * tv - t[it] > curpeak->t; it++);
			for (it = it; it < npc - 1; it++) {
				t1 = 2 * tv - t[it];
				while (t1 < curpeak->t && curpeak->prev) {
					curpeak = curpeak->prev;
				}
				if (!curpeak->prev) break;
				curdev = (y[it] - curpeak->y - (curpeak->next->y - curpeak->y) / (curpeak->next->t - curpeak->t) * (t1 - curpeak->t)) / err[it];
				if (curdev > asydev) {
					asydev = curdev;
					tasy = t[it];
				}
			}
			printf("\nasy: %lf %lf", tasy, asydev);

			// Find all other peaks
			curpeak = cpeaks->first->next;
			while (curpeak && curpeak->sig <= 0) curpeak = curpeak->next; //Jump to next concave section
			while (curpeak) {
				//Find start and end of the concave section
				startsection = curpeak->prev;
				while (curpeak->next && curpeak->sig > 0) curpeak = curpeak->next;
				endsection = curpeak;
				// Find peak in this concave section
				sectionpeak = 0;
				curpeak = startsection->next;
				while (sectionpeak == 0 && curpeak != endsection) {
					if (curpeak->y > curpeak->prev->y && curpeak->y > curpeak->next->y) {
						sectionpeak = curpeak;
					}
					else {
						curpeak = curpeak->next;
					}
				}
				// Find minimum between this peak and highest peak
				if (sectionpeak != 0 && sectionpeak != highestpeak) {
					if (sectionpeak->t < highestpeak->t) {
						minimum = sectionpeak->next;
						curpeak = sectionpeak->next;
						while (curpeak != highestpeak) {
							if (curpeak->y < minimum->y) minimum = curpeak;
							curpeak = curpeak->next;
						}
					}
					else {
						minimum = sectionpeak->prev;
						curpeak = sectionpeak->prev;
						while (curpeak != highestpeak) {
							if (curpeak->y < minimum->y) minimum = curpeak;
							curpeak = curpeak->prev;
						}
					}
					// Calculate prominence of this peak with respect to minimum with highest peak
					sectionpeak->sig = (sectionpeak->y - minimum->y) / sqrt(sectionpeak->yerr * sectionpeak->yerr + minimum->yerr * minimum->yerr);
				}
				// Calculate alternative prominence with respect to ends of the concave section.
				maxdev = 0;
				if (startsection->prev == 0) startsection = startsection->next;
				if (endsection->next == 0) endsection = endsection->prev;
				curpeak = startsection->next;
				while (curpeak != endsection) {
					w1 = (curpeak->t - startsection->t) / (endsection->t - startsection->t);
					w2 = 1 - w1;
					curdev = curpeak->y - startsection->y - (endsection->y - startsection->y) * w1;
					curdev /= sqrt(curpeak->yerr * curpeak->yerr + endsection->yerr * endsection->yerr * w1 * w1 + startsection->yerr * startsection->yerr * w2 * w2);
					if (curdev > maxdev) {
						if (sectionpeak == 0) sectionpeak = curpeak;
						maxdev = curdev;
					}
					curpeak = curpeak->next;
				}
				// Choose higher prominence between the two calculations
				if (sectionpeak && maxdev > sectionpeak->sig) sectionpeak->sig = maxdev;
				if (sectionpeak == 0) {
					startsection->sig = endsection->sig = 0;
				}

				// Remove other points in the same concave section of this peak
				curpeak = startsection->next; 
				while (curpeak != endsection) {
					curpeak = curpeak->next;
					if (curpeak->prev != sectionpeak) cpeaks->remove(curpeak->prev);
				}
				/*
				curpeak = sectionpeak->prev;
				while (curpeak != startsection) {
					curpeak = curpeak->prev;
					cpeaks->remove(curpeak->next);
				}
				curpeak = sectionpeak->next;
				while (curpeak->next && curpeak != endsection) {
					curpeak = curpeak->next;
					cpeaks->remove(curpeak->prev);
				}*/
				curpeak = endsection->next;
				while (curpeak && curpeak->sig <= 0) curpeak = curpeak->next; //Jump to next concave section
			}



			//	
			//// Identify most representative point (peak) in concave sections and remove other points
			//	for (curpeak = cpeaks->first; curpeak; curpeak = curpeak->next) {
			//		// Every time we meet two consecutive points with sig>0, we save only one and remove the other
			//		if (curpeak->sig > 1.e-100 && curpeak->prev->sig > 1.e-100) {
			//			maxdev = curpeak->prev->y - curpeak->prev->prev->y - (curpeak->next->y - curpeak->prev->prev->y) / (curpeak->next->t - curpeak->prev->prev->t)*(curpeak->prev->t - curpeak->prev->prev->t);
			//			curdev = curpeak->y - curpeak->prev->prev->y - (curpeak->next->y - curpeak->prev->prev->y) / (curpeak->next->t - curpeak->prev->prev->t)*(curpeak->t - curpeak->prev->prev->t);
			//			// maxdev is the deviation of the previous point, curdev is for the current point
			//			flag = 0;
			//			if (curpeak->prev->y > curpeak->prev->prev->y && curpeak->y > curpeak->next->y) {
			//				flag = (curpeak->prev->y > curpeak->y) ? 1 : 2;
			//			}
			//			// flag is 0 for concave sections without a true peak, 1 if previous point is a peak, 2 if current point is a peak
			//			// Current point is removed if flag=1 or flag=0 and the previous point has higher deviation
			//			// Otherwise, previous point is removed
			//			if ((flag==0 && curpeak->prev->sig > curpeak->sig) || (flag==1)){
			//				curpeak = curpeak->prev;						
			//				cpeaks->remove(curpeak->next);
			//				curdev = maxdev;
			//			}
			//			else {
			//				cpeaks->remove(curpeak->prev);
			//			}
			//			// Significance is recalculated correspondingly
			//			curpeak->sig = curdev / curpeak->yerr;
			//		}
			//	}

			// Calculate uncertainty interval of peaks
			for (curpeak = cpeaks->first->next; curpeak->next; curpeak = curpeak->next) {

				if (curpeak->sig > 1.e-100) {

					if (curpeak->prev == cpeaks->first) {
						curpeak->tl = -1.e100; // This peak is at the beginning of light curve
					}
					else {
						//For high significance peak we find the first data point that would give a deviation lower than sigmathr and use this as tl.
						int i = curpeak->prev->i;
						curdev = 1.e100;
						while (curdev > sigmathr) {
							i++;
							w1 = (curpeak->t - t[i]) / (curpeak->next->t - t[i]);
							w2 = 1 - w1;
							curdev = curpeak->y - y[i] - (curpeak->next->y - y[i]) * w1;
							curdev /= sqrt(curpeak->yerr * curpeak->yerr + curpeak->next->yerr * curpeak->next->yerr * w1 * w1 + err[i] * err[i] * w2 * w2);
						}
						curpeak->tl = t[i - 1];
					}

					// Right side of the uncertainty interval. Analogous to left side.
					if (curpeak->next == cpeaks->last) {
						curpeak->tr = 1.e100;
					}
					else {
						int i = curpeak->next->i;
						curdev = 1.e100;
						while (curdev > sigmathr) {
							i--;
							w1 = (curpeak->t - t[i]) / (curpeak->prev->t - t[i]);
							w2 = 1 - w1;
							curdev = curpeak->y - y[i] - (curpeak->prev->y - y[i]) * w1;
							curdev /= sqrt(curpeak->yerr * curpeak->yerr + curpeak->prev->yerr * curpeak->prev->yerr * w1 * w1 + err[i] * err[i] * w2 * w2);
						}
						curpeak->tr = t[i + 1];
					}
				}
			}

			// Removing anything that is not a maximum
			for (curpeak = cpeaks->first; curpeak; curpeak = p) {
				p = curpeak->next;
				if (curpeak->sig <=peakthr*0.01) cpeaks->remove(curpeak);
			}

			// cpeaks now only contains a list of maxima with their significance and time interval
			printf("\nAll peaks in this dataset");
			for (p = cpeaks->first; p; p = p->next) {
				printf("\n/^\\ %lf %lf %lf %le %lf", p->t, p->tl, p->tr, p->y, p->sig);
			}

		}
		// At the end of this cycle, peaklist is an array with the lists of peaks for each dataset
		// including their significance and time interval.


	// Joining all peak lists

		printf("\n\n- Joining all peak lists\n");
		// All peaks from different datasets are put together in the same list

		cpeaks = new dataset;
		for (int ic = 0; ic < nfil; ic++) {
			for (p = peaklist[ic]->first; p; p = p->next) {
				cpeaks->clone(p);
			}
		}

		for (p = cpeaks->first; p; p = p->next) {
			printf("\n%lf %lf %lf %le %lf", p->t, p->tl, p->tr, p->y, p->sig);
		}

		// Ordering peaks by significance

		printf("\n\n- Ordering peaks by significance \n");

		p = cpeaks->first;
		while (p) {
			datapoint* q;
			pm = p->next;
			flag = 0;
			while (pm) {
				q = pm->next;
				if (pm->sig > p->sig) {
					pm->prev->next = pm->next;
					if (pm->next) {
						pm->next->prev = pm->prev;
					}
					else cpeaks->last = pm->prev;
					if (p->prev) {
						p->prev->next = pm;
					}
					else cpeaks->first = pm;
					pm->prev = p->prev;
					p->prev = pm;
					pm->next = p;
					p = pm;
				}
				pm = q;
			}
			p = p->next;
		}

		for (p = cpeaks->first; p; p = p->next) {
			printf("\n%lf %lf %lf %le %lf", p->t, p->tl, p->tr, p->y, p->sig);
		}

		// Reducing list to independent peaks (necessary if more than one dataset exists)

		printf("\n\n- Reducing peak lists\n");

		// Check for peak pairs falling in both uncertainty ranges or partial peaks

		for (pm = cpeaks->first; pm; pm = pm->next) {
			for (p = pm->next; p; p = p->next) {
				if ((p->t<pm->tr && p->t>pm->tl && pm->t<p->tr && pm->t>p->tl) || (p->tl < -1.e99 && pm->t < p->tr) || (p->tr > 1.e99 && pm->t > p->tl) || (pm->tl < -1.e99 && p->t < pm->tr) || (pm->tr > 1.e99 && p->t > pm->tl)) {
					w1 = 1 / (pm->tr - pm->tl);
					w2 = 1 / (p->tr - p->tl);
					//if (p->sig > pm->sig) pm->sig = p->sig; //pm->sig = (pm->sig*w1+p->sig*w2) / (w1 + w2);
					if ((p->tr - p->tl < pm->tr - pm->tl) || (p->tr > 1.e99 && pm->tr > 1.e99 && p->tl > pm->tl) || (p->tl < -1.e99 && pm->tl < -1.e99 && p->tr < pm->tr)) {
						pm->t = p->t;
						//pm->sig = p->sig;
					}
					pm->sig += p->sig;
					if (p->tl > pm->tl) pm->tl = p->tl;
					if (p->tr < pm->tr) pm->tr = p->tr;
					pmm = p->prev;
					cpeaks->remove(p);
					p = pmm;
				}
			}
		}

		for (p = cpeaks->first; p; p = p->next) {
			printf("\n%lf %lf %lf %le %lf", p->t, p->tl, p->tr, p->y, p->sig);
		}

		// Check for peak pairs where only one falls in the other's uncertainty range

		printf("\n\n- Reducing peak lists 2\n");

		for (pm = cpeaks->first; pm; pm = pm->next) {
			for (p = pm->next; p; p = p->next) {
				if ((p->t<pm->tr && p->t>pm->tl) || (pm->t<p->tr && pm->t>p->tl)) {
					w1 = 1 / (pm->tr - pm->tl);
					w2 = 1 / (p->tr - p->tl);
					pm->sig += p->sig;//if (p->sig > pm->sig) pm->sig = p->sig;//pm->sig= (pm->sig*w1 + p->sig*w2) / (w1 + w2);
					if (p->tl > pm->tl) pm->tl = p->tl;
					if (p->tr < pm->tr) pm->tr = p->tr;
					if (p->t < pm->tr && p->t>pm->tl) pm->t = p->t;
					pmm = p->prev;
					cpeaks->remove(p);
					p = pmm;
				}
			}
		}

		for (p = cpeaks->first; p; p = p->next) {
			printf("\n%lf %lf %lf %le %lf", p->t, p->tl, p->tr, p->y, p->sig);
		}


		newpeaks = cpeaks;

		// Delete used data
		for (int ic = 0; ic < nfil; ic++) {
			for (p = peaklist[ic]->first; p; p = pm) {
				pm = p->next;
				delete p;
			}
		}
		free(peaklist);
		free(tt);
		free(yy);
		free(eerr);
		free(nnp);

		// Removing peaks less significant than peakthr and beyond nobspeaks 
		// Note that we may decide to use more than 2 peaks in the observed light curve to check for more initial conditions.
		// nobspeaks is set to 2 at the beginning of the code, but can be changed if necessary

		p = newpeaks->first->next;
		for (int i = 1; i < nobspeaks && p && p->sig > peakthr; i++) {
			p = p->next;
		}
		while (p) {
			pm = p;
			p = p->next;
			newpeaks->remove(pm);
			delete pm;
		}

		printf("\nfint: %lf", fint);

		// If not enough relevant peaks, use largest asymmetry as second peak
		if (newpeaks->first->next == 0 && tasy > 1) {
			printf("\n\nUsing maximal asymmetry %lf", tasy);
			newpeaks->addpoint(0, tasy, asydev, 1);
			datapoint* p;
			p = newpeaks->first;
			newpeaks->first = p->next;
			newpeaks->first->prev = 0;
			newpeaks->first->next = p;
			p->prev = newpeaks->first;
			p->next = 0;
		}

	}
// At this point, newpeaks contains the list of peaks to be used to generate initial conditions ordered by significance
// We are ready for matching to the template library

// Reading old best models, if any
// This block is used by to include the best models from previous run in the set of initial conditions for the next run.
// Conventions
// First letter: P for single-source-single-lens, L for binary lens, B for binary source
// Second letter: S for static, X for parallax, O for orbital motion
// Three numbers for models sorted by increasing chi square
// .txt

	current_path(eventname);
	if(!exists("InitCond"))	create_directory("InitCond");
	current_path("InitCond");
	
	searchstring = regex(".*Init.*\\.txt");
	for (auto const& itr : directory_iterator(".")) {
		string curfile = (itr).path().filename().string();
		if (regex_match(curfile, searchstring)) {
			remove(curfile);
		}
	}


// Single lens - Single Source initial conditions

	current_path(eventname);

	dn = 0;
	if (nostatic) {
		filebest = regex("PX.*\\.txt");
		strcpy(fileinit, "InitCondPX.txt");
		nps = 6;
	}
	else {
		filebest = regex("PS.*\\.txt");
		strcpy(fileinit, "InitCondPS.txt");
		nps = 4;
	}
	tt = (double*)malloc(sizeof(double) * nps * maxoldmodels);
	if (exists(runstring + "\\Models")) {
		current_path(runstring + "\\Models");
		for (auto const& itr : directory_iterator(".")) {
			string curfile = (itr).path().filename().string();
			if (regex_match(curfile, filebest)) {
				f = fopen(curfile.c_str(), "r");
				for (int j = 0; j < nps; j++) {
					fscanf(f, "%le", &tt[dn * nps + j]);
				}
				fclose(f);
				dn++;
			}
		}
	}
	printf("\n- Writing initial conditions for fitting to %s\n\n",fileinit);
	current_path(eventname);
	current_path("InitCond");
	f = fopen(fileinit, "w");
	int nu0 = 3, ntE = 5, nrho = 4;
	// First we write the number of peaks used (only 1) and the number of initial conditions that are going to be generated
	if (nostatic) {
		fprintf(f, "%d %d\n", 1, nu0 * ntE * nrho*2 + dn);
	}
	else {
		fprintf(f, "%d %d\n", 1, nu0* ntE* nrho + dn);
	}
	// Then we write the characteristics of the peaks used
	p = newpeaks->first;
	fprintf(f, "%le %le %le %le %le\n", p->t, p->tl, p->tr, p->y, p->sig);
	// First we write the initial conditions from previous best models
	for (int i = 0; i < dn; i++) {
		for (int j = 0; j < nps; j++) {
			fprintf(f, "%le ", tt[i * nps + j]);
		}
		fprintf(f, "\n");
	}
	// Here we write the initial conditions by matching the newpeaks to the peaks recorded in the template library
	for (int iu = 0; iu < nu0; iu++) {
		for (int itE = 0; itE < ntE; itE++) {
			for (int ir = 0; ir < nrho; ir++) {
				//			{u0, tE, t0, Rs}
				if (nostatic) {
					fprintf(f, "%le %le %le %le 0.0 0.0\n", pow(10., -2. + iu), pow(10., -1. + itE), p->t, pow(10., -3. + 1. * ir));
					fprintf(f, "%le %le %le %le 0.0 0.0\n", -pow(10., -2. + iu), pow(10., -1. + itE), p->t, pow(10., -3. + 1. * ir));
				}
				else {
					fprintf(f, "%le %le %le %le\n", pow(10., -2. + iu), pow(10., -1. + itE), p->t, pow(10., -3. + 1. * ir));
				}
			}
		}
	}
	fclose(f);
	free(tt);
	current_path(eventname);

	if (!(nostatic)) {
		printf("\n- Writing initial conditions for fitting with Parallax to PreInitCondPX.txt\n\n");

		dn = 0;
		tt = (double*)malloc(sizeof(double) * 6 * maxoldmodels);
		filebest = regex("PX.*\\.txt");
		if (exists(runstring + "\\Models")) {
			current_path(runstring + "\\Models");
			for (auto const& itr : directory_iterator(".")) {
				if (dn >= maxoldmodels) break;
				string curfile = (itr).path().filename().string();
				if (regex_match(curfile, filebest)) {
					f = fopen(curfile.c_str(), "r");
					for (int j = 0; j < 6; j++) {
						fscanf(f, "%le", &tt[dn * 6 + j]);
					}
					fclose(f);
					dn++;
				}
			}
		}
		current_path(eventname);
		current_path("InitCond");

		f = fopen("PreInitCondPX.txt", "w");
		fprintf(f, "%d\n", dn);
		for (int i = 0; i < dn; i++) {
			for (int j = 0; j < 6; j++) {
				fprintf(f, "%le ", tt[i * 6 + j]);
			}
			fprintf(f, "\n");
		}
		fclose(f);
		free(tt);
	}

// Binary Source initial conditions

	current_path(eventname);

	dn = 0;
	if (nostatic) {
		filebest = regex("BO.*\\.txt");
		strcpy(fileinit, "InitCondBO.txt");
		nps = 10;
	}
	else {
		filebest = regex("BS.*\\.txt");
		strcpy(fileinit, "InitCondBS.txt");
		nps = 7;
	}
	tt = (double*)malloc(sizeof(double) * nps * maxoldmodels);
	if (exists(runstring + "\\Models")) {
		current_path(runstring + "\\Models");
		for (auto const& itr : directory_iterator(".")) {
			if (dn >= maxoldmodels) break;
			string curfile = (itr).path().filename().string();
			if (regex_match(curfile, filebest)) {
				f = fopen(curfile.c_str(), "r");
				for (int j = 0; j < nps; j++) {
					fscanf(f, "%le", &tt[dn * nps + j]);
				}
				fclose(f);
				dn++;
			}
		}
	}
	printf("\n- Writing initial conditions for fitting to %s\n\n",fileinit);

	current_path(eventname);
	current_path("InitCond");

	f = fopen(fileinit, "w");
	nu0 = 3, ntE = 5;
	int nFR = 3;
	// First we write the number of peaks used and the number of initial conditions that are going to be generated
	if (nostatic) {
		fprintf(f, "%d %d\n", newpeaks->length, nu0 * nu0 * ntE * nFR * (newpeaks->length * (newpeaks->length - 1) / 2)*2 + dn);
	}
	else {
		fprintf(f, "%d %d\n", newpeaks->length, nu0* nu0* ntE* nFR* (newpeaks->length* (newpeaks->length - 1) / 2) + dn);
	}
	// Then we write the characteristics of the peaks used
	for (p = newpeaks->first; p; p = p->next) {
		fprintf(f, "%le %le %le %le %le\n", p->t, p->tl, p->tr, p->y, p->sig);
	}
	// First we write the initial conditions from previous best models
	for (int i = 0; i < dn; i++) {
		for (int j = 0; j < nps; j++) {
			fprintf(f, "%le ", tt[i * nps + j]);
		}
		fprintf(f, "\n");
	}
	// Here we write the initial conditions by matching the newpeaks
	for (pl = newpeaks->first; pl->next; pl = pl->next) {
		for (pr = pl->next; pr; pr = pr->next) {
			for (int iu = 0; iu < nu0; iu++) {
				for (int iu2 = 0; iu2 < nu0; iu2++) {
					for (int itE = 0; itE < ntE; itE++) {
						for (int iFR = 0; iFR < nFR; iFR++) {
							if (nostatic) {
								//{u0, t0, log_tE, log_Rs, xi1, xi2, omega, inc, phi, log_qs}
								double u01 = pow(10., -2. +  iu);
								double u02 = pow(10., -2. +  iu2);
								double tE = pow(10., -1. + itE);
								double qs = pow(10., (-1. + iFR) / 4.0);
								fprintf(f, "%le %le %le %le %le %le %le %le %le %le\n", u01, pl->t, tE, 0.0001,(pr->t - pl->t)/ tE/(1+qs)*qs, (-u02+u01)/ (1 + qs) * qs, 0.000001, 0.0001, 0.00001, qs);
								fprintf(f, "%le %le %le %le %le %le %le %le %le %le\n", u01, pl->t, tE, 0.0001, (pr->t - pl->t) / tE / (1 + qs) * qs, (u02 + u01) / (1 + qs) * qs, 0.000001, 0.0001, 0.00001, qs);
							}
							else {
								//			{tE, fluxratio, u0_1, u0_2, t0_1, t0_2, rho}
								fprintf(f, "%le %le %le %le %le %le %le\n", pow(10., -1. + itE), pow(10., -1. + iFR), pow(10., -2. + 2 * iu), pow(10., -2. + 2 * iu2), pl->t, pr->t, 0.0001);
							}
						}
					}
				}
			}
		}
	}
	fclose(f);
	free(tt);
	current_path("..");

	if (!(nostatic)) {
		printf("\n- Writing initial conditions for fitting with Xallarap to PreInitCondBO.txt\n\n");

		dn = 0;
		tt = (double*)malloc(sizeof(double) * 10 * maxoldmodels);
		filebest = regex("BO.*\\.txt");
		if (exists(runstring + "\\Models")) {
			current_path(runstring + "\\Models");
			for (auto const& itr : directory_iterator(".")) {
				if (dn >= maxoldmodels) break;
				string curfile = (itr).path().filename().string();
				if (regex_match(curfile, filebest)) {
					f = fopen(curfile.c_str(), "r");
					for (int j = 0; j < 10; j++) {
						fscanf(f, "%le", &tt[dn * 10 + j]);
					}
					fclose(f);
					dn++;
				}
			}
		}

		current_path(eventname);
		current_path("InitCond");

		f = fopen("PreInitCondBO.txt", "w");
		fprintf(f, "%d\n", dn);
		for (int i = 0; i < dn; i++) {
			for (int j = 0; j < 10; j++) {
				fprintf(f, "%le ", tt[i * 10 + j]);
			}
			fprintf(f, "\n");
		}
		fclose(f);
		free(tt);
	}

		
// Binary lens initial conditions

	current_path(eventname);
	current_path("InitCond");

	dn = 0;
	if (noparallax) {
		filebest = regex("LO.*\\.txt");
		strcpy(fileinit, "PreInitCondLO.txt");
		nps = 12;
	}
	else {
		if (nostatic) {
			filebest = regex("LX.*\\.txt");
			strcpy(fileinit, "PreInitCondLX.txt");
			nps = 9;
		}
		else {
			filebest = regex("LS.*\\.txt");
			strcpy(fileinit, "PreInitCondLS.txt");
			nps = 7;
		}
	}

	current_path(eventname);

	dn=0;
	tt=(double *) malloc(sizeof(double)*nps*maxoldmodels);
	if (exists(runstring + "\\Models")) {
		current_path(runstring + "\\Models");
		for (auto const& itr : directory_iterator(".")) {
			if (dn >= maxoldmodels) break;
			string curfile = (itr).path().filename().string();
			if (regex_match(curfile, filebest)) {
				f = fopen(curfile.c_str(), "r");
				for (int j = 0; j < nps; j++) {
					fscanf(f, "%le", &tt[dn * nps + j]);
				}
				fclose(f);
				dn++;
			}
		}
	}
	printf("\n- Writing initial conditions for fitting to %s\n\n",fileinit);

	// Reading parameters for initial conditions from the template library
	current_path(exedir);
	current_path("..");
	current_path("data");
	f=fopen("TemplateLibrary.txt","r");
	fscanf(f,"%d",&np);
	printf("\nTemplates in library: %d",np);
	yy=(double *) malloc(sizeof(double)*np*7); // yy will contain all the information read from TemplateLibrary.txt
	for(int i=0;i<np*7;i++){
		fscanf(f,"%lf",&yy[i]);
	}
	fclose(f);

	current_path(eventname);
	current_path("InitCond");

	f=fopen(fileinit,"w");
	// First we write the number of peaks used and the number of initial conditions that are going to be generated
	if (nostatic) {
		fprintf(f, "%d %d\n", newpeaks->length, np * (newpeaks->length * (newpeaks->length - 1))*2 + dn);
	}
	else {
		fprintf(f, "%d %d\n", newpeaks->length, np * (newpeaks->length * (newpeaks->length - 1)) + dn);
	}
	// Then we write the characteristics of the peaks used
	for(p=newpeaks->first;p;p=p->next){
		fprintf(f,"%le %le %le %le %le\n",p->t,p->tl,p->tr,p->y,p->sig);
	}
	// First we write the initial conditions from previous best models
	for(int i=0;i<dn;i++){
		for(int j=0;j<nps;j++){
			fprintf(f,"%le ",tt[i*nps+j]);
		}
		fprintf(f,"\n");
	}
	// Here we write the initial conditions by matching the newpeaks to the peaks recorded in the template library
	for(int i=0;i<np;i++){
		for(pl=newpeaks->first;pl->next;pl=pl->next){
			for(pr=pl->next;pr;pr=pr->next){
				t1=(pl->t<pr->t)? pl->t : pr->t;
				t2=(pl->t<pr->t)? pr->t : pl->t;
				tE=(t2-t1)/(yy[i*7+6]-yy[i*7+5]);  // yy[i*7+6] and yy[i*7+5] are the peak times of the ith template
				t0=t2-tE*yy[i*7+6];
				for(int j=0;j<5;j++){
					fprintf(f,"%le ",yy[i*7+j]); // We use the s,q,u0,alpha,rho parameters from the template
				}
				fprintf(f,"%le %le",tE,t0); // and use tE and t0 from the time matching
				if (nostatic) {
					fprintf(f, " 0.0 0.0"); // parallax for nostatic
					if (noparallax) {
						fprintf(f, " 0.0 0.0 0.0001"); // starting parameters for orbital motion
					}
				}
				fprintf(f, "\n");

				if (nostatic) { // Reflected initial condition for nostatic
					yy[i * 7 + 2] = -yy[i * 7 + 2];
					yy[i * 7 + 3] = -yy[i * 7 + 3];
					for (int j = 0; j < 5; j++) {
						fprintf(f, "%le ", yy[i * 7 + j]); 
					}
					fprintf(f, "%le %le", tE, t0);
					fprintf(f, " 0.0 0.0"); // parallax for nostatic
					if (noparallax) {
						fprintf(f, " 0.0 0.0 0.0001"); // starting parameters for orbital motion
					}
					fprintf(f, "\n");
				}

				tE=(t1-t2)/(yy[i*7+6]-yy[i*7+5]); // We also include the time-reverse matching
				t0=t1-tE*yy[i*7+6];
				yy[i*7+2]=-yy[i*7+2];  //u0 and alpha are reversed
				yy[i*7+3]=yy[i*7+3]+M_PI;
				tE=-tE;
				for(int j=0;j<5;j++){
					fprintf(f,"%le ",yy[i*7+j]);
				}
				fprintf(f, "%le %le", tE, t0); // and use tE and t0 from the time matching
				if (nostatic) {
					fprintf(f, " 0.0 0.0"); // parallax for nostatic
					if (noparallax) {
						fprintf(f, " 0.0 0.0 0.0001"); // starting parameters for orbital motion
					}
				}
				fprintf(f, "\n");

				if (nostatic) { // Reflected initial condition for nostatic
					yy[i * 7 + 2] = -yy[i * 7 + 2];
					yy[i * 7 + 3] = -yy[i * 7 + 3];
					for (int j = 0; j < 5; j++) {
						fprintf(f, "%le ", yy[i * 7 + j]);
					}
					fprintf(f, "%le %le", tE, t0);
					fprintf(f, " 0.0 0.0"); // parallax for nostatic
					if (noparallax) {
						fprintf(f, " 0.0 0.0 0.0001"); // starting parameters for orbital motion
					}
					fprintf(f, "\n");
				}
			}
		}
	}
	fclose(f);
	free(tt);

	
// Reading old best parallaxmodels, if any

	current_path(eventname);

	if (!nostatic) {
		printf("\n- Writing initial conditions for fitting with Parallax to PreInitCondLX.txt\n\n");

		dn = 0;
		tt = (double*)malloc(sizeof(double) * 9 * maxoldmodels);
		filebest = regex("LX.*\\.txt");
		if (exists(runstring + "\\Models")) {
			current_path(runstring + "\\Models");
			for (auto const& itr : directory_iterator(".")) {
				if (dn >= maxoldmodels) break;
				string curfile = (itr).path().filename().string();
				if (regex_match(curfile, filebest)) {
					f = fopen(curfile.c_str(), "r");
					for (int j = 0; j < 9; j++) {
						fscanf(f, "%le", &tt[dn * 9 + j]);
					}
					fclose(f);
					dn++;
				}
			}
		}

		current_path(eventname);
		current_path("InitCond");

		f = fopen("PreInitCondLX.txt", "w");
		fprintf(f, "%d\n", dn);
		for (int i = 0; i < dn; i++) {
			for (int j = 0; j < 9; j++) {
				fprintf(f, "%le ", tt[i * 9 + j]);
			}
			fprintf(f, "\n");
		}
		fclose(f);
		free(tt);
	}
	
// Reading old best orbitalmodels, if any

	current_path(eventname);

	if (!noparallax) {
		printf("\n- Writing initial conditions for fitting with Orbital Motion to PreInitCondLO.txt\n\n");

		dn = 0;
		tt = (double*)malloc(sizeof(double) * 12 * maxoldmodels);
		filebest = regex("LO.*\\.txt");
		if (exists(runstring + "\\Models")) {
			current_path(runstring + "\\Models");
			for (auto const& itr : directory_iterator(".")) {
				if (dn >= maxoldmodels) break;
				string curfile = (itr).path().filename().string();
				if (regex_match(curfile, filebest)) {
					f = fopen(curfile.c_str(), "r");
					for (int j = 0; j < 12; j++) {
						fscanf(f, "%le", &tt[dn * 12 + j]);
					}
					fclose(f);
					dn++;
				}
			}
		}

		current_path(eventname);
		current_path("InitCond");

		f = fopen("PreInitCondLO.txt", "w");
		fprintf(f, "%d\n", dn);
		for (int i = 0; i < dn; i++) {
			for (int j = 0; j < 12; j++) {
				fprintf(f, "%le ", tt[i * 12 + j]);
			}
			fprintf(f, "\n");
		}
		fclose(f);
		free(tt);
	}


	printf("\n---- Done");
	//Sleep(5000l);


	free(yy);
	delete newpeaks;

    return 0;
}

////////////////////////////////
/////////////////////////////
///////////////////////////
// Methods in dataset and datapoint structures
///////////////////////////
////////////////////////////

dataset::dataset(){
	length=0;
	prev=next=0;
	first=last=0;
}

dataset::~dataset(){
	datapoint *p,*q;

	delete next;
	for(p=first;p;p=q){
		q=p->next;
		delete p;
	}
}

datapoint *dataset::addpoint(int i,double t, double y,double err){
	datapoint *p;
	p=new datapoint;
	p->i = i;
	p->t=t;
	p->y=y;
	p->yerr = err;
	if(length){
		if (p->t < first->t) {
			first->prev = p;
			p->next = first;
			p->prev = 0;
			first = p;
		}else{
			if (p->t >= last->t) {
				last->next = p;
				p->prev = last;
				p->next = 0;
				last = p;
			}		
			else {
				datapoint *scan = first;
				while (scan->next->t < t) scan = scan->next;
				scan->next->prev = p;
				p->next = scan->next;
				scan->next = p;
				p->prev = scan;
			}
		}
	}else{
		first=last=p;
		p->prev=0;
		p->next = 0;
	}
	p->tl=p->tr=0;
	p->sig=0;
	length++;
	return p;
}

void dataset::addpoint(double *t, double *y,double *err,int i,int n,int sign){
	datapoint *p;
	int flag,j;
	p=new datapoint;
	if(length){
		p->prev=last;
		last->next=p;
		last=p;
	}else{
		first=last=p;
		p->prev=0;
	}
	p->next=0;
	length++;

	p->t=t[i];
	p->y=y[i];
	p->yerr=err[i];
	p->sig=0;

	flag=0;
	j=i-1;
	while(flag<2){
		if((j<0)||(sign*(p->y-y[j])>sigmathr*(err[j]+p->yerr))){
			flag++;
		}else flag=0;
		j--;
	}
	if(j<0){
		p->tl=-1.e100;
	}else{
		p->tl=t[j+1];
	}

	flag=0;
	j=i+1;
//	double xx,tt;
	while(flag<2){
		//xx=(sign*(p->y-y[j])/(err[j]+p->yerr));
		//tt=t[j];
		if((j>=n) || (sign*(p->y-y[j])>sigmathr*(err[j]+p->yerr))){
			flag++;
		}else flag=0;
		j++;
	}
	if(j>=n){
		p->tr=1.e100;
	}else{
		p->tr=t[j-1];
	}


	//if(i<=1){
	//	p->terr=-1.e100;
	//}else{
	//	if(i>=n-2){
	//		p->terr=1.e100;
	//	}else{
	//		p->terr=(t[i+2]-t[i-2])/2;
	//	}
	//}
}

void dataset::remove(datapoint *p){
	if(first==p){
		first=p->next;
	}else{
		 p->prev->next=p->next;
	}
	if(last==p){
		last=p->prev;
	}else{
		p->next->prev=p->prev;
	}
	length--;
}

void dataset::clone(datapoint *p){
	datapoint *q;
	q=addpoint(p->i,p->t,p->y,p->yerr);
	q->sig=p->sig;
	q->tl=p->tl;
	q->tr=p->tr;
}