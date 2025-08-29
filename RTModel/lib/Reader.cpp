// Reader.cpp : main project file.
// This program formats data for a specific event for following work.
// for following work.

#define _CRT_SECURE_NO_WARNINGS
#define _USE_MATH_DEFINES
#include <cstdio>
#include <cstring>
#include <cmath>
#include <cstdlib>
#include <regex>
#include <filesystem>
#include <list>

using namespace std;
using namespace std::filesystem;

#ifdef _WIN32
char systemslash = '\\';
#else
char systemslash = '/';
#endif

double tau = 1.0; // time-scale to trigger scatter assessment
int npmax = 4000; // maximum number of points left after re-binning
int otherseasons = 100; // How to use other seasons
int renormalize = 1; // Re-normalize error bars
double thresholdoutliers = 10; // Threshold for removing outliers

struct datapoint {
	datapoint* prev, * next;
	double t, y, err, sig, basesig, Dec, errDec, RA, errRA;
	int dn;
};

struct dataset {
	dataset* prev, * next;
	datapoint* first, * last;
	int length;
	char label[100];
	bool input_in_mags;

	dataset();
	~dataset();
	void addpoint(double, double, double, double, double, double, double);
	void deletepoint(datapoint*);
	void swappoints(datapoint*,datapoint*);
};

#define _computesig\
	if(p->err>0){\
		p->sig=(p->y-p->prev->y);\
		pc=(p->t-p->prev->t)/tau;\
		p->sig*=p->sig/(p->err*p->err + p->prev->err*p->prev->err);\
		p->sig+=pc*pc;\
		p->sig *= p->basesig;\
	}else p->sig = 1.e100;

int main(int argc, char* argv[])
{
	char exedir[256] = "";
	char eventname[512] = "";
	char filename[256] = "";
	char titstring[256] = "", nostr[2], * undersc;
	char command[256], buffer[256];
	double value;
	double t, y, err, Dec, errDec, RA, errRA, yr, ys, tr, ts, errr, errs, w1, w2;
	FILE* f;
	int ifile, flag, nps, normalized = 0, satellite;
	double pc, residual, residual1, residual2, residual3, outlier, crosscheck, weight, minfac, maxlength;
	dataset* datalist = 0, * curdataset, * pmaxdataset;
	datapoint* p, * pmax, * p1, * p2;

	// Directory preliminaries. Reads event name from arguments.

	setbuf(stdout, nullptr);
	printf("******************************************\n");
	printf("*************      Reader      **********\n");
	printf("******************************************\n\n");
	printf("*** This program formats data for a specific event\n\n");


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

	if (exists("LCToFit.txt")) {
		printf("\n\n- Data already processed");
		return 0;
	}

	// Reading Reader.ini and set parameters accordingly

	if (exists("ini")) {
		current_path("ini");
		f = fopen("Reader.ini", "r");
		if (f != 0) {
			printf("\n\n- Reading options in Reader.ini");
			while (!feof(f)) {
				int red = fscanf(f, "%s %s %lf", command, buffer, &value);
				if (red != 3) {
					command[0] = 0;
					//if (red != 0) {
					//	printf("\n\n!!! Bad command in Reader.ini");
					//	return -1;
					//};
				}
				if (strcmp(command, "binning") == 0) {
					npmax = (int)value;
				}
				if (strcmp(command, "tau") == 0) {
					tau = value;
				}
				if (strcmp(command, "otherseasons") == 0) {
					otherseasons = (int)value;
				}
				if (strcmp(command, "renormalize") == 0) {
					renormalize = (int)value;
				}
				if (strcmp(command, "thresholdoutliers") == 0) {
					thresholdoutliers = value;
				}
			}
			fclose(f);
		}
		else {
			printf("\n\n- Default options:");
		}
		current_path(eventname);
	}



	printf("\nUse of other seasons %d", otherseasons);
	printf("\nBinning down to %d", npmax);
	printf("\nThreshold for outliers removal: %lf", thresholdoutliers);
	printf("\nTimescale for scatter evaluation %lf", tau);
	printf("\nRenormalize error bars: %s", (renormalize > 0) ? "yes" : "no");
	if (renormalize == 0) normalized = 1;


	// Read data from files and create structured lists



	printf("\n\n- Reading data from files\n");

	curdataset = datalist = 0;
	current_path("Data");
	auto searchstring = regex(".*\\.dat");
	for (auto const& itr : directory_iterator(".")) {
		string curfile = itr.path().filename().string();
		if (regex_match(curfile, searchstring)) {
			int ncolumns;
			strcpy(filename, curfile.c_str());
			if (!datalist) {
				curdataset = datalist = new dataset;
			}
			else {
				curdataset->next = new dataset;
				curdataset->next->prev = curdataset;
				curdataset = curdataset->next;
			}
			printf("\n%s", filename);
			strcpy(curdataset->label, filename);
			f = fopen(filename, "r");
			fscanf(f, "%[^\n]", titstring);
			fscanf(f, "%[\n]", nostr);
			if (strstr(titstring, "Mag") != 0) curdataset->input_in_mags = true;
			while (fscanf(f, "%[^\n]", titstring) == 1) {
				fscanf(f, "%[\n]", nostr);
				ncolumns = sscanf(titstring, "%lg %lg %lg %lg %lg %lg %lg", &y, &err, &t, &Dec, &errDec, &RA, &errRA);
				if (y>-1.e100 && err>-1.e100 && t>-1.e100) {
					if (ncolumns < 7) {
						Dec = errDec = RA = errRA = -1;
					}
					// Sanity checks on individual point
					//if (y < 90 && err >0) {
					curdataset->addpoint(t, y, err, Dec, errDec, RA, errRA);
					//}
				}
			}
			fclose(f);
			printf("\npoints: %d", curdataset->length);
		}
	}

	current_path(eventname);

	// Sort lists

	printf("\n\n- Time ordering of data\n");

	curdataset = datalist;
	while (curdataset) {
		printf("\n%s", curdataset->label);
		p = curdataset->first;
		while (p) {
			pmax = p;
			p1 = pmax->next;
			while (p1) {
				if (p1->t == pmax->t) {
					p2 = p1->prev;
					curdataset->deletepoint(p1);
					p1 = p2;
				}
				else {
					if (p1->t < pmax->t) {
						pmax = p1;
					}				
				}
				p1 = p1->next;
			}
			if (p != pmax) {
				curdataset->swappoints(p, pmax);
				p = pmax;
			}
			p = p->next;
		}

		curdataset = curdataset->next;

	}

	// Estimating error bars for each dataset

	double ft1, ft2;
	printf("\n\n- Assessing error bars\n");

	minfac = 1.e100;
	maxlength = -minfac;
	for (curdataset = datalist; curdataset; curdataset = curdataset->next) {
		std::list<double> reslist;
		residual = weight = tr = ts = 0;
		if (curdataset->first->err > 0) {
			for (p = curdataset->first->next; p; p = p->next) {
				residual1 = residual2 = residual3 = 0;
				if (p->next) {
					p1 = p->next;
					if (p1->next) {
						p2 = p1->next;
						ft2 = (p->t - p1->t) / (p2->t - p1->t);
						ft1 = 1 - ft2;
						y = p1->y + ((p2->t == p1->t) ? 0 : (p2->y - p1->y) * ft2);
						pc = (p2->t - p->t);
						err = (p1->err * p1->err * ft1 * ft1 + p2->err * p2->err * ft2 * ft2 + p->err * p->err);// *exp(pc);
						//w1 = 1 / exp(pc);
						pc = (p->y - y);
						residual1 = pc / sqrt(err);
						yr = y;
						tr = p2->t;
						errr = err;
					}
				}
				if (p->prev) {
					p1 = p->prev;
					if (p1->prev) {
						p2 = p1->prev;
						ft2 = (p->t - p1->t) / (p2->t - p1->t);
						ft1 = 1 - ft2;
						y = p1->y + ((p2->t == p1->t) ? 0 : (p2->y - p1->y) / (p2->t - p1->t) * (p->t - p1->t));
						pc = (p2->t - p->t);
						err = (p1->err * p1->err * ft1 * ft1 + p2->err * p2->err * ft2 * ft2 + p->err * p->err);// *exp(-pc);
						ys = y;
						ts = p2->t;
						errs = err;
						//w2 = 1 / exp(-pc );
						pc = (p->y - y);
						residual2 = pc / sqrt(err);
					}
				}
				if (tr - ts < tau) {
					if (p->prev) {
						p1 = p->prev;
						if (p->next) {
							p2 = p->next;
							ft2 = (p->t - p1->t) / (p2->t - p1->t);
							ft1 = 1 - ft2;
							y = p1->y + ((p2->t == p1->t) ? 0 : (p2->y - p1->y) / (p2->t - p1->t) * (p->t - p1->t));
							pc = (p2->t - p->t);
							err = (p1->err * p1->err * ft1 * ft1 + p2->err * p2->err * ft2 * ft2 + p->err * p->err);// *exp(-pc);
							//w2 = 1 / exp(-pc);
							pc = (p->y - y);
							residual3 = pc / sqrt(err);
						}
					}
					outlier = residual1 + residual2 + 2 * residual3;
					crosscheck = (ys - yr);
					crosscheck *= crosscheck;
					crosscheck /= (err + errr);
					if (/*residual1*residual2>0 && */crosscheck < 9 && fabs(outlier) >4 * thresholdoutliers) {
						printf("\nOutlier found: %lf %lf %lf", p->t, p->y, outlier);
						p1 = p->prev;
						curdataset->deletepoint(p);
						/*p2 = p->next;
						p1->next = p2;
						p2->prev = p1;
						delete p;
						curdataset->length--;*/
						p = p1;
					}
					else {
						//residual += outlier;
						//weight += w1 + w2;
						reslist.push_back(outlier * outlier);
					}
				}
			}
			// 		residual*=0.5*curdataset->length/weight;
			if (!reslist.empty()) {
				reslist.sort();
				auto it = reslist.begin();
				int imed = (reslist.size() - 1) / 2;
				residual = *std::next(it, imed) * 0.186914;
			}
			else residual = 1;
			//residual = (residual + 1.) / (weight + 1.);
			pc = sqrt(residual);
			printf("\n%s", curdataset->label);
			//printf("\nResidual: %le      Length: %d      Weight: %le       Normalization factor: %le", residual, curdataset->length, weight, pc);
			printf("\nLength: %d      Res.list: %d     Residual: %lf      Normalization factor: %lf", curdataset->length, reslist.size(), residual, pc);
			curdataset->first->sig = pc;
			if (curdataset->length > maxlength) {
				minfac = pc;
				maxlength = curdataset->length;
			}
		}
	}

	// Re-normalizing error bars for each dataset

	if (normalized < 0.5) {
		printf("\n\n- Re-normalizing error bars");
		for (curdataset = datalist; curdataset; curdataset = curdataset->next) {
			if (curdataset->first->err > 0) {
				pc = curdataset->first->sig;// / minfac;
				for (p = curdataset->first; p; p = p->next) {
					p->err *= pc;
				}
				curdataset->first->sig = 0;
			}
		}
	}


	// Calculate peak season

	if (otherseasons != 0) {

		printf("\n\n- Calculate peak season\n");

		for (curdataset = datalist; curdataset; curdataset = curdataset->next) {
			datapoint** chunkfirst, ** chunklast;
			int nchunks, ns, imaxdev;
			double* devs;
			double sy, sy2, scury2, mean, dev, ss, maxdev;
			printf("\n%s", curdataset->label);
			// Divide dataset in chunks separated by gaps longer than 30 days
			chunkfirst = (datapoint**)malloc(sizeof(chunkfirst));
			chunklast = (datapoint**)malloc(sizeof(chunkfirst));
			chunkfirst[0] = curdataset->first;
			nchunks = 0;
			do {
				p = chunkfirst[nchunks];
				while (p && p->next && p->next->t - p->t < 30) p = p->next;
				chunklast[nchunks] = p;
				if (p->next) {
					nchunks++;
					chunkfirst = (datapoint**)realloc(chunkfirst, sizeof(chunkfirst) * (nchunks + 1));
					chunklast = (datapoint**)realloc(chunklast, sizeof(chunkfirst) * (nchunks + 1));
					chunkfirst[nchunks] = p->next;
				}
			} while (p->next);
			nchunks++;

			// Calculate standard deviation from flat light curve for each chunk
			devs = (double*)malloc(sizeof(double) * nchunks);
			imaxdev = 0;
			maxdev = -1;
			for (int ichunk = 0; ichunk < nchunks; ichunk++) {
				sy = ss = sy2 = 0;
				ns = 0;
				for (p = chunkfirst[ichunk]; p != chunklast[ichunk]->next; p = p->next) {
					sy += p->y / (p->err * p->err);
					sy2 += p->y * p->y / (p->err * p->err);
					ss += 1 / (p->err * p->err);
					ns++;
				}
				//				devs[ichunk] = sqrt((sy2 - sy * sy / ss) / ns);
				mean = sy / ss;
				scury2 = devs[ichunk] = 0;
				for (p = chunkfirst[ichunk]; p != chunklast[ichunk]->next; p = p->next) {
					dev = p->y - mean;
					if (dev>0) {
						dev /= p->err;
						scury2 += dev * dev;
					}
					else {
						if (scury2 > devs[ichunk]) devs[ichunk] = scury2;
						scury2 = 0;
					}
				}
				if (scury2 > devs[ichunk]) devs[ichunk] = scury2;
				// Note chunk with maximal deviation (should be peak season)
				if (devs[ichunk] > maxdev) {
					maxdev = devs[ichunk];
					imaxdev = ichunk;
				}
				printf("\n%lf --- %lf : dev: %lf", chunkfirst[ichunk]->t, chunklast[ichunk]->t, devs[ichunk]);
			}

			if (otherseasons > 0) { // Diminish significance of seasons other than peak season
				printf("\nSeasons other than peak season considered less significant in re-binning");
				for (int ichunk = 0; ichunk < nchunks; ichunk++) {
					if (ichunk != imaxdev) {
						double fac;
//						double fac = (devs[ichunk]+1.e-10) /(maxdev + 1.e-10);
//						fac *= fac;
						for (p = chunkfirst[ichunk]; p != chunklast[ichunk]->next; p = p->next) {
							fac = (p->t > chunklast[imaxdev]->t) ? (p->t > chunklast[imaxdev]->t) : (chunkfirst[imaxdev]->t - p->t);
							p->basesig = 1.0/otherseasons/fac;// fac;
						}
					}
				}
			}
			else { // Otherwise remove seasons other than peak season
				printf("\nSeasons other than peak season are removed");
				for (p = curdataset->first; p != chunkfirst[imaxdev]; p = p1) {
					p1 = p->next;
					delete p;
					curdataset->length--;
				}
				curdataset->first = chunkfirst[imaxdev];
				curdataset->first->prev = 0;
				for (p = chunklast[imaxdev]->next; p; p = p1) {
					p1 = p->next;
					delete p;
					curdataset->length--;
				}
				curdataset->last = chunklast[imaxdev];
				curdataset->last->next = 0;
			}

			free(devs);
			free(chunkfirst);
			free(chunklast);
		}
	}



	// Calculate significance of each point with respect to the previous one

	printf("\n- Calculate significance of each data point\n");

	for (curdataset = datalist; curdataset; curdataset = curdataset->next) {
		for (p = curdataset->first->next; p; p = p->next) {
			_computesig
		}
	}


	// Eliminating short datasets

	printf("\n\n- Removing short datasets\n");

	pmaxdataset = datalist;
	nps = 0;
	while (pmaxdataset) {
		curdataset = pmaxdataset->next;
		if (pmaxdataset->length < 2) {
			printf("\n- Discarding: %s", pmaxdataset->label);
			if (pmaxdataset->prev) {
				pmaxdataset->prev->next = pmaxdataset->next;
			}
			else {
				datalist = pmaxdataset->next;
			}
			if (pmaxdataset->next) pmaxdataset->next->prev = pmaxdataset->prev;
			pmaxdataset->next = 0;
			delete pmaxdataset;
		}
		else {
			nps += pmaxdataset->length;
		}
		pmaxdataset = curdataset;
	}


	// Rebin useless points, reducing the number of points to npmax

	printf("\n- Rebinning data down to %d\n", npmax);
	while (nps > npmax) {
		pmax = datalist->first->next;
		pmaxdataset = datalist;
		ifile = flag = 0;
		for (curdataset = datalist; curdataset; curdataset = curdataset->next) {
			for (p = curdataset->first->next; p; p = p->next) {
				if (p->sig < pmax->sig) {
					pmax = p;
					pmaxdataset = curdataset;
					flag = ifile;
				}
			}
			ifile++;
		}
		p = pmax->prev;
		//			printf("\n%d %d %lf\n%lf %le %le\n%lf %le %le",nps,flag,pmax->sig,p->t,p->y,p->err,pmax->t,pmax->y,pmax->err);
		w1 = 1 / (p->err * p->err);
		w2 = 1 / (pmax->err * pmax->err);
		p->t = (p->t * w1 + pmax->t * w2) / (w1 + w2);
		p->y = (p->y * w1 + pmax->y * w2) / (w1 + w2);
		p->err = 1 / sqrt(w1 + w2);
		if (p->errDec > 0) {
			w1 = 1 / (p->errDec * p->errDec);
			w2 = 1 / (pmax->errDec * pmax->errDec);
			p->Dec = (p->Dec * w1 + pmax->Dec * w2) / (w1 + w2);
			p->errDec = 1 / sqrt(w1 + w2);
			w1 = 1 / (p->errRA * p->errRA);
			w2 = 1 / (pmax->errRA * pmax->errRA);
			p->RA = (p->RA * w1 + pmax->RA * w2) / (w1 + w2);
			p->errRA = 1 / sqrt(w1 + w2);
		}
		pmaxdataset->deletepoint(pmax);
		nps--;

		if (p->t > -1.e99 && p->y > -1.e99) {
		}
		else {
			pmaxdataset->deletepoint(p);
			nps--;
		}


		if (pmaxdataset->length < 2) {
			if (pmaxdataset->prev) {
				pmaxdataset->prev->next = pmaxdataset->next;
			}
			else {
				datalist = pmaxdataset->next;
			}
			if (pmaxdataset->next) pmaxdataset->next->prev = pmaxdataset->prev;
			nps -= pmaxdataset->length;
			pmaxdataset->next = 0;
			delete pmaxdataset;
		}
		else {
			if (p->prev) {
				_computesig
			}
			p = p->next;
			if (p) {
				_computesig
			}
		}
	}



	// Write rebinned datasets to file

	printf("\n- Writing final data to LCToFit.txt\n");

	remove("LCToFit.bkp");
	rename("LCToFit.txt", "LCToFit.bkp");

	f = fopen("LCToFit.txt", "w");
	ifile = 0;
	fprintf(f, "%d\n", nps);
	for (curdataset = datalist; curdataset; curdataset = curdataset->next) {
		printf("\n%d", curdataset->length);
		undersc = curdataset->label - 5 + strlen(curdataset->label);
		satellite = (*undersc <= '9' && *undersc > '0') ? (*undersc) - '0' : 0; // satellite data should have names like ZOB1501241.dat where 1.dat distinguishes satellites data
		for (p = curdataset->first; p; p = p->next) {
			if (curdataset->input_in_mags) {
				y = pow(10., -0.4 * p->y);
				err = (p->err > 0) ? p->err * y * 0.9210340371976184 : p->err;
			}
			else {
				y = p->y;
				err = p->err;
			}
			fprintf(f, "%d %.10le %.10le %.10le %d %.10le %.10le %.10le %.10le\n", ifile, p->t, y, err, satellite, p->Dec, p->errDec,p->RA, p->errRA);
		}
		ifile++;
	}
	fclose(f);

	f = fopen("FilterToData.txt", "w");
	for (curdataset = datalist; curdataset; curdataset = curdataset->next) {
		fprintf(f, "%s\n", curdataset->label);
	}
	fclose(f);

	// Exiting

	delete datalist;

	//printf("\nHello!");
//	getchar();

	return 0;
}


dataset::dataset() {
	length = 0;
	prev = next = 0;
	first = last = 0;
	input_in_mags = false;
}

dataset::~dataset() {
	datapoint* p, * q;

	delete next;
	for (p = first; p; p = q) {
		q = p->next;
		delete p;
	}
}

void dataset::addpoint(double t, double y, double err, double Dec, double errDec, double RA, double errRA) {
	datapoint* p;
	p = new datapoint;
	p->t = t;
	p->y = y;
	p->err = err;
	p->Dec = Dec;
	p->errDec = errDec;
	p->RA = RA;
	p->errRA = errRA;
	if (length) {
		p->prev = last;
		last->next = p;
		last = p;
	}
	else {
		first = last = p;
		p->prev = 0;
	}
	p->next = 0;
	p->sig = 0;
	p->basesig = 1;
	length++;
}

void dataset::deletepoint(datapoint* pdel) {
	datapoint* p;
	if (pdel == first) {
		first = pdel->next;
		if (first) {
			first->prev = 0;
		}
		else {
			last = 0;
		}
	}
	else {
		if (pdel == last) {
			last = pdel->prev;
			last->next = 0;
		}
		else {
			pdel->prev->next = pdel->next;
			pdel->next->prev = pdel->prev;
		}
	}
	delete pdel;
	length--;
}

void dataset::swappoints(datapoint* p1, datapoint* p2) {
	datapoint* p1p, * p1n, * p2p, * p2n;
	if (first == p1) {
		first = p2;
	}
	else {
		if (first == p2) {
			first = p1;
		}
	}

	if (last == p1) {
		last = p2;
	}
	else {
		if (last == p2) {
			last = p1;
		}
	}
	p1p = p1->prev;
	p1n = p1->next;
	p2p = p2->prev;
	p2n = p2->next;
	if (p1n == p2) {
		if(p1p) p1p->next = p2;
		p2->prev = p1p;
		p2->next = p1;
		p1->prev = p2;
		p1->next = p2n;
		if(p2n) p2n->prev = p1;
	}
	else {
		if (p2n == p1) {
			if(p2p) p2p->next = p1;
			p1->prev = p2p;
			p1->next = p2;
			p2->prev = p1;
			p2->next = p1n;
			if(p1n) p1n->prev = p2;
		}
		else {
			if (p1p) p1p->next = p2;
			p2->prev = p1p;
			p2->next = p1n;
			if (p1n) p1n->prev = p2;
			if (p2p) p2p->next = p1;
			p1->prev = p2p;
			p1->next = p2n;
			if (p2n) p2n->prev = p1;
		}
	}

}
