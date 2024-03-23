// Reader.cpp : main project file.
// This program formats data for a specific event for following work.
// for following work.

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

double tau=0.1; // conventional correlation time for consecutive points
int npmax=4000; // maximum number of points left after re-binning
int otherseasons = 1; // How to use other seasons
int renormalize = 1; // Re-normalize error bars
double thresholdoutliers = 10; // Threshold for removing outliers

struct datapoint{
	datapoint *prev,*next;
	double t,y,err,sig,basesig;
	int dn;
};

struct dataset{
	dataset *prev,*next;
	datapoint *first,*last;
	int length;
	char label[100];

	dataset();
	~dataset();
	void addpoint(double,double,double);
};

#define _computesig\
	p->sig=(p->y-p->prev->y);\
	pc=(p->t-p->prev->t)/tau;\
	p->sig*=p->sig/(p->err*p->err + p->prev->err*p->prev->err);\
	p->sig+=pc*pc;\
	p->sig *= p->basesig;



int main(int argc, char *argv[])
{
	char exedir[256]="";
	char eventname[512] = "";
	char filename[256] = "";
	char titstring[256]="",nostr[2],*undersc;
	char command[256],buffer[256];
	double value;
	double t,y,err,yr,errr,w1,w2;
	FILE *f;
	int ifile,flag,nps,normalized=0,satellite;
	double pc,residual,residual1,residual2,outlier,crosscheck,weight,minfac,maxlength;
	dataset *datalist=0,*curdataset,*pmaxdataset;
	datapoint *p,*pmax,*p1,*p2;
	
// Directory preliminaries. Reads event name from arguments.

	printf("******************************************\n");
	printf("*************      Reader      **********\n");
	printf("******************************************\n\n");
	printf("*** This program formats data for a specific event\n\n");


	if(argc>1){
		strcpy(eventname,argv[1]);
	}
	else {
		printf("\n\nEvent name? ");
		scanf("%s", eventname);
		//sprintf(eventname, "WDC10193");
	}

	printf("\n\n- Event: %s\n",eventname);

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
					npmax = (int) value;
				}
				if (strcmp(command, "tau") == 0) {
					tau = value;
				}
				if (strcmp(command, "otherseasons") == 0) {
					otherseasons = (int) value;
				}
				if (strcmp(command, "renormalize") == 0) {
					renormalize = (int) value;
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
	printf("\nRenormalize error bars: %s", (renormalize>0)? "yes" : "no");
	if (renormalize == 0) normalized = 1;


// Read data from files and create structured lists



	printf("\n\n- Reading data from files\n");

	curdataset = datalist = 0;
	current_path("Data");
	auto searchstring = regex(".*\\.dat");
	for (auto const& itr : directory_iterator(".")) {
		string curfile = itr.path().filename().string();
		if (regex_match(curfile, searchstring)) {
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
			while (fscanf(f, "%[^\n]", titstring) == 1) {
				fscanf(f, "%[\n]", nostr);
				sscanf(titstring, "%lf %lf %lf", &y, &err, &t);
				// Sanity checks on individual point
				if (y < 90 && err >0) {
					curdataset->addpoint(t, y, err);
				}
			}
			fclose(f);
			printf("\npoints: %d", curdataset->length);
		}
	}
	//auto itr = directory_iterator(".");
	//while (!itr._At_end()) {
	//	string curfile = (*itr).path().filename().string();
	//	if (regex_match(curfile, searchstring)) {
	//		strcpy(filename, curfile.c_str());
	//		if (!datalist) {
	//			curdataset = datalist = new dataset;
	//		}
	//		else {
	//			curdataset->next = new dataset;
	//			curdataset->next->prev = curdataset;
	//			curdataset = curdataset->next;
	//		}
	//		printf("\n%s", filename);
	//		strcpy(curdataset->label,filename);
	//		f = fopen(filename, "r");
	//		fscanf(f, "%[^\n]", titstring);
	//		fscanf(f, "%[\n]", nostr);
	//		while (fscanf(f, "%[^\n]", titstring) == 1) {
	//			fscanf(f, "%[\n]", nostr);
	//			sscanf(titstring, "%lf %lf %lf", &y, &err, &t);
	//			// Sanity checks on individual point
	//			if (y < 90 && err >0) {
	//				curdataset->addpoint(t, y, err);
	//			}
	//		}
	//		fclose(f);
	//		printf("\npoints: %d", curdataset->length);
	//		fclose(f);
	//	}
	//	itr++;
	//}
	current_path(eventname);

// Sort lists

	printf("\n\n- Time ordering of data\n");

	curdataset = datalist;
	while(curdataset){
		printf("\n%s", curdataset->label);
		p = curdataset->first;
		while (p) {
			pmax = p;
			p1 = pmax->next;
			while (p1) {
				if (p1->t == pmax->t) {
					p2 = p1->prev;
					if (p1->prev)	p1->prev->next = p1->next;
					if (p1->next)	p1->next->prev = p1->prev;
					delete p1;
					curdataset->length--;
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
				pc = pmax->t;
				pmax->t = p->t;
				p->t = pc;
				pc = pmax->y;
				pmax->y = p->y;
				p->y = pc;
				pc = pmax->err;
				pmax->err = p->err;
				p->err = pc;
			}
			p = p->next;
		}

		curdataset = curdataset->next;
		//if (curdataset->length <= 3) {
		//	printf(" discarded!");
		//	if (curdataset->prev) {
		//		dataset* predataset = curdataset->prev;
		//		curdataset = curdataset->next;
		//		delete predataset->next;
		//		predataset->next = curdataset;
		//		if (curdataset) curdataset->prev = predataset;
		//	}
		//	else {
		//		datalist = curdataset->next;
		//		if (datalist) datalist->prev = 0;
		//		delete curdataset;
		//		curdataset = datalist;
		//	}
		//}
		//else {
		// 			curdataset = curdataset->next;
		//}
		//pmax=curdataset->first->next;
		//while(pmax){
		//	flag=0;
		//	p=curdataset->first;
		//	while(p->t<=pmax->t && p!=pmax) p=p->next;
		//	if(p!=pmax){
		//		if(p->t!=pmax->t){
		//			pc=pmax->t;
		//			pmax->t=p->t;
		//			p->t=pc;
		//			pc=pmax->y;
		//			pmax->y=p->y;
		//			p->y=pc;
		//			pc=pmax->err;
		//			pmax->err=p->err;
		//			p->err=pc;
		//		}else{
		//			if(p->prev)	p->prev->next=p->next;
		//			if(p->next)	p->next->prev=p->prev;
		//			delete p;
		//			curdataset->length--;
		//		}
		//		flag=1;
		//	}
		//	if(!flag)	pmax=pmax->next;
		//}
	}


// Estimating error bars for each dataset
	
	double ft1, ft2;
	printf("\n\n- Assessing error bars\n");

	minfac=1.e100;
	maxlength=-minfac;
	for(curdataset=datalist;curdataset;curdataset=curdataset->next){
		residual=weight=0;
		for(p=curdataset->first->next;p;p=p->next){
			residual1=residual2 = 0;
			if(p->next){
				p1=p->next;
				if(p1->next){
					p2=p1->next;
					ft2 = (p->t - p1->t) / (p2->t - p1->t);
					ft1 = 1 - ft2;
					y=p1->y+((p2->t==p1->t)? 0 : (p2->y-p1->y)*ft2);
					pc=(p2->t-p->t)/tau;
					err=(p1->err*p1->err*ft1*ft1+p2->err*p2->err*ft2*ft2 +p->err*p->err)*exp(pc/tau);
					w1=1/exp(pc/tau);
					pc=(p->y-y);
					residual1= pc /sqrt(err);
					yr = y;
					errr = err;
				}
			}
			if(p->prev){
				p1=p->prev;
				if(p1->prev){
					p2=p1->prev;
					ft2 = (p->t - p1->t) / (p2->t - p1->t);
					ft1 = 1 - ft2;
					y=p1->y+((p2->t==p1->t)? 0 : (p2->y-p1->y)/(p2->t-p1->t)*(p->t-p1->t));
					pc=(p2->t-p->t)/tau;
					err=(p1->err*p1->err*ft1*ft1+p2->err*p2->err*ft2*ft2+p->err*p->err)*exp(-pc/tau);
					w2=1/exp(-pc/tau);
					pc=(p->y-y);
					residual2 = pc / sqrt(err);
				}
			}
			outlier = residual1 * residual1 + residual2 * residual2;
			crosscheck = (y - yr);
			crosscheck *= crosscheck;
			crosscheck /= (err + errr);
			if (residual1!=0 && residual2!=0  && crosscheck <9 && sqrt(outlier) >thresholdoutliers) {
				printf("\nOutlier found: %lf %lf %lf", p->t, p->y, sqrt(outlier));
				p1 = p ->prev;
				p2 = p->next;
				p1->next = p2;
				p2->prev = p1;
				delete p;
				curdataset->length--;
				p = p1;
			}
			else {
				residual += outlier;
				weight += w1 + w2;
			}
		}
		// 		residual*=0.5*curdataset->length/weight;
		residual=(residual+1.)/(weight+1.);
		pc=sqrt(residual);
		printf("\n%s", curdataset->label);
		printf("\nResidual: %le      Length: %d      Weight: %le       Normalization factor: %le",residual,curdataset->length,weight,pc);
		curdataset->first->sig=pc;
		if(curdataset->length>maxlength){
			minfac=pc;
			maxlength=curdataset->length;
		}
	}

	// Re-normalizing error bars for each dataset

	if (normalized < 0.5) {
		printf("\n\n- Re-normalizing error bars");
		for (curdataset = datalist; curdataset; curdataset = curdataset->next) {
			pc = curdataset->first->sig;// / minfac;
			for (p = curdataset->first->next; p; p = p->next) {
				p->err *= pc;
			}
			curdataset->first->sig = 0;
		}
	}


	// Eliminating short datasets

	printf("\n\n- Removing short datasets\n");

	pmaxdataset = datalist;
	nps = 0;
	while (pmaxdataset) {
		curdataset = pmaxdataset->next;
		if (pmaxdataset->length < 4) {
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

// Calculate peak season

	if (otherseasons > 0) {

		printf("\n\n- Calculate peak season\n");
		
		for (curdataset = datalist; curdataset; curdataset = curdataset->next) {
			datapoint** chunkfirst, ** chunklast;
			int nchunks, ns, imaxdev;
			double* devs;
			double sy, sy2, ss, maxdev;
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
				devs[ichunk] = sqrt((sy2 - sy * sy / ss) / ns);
				// Note chunk with maximal deviation (should be peak season)
				if (devs[ichunk] > maxdev) {
					maxdev = devs[ichunk];
					imaxdev = ichunk;
				}
				printf("\n%lf --- %lf : dev: %lf", chunkfirst[ichunk]->t, chunklast[ichunk]->t, devs[ichunk]);
			}

			if (otherseasons == 1) { // Diminish significance of seasons other than peak season
				printf("\nSeasons other than peak season considered less significant in re-binning");
				for (int ichunk = 0; ichunk < nchunks; ichunk++) {
					if (ichunk != imaxdev) {
						double fac = devs[ichunk] / maxdev;
						fac *= fac;
						for (p = chunkfirst[ichunk]; p != chunklast[ichunk]->next; p = p->next) {
							p->basesig = fac;
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

	for(curdataset=datalist;curdataset;curdataset=curdataset->next){
		for(p=curdataset->first->next;p;p=p->next){
			_computesig
		}
	}


// Rebin useless points, reducing the number of points to npmax

	printf("\n- Rebinning data down to %d\n",npmax);
	while(nps>npmax){
		pmax=datalist->first->next;
		pmaxdataset=datalist;
		ifile=flag=0;
		for(curdataset=datalist;curdataset;curdataset=curdataset->next){
			for(p=curdataset->first->next;p;p=p->next){
				if(p->sig<pmax->sig){
					pmax=p;
					pmaxdataset=curdataset;
					flag=ifile;
				}
			}
			ifile++;
		}
		p=pmax->prev;
//			printf("\n%d %d %lf\n%lf %le %le\n%lf %le %le",nps,flag,pmax->sig,p->t,p->y,p->err,pmax->t,pmax->y,pmax->err);
		w1=1/(p->err*p->err);
		w2=1/(pmax->err*pmax->err);
		p->t=(p->t*w1+pmax->t*w2)/(w1+w2);
		p->y=(p->y*w1+pmax->y*w2)/(w1+w2);
		p->err=1/sqrt(w1+w2);
		p->next=pmax->next;
//			printf("\n%lf %le %le",p->t,p->y,p->err);

		if(pmax->next){
			pmax->next->prev=p;
		}else{
			pmaxdataset->last=p;
		}
		pmaxdataset->length--;
		delete pmax;
		nps--;

		if(p->t>-1.e99 && p->y>-1.e99){
		}else{
			if(p->prev){
				p->prev->next=p->next;
			}else{
				datalist->first=p->next;
			}
			if(p->next){
				p->next->prev=p->prev;
			}else{
				datalist->last=p->prev;
			}
			delete p;
			pmaxdataset->length--;
			nps --;
		}


		if(pmaxdataset->length<4){
			if(pmaxdataset->prev){
				pmaxdataset->prev->next=pmaxdataset->next;
			}else{
				datalist=pmaxdataset->next;
			}
			if(pmaxdataset->next) pmaxdataset->next->prev=pmaxdataset->prev;
			nps-=pmaxdataset->length;
			pmaxdataset->next=0;
			delete pmaxdataset;
		}else{
			if(p->prev){
				_computesig
			}
			p=p->next;
			if(p){
				_computesig
			}
		}
	}

// Write rebinned datasets to file

	printf("\n- Writing final data to LCToFit.txt\n");

	remove("LCToFit.bkp");
	rename("LCToFit.txt","LCToFit.bkp");

	f=fopen("LCToFit.txt","w");
	ifile=0;
	fprintf(f,"%d\n",nps);
	for(curdataset=datalist;curdataset;curdataset=curdataset->next){
		printf("\n%d",curdataset->length);
		undersc=curdataset->label-5+strlen(curdataset->label);
		satellite=(*undersc<'A')? (*undersc)-'0' : 0; // satellite data should have names like ZOB1501241.dat where 1.dat distinguishes satellites data
		for(p=curdataset->first;p;p=p->next){
			y=pow(10.,-0.4*p->y);
			err=p->err*y*0.9210340371976184;
			fprintf(f,"%d %.10le %.10le %.10le %d\n",ifile,p->t,y,err,satellite);
		}
		ifile++;
	}
	fclose(f);

	f=fopen("FilterToData.txt","w");
	for(curdataset=datalist;curdataset;curdataset=curdataset->next){
		fprintf(f,"%s\n",curdataset->label);
	}
	fclose(f);

// Exiting

	delete datalist;

	//printf("\nHello!");
//	getchar();

    return 0;
}


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

void dataset::addpoint(double t, double y, double err){
	datapoint *p;
	p=new datapoint;
	p->t=t;
	p->y=y;
	p->err=err;
	if(length){
		p->prev=last;
		last->next=p;
		last=p;
	}else{
		first=last=p;
		p->prev=0;
	}
	p->next=0;
	p->sig=0;
	p->basesig = 1;
	length++;
}