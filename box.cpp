/*Run Photon Maps with uniform distribution*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <vector>
#include <limits>

#ifndef M_PI
#define M_PI 3.1415926535
#endif

typedef enum { X, Y, Z, MX, MY, MZ, UND
} normal_t;

void swap(double *x,double *y)
{
   double temp;
   temp = *x;
   *x = *y;
   *y = temp;
}

double selection_sort(std::vector< std::vector<float> > &dists, int k)
{
  int i, j, min;
  if (k>dists.size()) {
	fprintf (stderr, "Not enough photon impacts.");
	exit (-1);
  }
  for (i = 0; i < k ; i++)
  {
    min = i;
    for (j = i+1; j < dists.size(); j++)
    {
       if (dists[j][3] < dists[min][3])
       {
          min = j;
       }
    }
    //swap
    std::vector<float> tmp;
    tmp=dists[min];
    dists[min]=dists[i];
    dists[i]=tmp;
  }

return dists[k-1][3];
}


/* cornell box:
Coordinates: X=depth,+=advance, back of box; Y=left(-), right(+);  Z=vertical (+:up)
Based on a [-1,+1] box, axis aligned, bottom center on 0,0,0.
Back is at
Lamp: 0.4 side square, at height 1.99, centred
Small box, axis aligned. 0.4 side, centred at (0, -0.3, 0.3)
Point I(P) at (0,0,0)

*/

void emit (float dir [3])
//uniform in xy plane, project to unit sphere in +z
{
	float distancesq;
	do {
		for (int j=0;j<2;j++) {
			//positions[i][j]=(double)random()/RAND_MAX*2*radius;
			dir[j]=(double)rand()/RAND_MAX*2-1;
		}
		distancesq= (dir[0]*dir[0] + dir[1]*dir[1]);
	} while (distancesq>1);
	dir[2]=sqrtf (1-distancesq);
}
void normaldown (float dir[3]) {
dir[0]=0;
dir[1]=0;
dir[2]=-1;
}
void normal (float dir[3], normal_t nt)
{
	float tmp[3];
	switch(nt) {
		case Z:
			emit(dir);
			break;
		case MZ:
			emit(dir);
			dir[2]=-dir[2];
			break;
		case X:
			emit(tmp);
			dir[0]=tmp[2];
			dir[1]=tmp[1];
			dir[2]=tmp[0];
			break;
		case MX:
			emit(tmp);
			dir[0]=-tmp[2];
			dir[1]=tmp[1];
			dir[2]=tmp[0];
			break;
		case Y:
			emit(tmp);
			dir[0]=tmp[0];
			dir[1]=tmp[2];
			dir[2]=tmp[1];
			break;
		case MY:
			emit(tmp);
			dir[0]=-tmp[0];
			dir[1]=-tmp[2];
			dir[2]=tmp[1];
			break;
	}
}
int numberPhotons = 120000;
int strata = 120;
int count = 1;
int ring = 1; //ringMax = 15
int ringMax = 15;
int stratumSize = numberPhotons/strata;
int ringSize = stratumSize*ring;
double thetaBasis = 2*M_PI/ring;
int stratum = 1;

//stratified  final version

void emitStra (float dir [3])
//uniform in xy plane, project to unit sphere in +z
{
	float distancesq;
    double thetaUpper = thetaBasis*stratum;
    double thetaLower = thetaBasis*(stratum-1);
    double theta = (thetaUpper-thetaLower)*((double)rand()/RAND_MAX)+thetaLower;
    double ringSumO = (double)ring*(ring+1)/(double)2;
    double ringSumI = (double)(ring-1)*(ring)/(double)2;
    double rOuter = (sqrtf (ringSumO/(double)strata));
    double rInner = (sqrtf ((ringSumI)/(double)strata));
    double radius = (rOuter*rOuter-rInner*rInner)*((double)rand()/RAND_MAX)+rInner*rInner;
    dir[0]=sqrtf (radius)*cos(theta);
    dir[1]=sqrtf (radius)*sin(theta),
	distancesq= (dir[0]*dir[0] + dir[1]*dir[1]);
	dir[2]=sqrtf (1-distancesq);
    if(count == stratumSize) {
          stratum ++;
          count = 0;
    }
    if(stratum > ring){
          ring++;
          thetaBasis = 2*M_PI/ring;
          count = 0;
          stratum = 1;
    }
    if(ring > ringMax){
        ring = 1;
        thetaBasis = 2*M_PI/ring;
    }
    count++;
}
void normalStra (float dir[3], normal_t nt)
{
	float tmp[3];
	switch(nt) {
		case Z:
			emitStra(dir);
			break;
		case MZ:
			emitStra(dir);
			dir[2]=-dir[2];
			break;
		case X:
			emitStra(tmp);
			dir[0]=tmp[2];
			dir[1]=tmp[1];
			dir[2]=tmp[0];
			break;
		case MX:
			emitStra(tmp);
			dir[0]=-tmp[2];
			dir[1]=tmp[1];
			dir[2]=tmp[0];
			break;
		case Y:
			emitStra(tmp);
			dir[0]=tmp[0];
			dir[1]=tmp[2];
			dir[2]=tmp[1];
			break;
		case MY:
			emitStra(tmp);
			dir[0]=-tmp[0];
			dir[1]=-tmp[2];
			dir[2]=tmp[1];
			break;
        case UND:
            break;
	}
}
void copy (float *a, float *b) // [3]
{
	for (int i=0;i<3;i++)
		a[i]=b[i];
}

std::vector< std::vector<float> > hits; // x, y, z, distance to P=0,0,0, power

float Bt (float p[3], float d[3], float res[3], normal_t& n) {
	float l;
	l=(0-p[2]) / d[2];
	if (l<0)
		return l;
	res[0]=p[0]+l*d[0];
	res[1]=p[1]+l*d[1];
	res[2]=0;
	if (res[0]<-1 || res[0]>1 || res[1]<-1 || res[1]>1)
		return -1;
	n=Z;
	return l;
}
//C=cube
float BtC (float p[3], float d[3], float res[3], normal_t& n) {
	float l;
	l=(0-p[2]) / d[2];
	if (l<0)
		return l;
	res[0]=p[0]+l*d[0];
	res[1]=p[1]+l*d[1];
	res[2]=0.1;
	if (res[0]<-0.2 || res[0]>0.2 || res[1]<-0.5 || res[1]>-0.1)
		return -1;
	n=Z;
	return l;
}

float T (float p[3], float d[3], float res[3], normal_t& n) {
	float l;
	l=(2-p[2]) / d[2];
	res[0]=p[0]+l*d[0];
	res[1]=p[1]+l*d[1];
	if (res[0]<-1 || res[0]>1 || res[1]<-1 || res[1]>1)
		return -1;
	res[2]=2;
	n=MZ;
	return l;
}

float TC (float p[3], float d[3], float res[3], normal_t& n) {
	float l;
	l=(2-p[2]) / d[2];
	res[0]=p[0]+l*d[0];
	res[1]=p[1]+l*d[1];
	if (res[0]<-0.2 || res[0]>0.2 || res[1]<-0.5 || res[1]>-0.1)
		return -1;
	res[2]=0.5;
	n=MZ;
	return l;
}
//Y=1
float R (float p[3], float d[3], float res[3], normal_t& n) {
	float l;
	l=(1-p[1]) / d[1];
	res[0]=p[0]+l*d[0];
	res[1]=1;
	res[2]=p[2]+l*d[2];
	if (res[0]<-1 || res[0]>1 || res[2]<0 || res[2]>2)
		return -1;
	n=MY;
	return l;
}

float RC (float p[3], float d[3], float res[3], normal_t& n) {
	float l;
	l=(1-p[1]) / d[1];
	res[0]=p[0]+l*d[0];
	res[1]=-0.1;
	res[2]=p[2]+l*d[2];
	if (res[0]<-0.2 || res[0]>0.2 || res[2]<0.1 || res[2]>0.5)
		return -1;
	n=MY;
	return l;
}
float L(float p[3], float d[3], float res[3], normal_t& n) {
	float l;
	l=(-1-p[1]) / d[1];
	res[0]=p[0]+l*d[0];
	res[1]=-1;
	res[2]=p[2]+l*d[2];
	if (res[0]<-1 || res[0]>1 || res[2]<0 || res[2]>2)
		return -1;
	n=Y;
	return l;
}
float LC(float p[3], float d[3], float res[3], normal_t& n) {
	float l;
	l=(-1-p[1]) / d[1];
	res[0]=p[0]+l*d[0];
	res[1]=-0.5;
	res[2]=p[2]+l*d[2];
	if (res[0]<-0.2 || res[0]>0.2  || res[2]<0.1 || res[2]>0.5)
		return -1;
	n=Y;
	return l;
}


float Bk (float p[3], float d[3], float res[3], normal_t& n) {
	float l;
	l=(1-p[0]) / d[0];
	res[0]=1;
	res[1]=p[1]+l*d[1];
	res[2]=p[2]+l*d[2];
	if (res[2]<-0 || res[2]>2 || res[1]<-1 || res[1]>1)
		return -1;
	n=MX;
	return l;
}

float BkC (float p[3], float d[3], float res[3], normal_t& n) {
	float l;
	l=(1-p[0]) / d[0];
	res[0]=0.2;
	res[1]=p[1]+l*d[1];
	res[2]=p[2]+l*d[2];
	if (res[2]<-0.5 || res[2]>-0.1 || res[1]<-0.2 || res[1]>0.2)
		return -1;
	n=MX;
	return l;
}

//front box
float FC (float p[3], float d[3], float res[3], normal_t& n) {
	float l;
	l=(1-p[0]) / d[0];
	res[0]=-0.2;
	res[1]=p[1]+l*d[1];
	res[2]=p[2]+l*d[2];
	if (res[2]<-0.5 || res[2]>-0.1 || res[1]<-0.2 || res[1]>0.2)
		return -1;
	n=MX;
	return l;
}

void render (int n, int square, int disk) {
	//use russian roulette, p=0.5
hits.clear();
const float reflectivity=0.5;
int decCount = 0;
for (int i=0; i< n; i++) {
	float p[3];
	int j;
	if(square==0){
	 for (j=0;j<2;j++)
		 //Lamp: 0.4 side square, at height 1.99, centred
		 p[j]=(double)rand()/RAND_MAX*0.4-0.2;
	}
	else{
//   Stratified assuming 10X10 Lamp and 100 000 photons.
      double stratDia = 0.4/(double)10;
      double horizSt = (i)%10;
      p[0]=((double)rand()/RAND_MAX*stratDia-stratDia/2.0)-(double)(0.2-stratDia/2.0-(horizSt*stratDia));
      p[1]=((double)rand()/RAND_MAX*stratDia-stratDia/2.0)-(double)(0.2-stratDia/2.0-(decCount*stratDia));
      if(horizSt == 9){
	      decCount++;
	    	}
      if(decCount == 10){
          decCount = 0;
            }
	}
	p[2]=1.99;
	float flux=1.0f;
	float dir1[3];
	if(disk==1){
	   normalStra (dir1, MZ);
	}
	else{
       normal (dir1, MZ);
	}
	//intersect with scene
	float res[3];
	float minl=std::numeric_limits<float>::infinity();
	normal_t nt=UND;
	float l;
	normal_t lnt;
	float lres[3];
	bool absorbed=false;
	do {
	minl=std::numeric_limits<float>::infinity();
	//enclosure: bottom, top, left, right back
	l=Bt(p, dir1, res, nt);
	if (l>0) {//in
		minl=l;
		lnt=nt;
		copy (lres, res);
	}
	l=T(p, dir1, res, nt);
	if (l>0) {//in
		minl=std::min (l, minl);
		if (minl == l)
		{
		lnt=nt;
		copy (lres, res);
		}
	}
	l=L(p, dir1, res, nt);
	if (l>0) {//in
		minl=std::min (l, minl);
		if (minl == l)
		{
		lnt=nt;
		copy (lres, res);
		}
	}
	l=R(p, dir1, res, nt);
	if (l>0) {//in
		minl=std::min (l, minl);
		if (minl == l)
		{
		lnt=nt;
		copy (lres, res);
		}
	}
	l=Bk(p, dir1, res, nt);
	if (l>0) {//in
		minl=std::min (l, minl);
		if (minl == l)
		{
		lnt=nt;
		copy (lres, res);
		}
	}
//small cube
	l=BtC(p, dir1, res, nt);
	if (l>0) {//in
		minl=std::min (l, minl);
		if (minl == l)
		{
		lnt=nt;
		copy (lres, res);
		}
	}
	l=TC(p, dir1, res, nt);
	if (l>0) {//in
		minl=std::min (l, minl);
		if (minl == l)
		{
		lnt=nt;
		copy (lres, res);
		}
	}
	l=LC(p, dir1, res, nt);
	if (l>0) {//in
		minl=std::min (l, minl);
		if (minl == l)
		{
		lnt=nt;
		copy (lres, res);
		}
	}
	l=RC(p, dir1, res, nt);
	if (l>0) {//in
		minl=std::min (l, minl);
		if (minl == l)
		{
		lnt=nt;
		copy (lres, res);
		}
	}
	l=BkC(p, dir1, res, nt);
	if (l>0) {//in
		minl=std::min (l, minl);
		if (minl == l)
		{
		lnt=nt;
		copy (lres, res);
		}
	}
	l=FC(p, dir1, res, nt);
	if (l>0) {//in
		minl=std::min (l, minl);
		if (minl == l)
		{
		lnt=nt;
		copy (lres, res);
		}
	}
	if (minl<std::numeric_limits<float>::infinity()) {
		//store
		std::vector<float> h;
		for (int m=0;m<3;m++)
			h.push_back(lres[m]);
		h.push_back(sqrtf(h[0]*h[0]+h[1]*h[1]+h[2]*h[2])); // distance to P=0,0,0
		h.push_back(flux);
		hits.push_back(h);
		float russian=(double)rand()/RAND_MAX;
		if (russian < reflectivity) { 	//reflect
			for (int m=0;m<3;m++)
				p[m]=h[m];
			normal (dir1, lnt);

		} else absorbed=true;
	} else // leaves the scene
		absorbed=true;
	} while (!absorbed);
} // for n
} // render

double cone(double ri, double rk){

double result;
double first = 1 - (ri/rk);
double second = 1-(2.0/3.0);
result = first/second;

return result;
}
double epanechnikov (double ri, double rk){

double result;
double first = 1- (ri/rk)*(ri/rk);
result = 2 * first;
return result;
}
double silverman (double ri, double rk){

double result;
double first = 1- (ri/rk)*(ri/rk);
result = 3 * first*first;
return result;
}
double gaussian (double ri, double rk){

double result;
double alfa = 1.728;
double beta = 1.953;
double first = -beta*(ri*ri/(2*rk*rk));
double second = 1-std::exp(first);
double third = 1-std::exp(-beta);
result = alfa*(1 - second/third);

return result;
}
int main () {

srand (time (NULL));

char strategy [20];
char kernel [20];
int square;
int disk;
int n=120000;
int k;
int iterations;
printf("Cornell Box in stratified Photon Mapping\n");
printf("Number of photons is 120 000.\n");
printf("Enter a value  for k:\n");
scanf("%d",&k);
printf("Enter the number of iterations:\n");
scanf("%d",&iterations);
printf("Select a strategy ->>>> SQUARE|DISK|COMBINED|UNSTRATIFIED:\n");
scanf("%S",strategy);
printf("Select a kernel ->>>> EPA|CONE|GAUSS|SILV|NONE:\n");
scanf("%S",kernel);
switch (strategy[0]){
 case 'S':
    square=1;
    disk=0;
    break;
 case 'D':
     square=0;
    disk=1;
    break;
 case 'C':
    square=1;
    disk=1;
    break;
 case 'U':
     square=0;
     disk=0;
     break;
}
double resPM=0;
double resPMF=0;
double res2PM=0;
double res2PMF=0;
double power=1.0/n;
int iter;
for (iter=1;iter<=iterations;iter++){
    double powerdensityF;
    double powerdensity;
	printf ("iteration %d\n", iter);
	render(n,square,disk);
	double r=selection_sort (hits, k);
	double areak=M_PI*r*r;
    float acc = 0;
    switch (kernel[0]){
        case 'N':
           //Constant:
	       for (int i=0;i<k;i++)
		      acc+=hits[i][4];
           powerdensityF = acc * power / areak;
	       acc+=hits[k][4];
	       powerdensity = acc * power / areak;
           break;
        case 'C':
           //Cone slope s=1:
           for (int i=0;i<k;i++)
              acc+=hits[i][4]*cone(hits[i][3],r);
           powerdensityF = acc * power / areak;
	       acc+=hits[k][4]*cone(r,r);
	       powerdensity = acc * power / areak;
           break;
        case 'S':
           //Silverman:
           for (int i=0;i<k;i++)
              acc+=hits[i][4]*silverman(hits[i][3],r);
           powerdensityF = acc * power / areak;
	       acc+=hits[k][4]*silverman(r,r);
	       powerdensity = acc * power / areak;
           break;
        case 'E':
          //Epanechnikov:
           for (int i=0;i<k;i++)
             acc+=hits[i][4]*epanechnikov(hits[i][3],r);
           powerdensityF = acc * power / areak;
	       acc+=hits[k][4]*epanechnikov(r,r);
	       powerdensity = acc * power / areak;
           break;
        case 'G':
          //Gaussian:
           for (int i=0;i<k;i++)
             acc+=hits[i][4]*gaussian(hits[i][3],r);
           powerdensityF = acc * power / areak;
	       acc+=hits[k][4]*cone(r,r);
	       powerdensity = acc * power / areak;
           break;
    }
	printf ("PM powerdensity %lf\n", powerdensity);

	resPM += powerdensity;
	resPMF += powerdensityF;

	res2PM +=powerdensity*powerdensity;
	res2PMF+=powerdensityF*powerdensityF;
}

resPM/=iterations;
resPMF/=iterations;

res2PM/=iterations;
res2PMF/=iterations;

double varPM =  res2PM - resPM*resPM;
double varPMF = res2PMF - resPMF*resPMF;

printf ("%i\n", k);
printf ("Averaged PM %lf\n", resPM);
printf ("Averaged PMF %lf\n", resPMF);

printf ("Standard deviation PM %lf\n", sqrt(varPM));
printf ("Standard deviation PMF %lf\n", sqrt(varPMF));

printf ("PM pseudo SNR %lf\n", resPM / sqrt(varPM));
printf ("PMF pseudo SNR %lf\n", resPMF / sqrt(varPMF));

system("pause");
return 0;
}

