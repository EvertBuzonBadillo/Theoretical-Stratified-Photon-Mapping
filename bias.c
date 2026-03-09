#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include<time.h>
#include <string.h>

/*Empirical study of the over- and under-estimation biases in stratified Photon Mapping with kernels*/

double alpha = 1.728;
double alphav = 1.97758;
double beta = 1.953;
double* sceneDistances (int n, int k, int dim, double distances[]);
double epanechnikov (double distances[],int k, int dim);
double cone (double distances[], int k, int dim);
double gaussian (double distances[], int k, int dim);
double silverman (double distances[], int k, int dim);

/*Constant terms of the irradiance estimate of the kernels*/

int main()
{
    srand(time(NULL));
    const double EpaConst2D = 2/(double)M_PI;
    const double EpaConst3D = 15/(double)(8*M_PI);
    const double ConeConst2D = 1/(double)(M_PI*(1-(2/3.0)));
    const double ConeConst3D = 3/(double)(M_PI);
    const double GaussConst2D = alpha/(double)M_PI;
    const double GaussConst3D = (3*alphav)/(double)(4*M_PI);
    const double SilvConst2D = 3/(double)M_PI;
    const double SilvConst3D = 105/(double)(32*M_PI);
    int numberPhotons;
    int k;
    int iter;
    int dim;
    char kernel[10];
    printf("Relative error of the irradiance estimate in stratified Photon Mapping \n");
    printf("Enter the number of photons: \n");
    scanf("%d",&numberPhotons);
    printf("Enter a value for k:\n");
    scanf("%d",&k);
    printf("Enter the number of iterations: \n");
    scanf("%d",&iter);
    printf("Select a kernel ->>>> EPA|CONE|GAUSS|SILV: \n ");
    scanf("%s",kernel);
    printf("Enter the dimensions of the kernel ->>>> 2|3: \n ");
    scanf("%d",&dim);
    double power = M_PI/(double)numberPhotons;
    double RPD = 1;
    double constant;
    double irradiance = 0;
    double distances [k];
    int i;
    switch (kernel[0]){
        case 'E':
            for(i=0;i<iter;i++){
              sceneDistances(numberPhotons,k,dim,distances);
              irradiance = irradiance+epanechnikov(distances,k,dim);
            }
            if(dim==2){
                constant = EpaConst2D;
            }
            else {
                constant = EpaConst3D;
                power=power*(4.0/3.0);
            }
            break;
        case 'C':
            for (i=0;i<iter;i++){
             sceneDistances(numberPhotons,k,dim, distances);
             irradiance = irradiance+cone(distances,k,dim);
            }
            if(dim==2){
                constant = ConeConst2D;
            }
            else {
                constant = ConeConst3D;
                power=power*(4.0/3.0);
            }
            break;
        case 'G':
            for(i=0;i<iter;i++){
             sceneDistances(numberPhotons,k,dim, distances);
             irradiance = irradiance+gaussian(distances,k,dim);
            }
            if(dim==2){
                constant = GaussConst2D;
            }
            else {
                constant = GaussConst3D;
                power=power*(4.0/3.0);
            }
            break;
        case 'S':
            for(i=0;i<iter;i++){
              sceneDistances(numberPhotons,k,dim, distances);
              irradiance = irradiance+silverman(distances,k,dim);
            }
            if(dim==2){
                constant = SilvConst2D;
            }
            else {
                constant = SilvConst3D;
                power=power*(4.0/3.0);
            }
            break;
    }
    irradiance = irradiance/(double)iter; //Average irradiance.
    irradiance = irradiance*power*constant;
    double relativeError = (irradiance-RPD)/(double)RPD;
    printf("The relative error is:  \n");
    printf("%.12f \n",relativeError);
    return 0;
}
/*Calculates distances of uniformly distributed random points in the scene*/
double* sceneDistances (int n, int k, int dim, double distances[])
{
    int i;
    for(i = 1;i<=k;i++){
        double uniform = rand()/(double)RAND_MAX;
        double r = (i - 1 +uniform)/(double)n;
        double dist = pow(r,1/(double)dim);
        distances[i-1] = (double)dist;
    }
    return distances;
}
/*Fast Epanechnikov Kernel*/
double epanechnikov (double distances[],int k, int dim)
{
  double firstTerm = (k-1)*distances[k-1]*distances[k-1];
  double secondTerm = 1/(double)(pow(distances[k-1],(double)(dim+2.0)));
  int i;
  double thirdTerm = 0;
  for (i = 0;i<(k-1);i++){
    thirdTerm = thirdTerm + (distances[i]*distances[i]);
  }
  double kernel = (firstTerm-thirdTerm)*secondTerm;
  return kernel;
}
/*Fast Cone Kernel for slope s=1*/
double cone (double distances[], int k, int dim)
{
  double firstTerm = (k-1)*distances[k-1];
  double secondTerm = 1/(double)(pow(distances[k-1],(double)(1.0+dim)));
  int i;
  double thirdTerm = 0;
  for (i = 0;i<(k-1);i++){
    thirdTerm = thirdTerm + distances[i];
  }
  double kernel = (firstTerm-thirdTerm)*secondTerm;
  return kernel;
}
/*Fast Gaussian Kernel*/
double gaussian (double distances[], int k, int dim)
{
    double firstTerm = (-beta/2.0)*(1/(double)(distances[k-1]*distances[k-1]));
    double secondTerm = 0;
    int i;
    for (i = 0; i< k;i++){
        secondTerm = secondTerm + exp((double)firstTerm*(distances[i]*distances[i]));
    }
    double thirdTerm = k * exp((double)(-beta));
    double fourthTerm = (1-exp((double)(-beta)))*pow(distances[k-1],(double)dim);
    double kernel = (secondTerm-thirdTerm)/(double)fourthTerm;
    return kernel;
}
/*Fast Silverman Kernel*/
double silverman (double distances[], int k, int dim)
{
    double firstTerm = -2*distances[k-1]*distances[k-1];
    int i;
    double secondTerm = 0;
    for (i=0;i<k;i++){
        secondTerm = secondTerm+(distances[i]*distances[i]);
    }
    int e;
    double thirdTerm = 0;
    for(e=0;e<(k-1);e++){
        thirdTerm = thirdTerm+pow(distances[e],4.0);
    }
    double fourthTerm = (k+1)*pow(distances[k-1],4.0);
    double fifthTerm = pow(distances[k-1],(double)(4.0+dim));
    double kernel = (firstTerm*secondTerm+thirdTerm+fourthTerm)/(double)fifthTerm;
    return kernel;
}
