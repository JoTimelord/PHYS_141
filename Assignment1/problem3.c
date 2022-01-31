#include <stdio.h>
#include <string.h>
#include <math.h>

#define PI 3.14159265
double a; // in AU
double e;
char planetname[] = "Mars";


void planetinfo();
void analytic();
void printstate();

void main()
{
    double theta[3], r[3], x[3], y[3];
    theta[0]=0;
    double npoint=300;
    double dtheta=2*PI/npoint;
    planetinfo(planetname);
    for (int i=0; i<npoint; i++){
        analytic(r,theta,x,y,dtheta);
        printstate(x,y,theta);
    }
}

void analytic(r, theta, x, y, dtheta)
double r[];
double theta[];
double x[];
double y[];
double dtheta;
{
    r[0]=a*(1-pow(e,2))/(1-e*cos(theta[0])); 
    x[0]=r[0]*cos(theta[0]);
    y[0]=r[0]*sin(theta[0]);
    theta[0]=theta[0]+dtheta;
}

void printstate(x,y)
double x[];                 /* positions of all points  */
double y[];                 /* velocities of all points */
{

    printf("%12.6f%12.6f\n", x[0], y[0]);
}


void planetinfo(planetname)
char planetname[];
{
    if (strcmp(planetname, "Earth")==0)
    {
        a=15.0*pow(10,10)/(1.496*pow(10,11));
        e=0.0167;
    }
    else if (strcmp(planetname, "Mars")==0)
    {
        a=22.8*pow(10,10)/(1.496*pow(10,11));
        e=0.0934;
    }
    else {printf("I should never get here. \n");}

}
















