#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define dot(a,b) ((a[0])*(b[0]) + (a[1])*(b[1]) + (a[2])*(b[2]))
#define lengthSquared(a) dot(a,a)
#define length(a) sqrt(lengthSquared(a))

#define cross(c,a,b) c[0] = (a[1])*(b[2]) - (a[2])*(b[1]); \
                     c[1] = (a[2])*(b[0]) - (a[0])*(b[2]); \
                     c[2] = (a[0])*(b[1]) - (a[1])*(b[0]);

#define kMaxParticles 100000
#define PI 3.141592654
#define kNumDims 3


// global variables
// define a new mass unit = solarmass/(4.302*pow(10,-6))
double M=pow(10.0,11)*4.302*pow(10,-6); /* in new mass */
double G=1; /* in kiloparsec/newmass*(km/s)^2 */ 
double R=1.5; /* in kiloparsec */
double E;

void position();
double rand_0_1();
double g();


int main(argc, argv)
int argc;
char *argv[];
{
    E=-3.0*PI/64.0*G*M*M/R; /* in newmass*(km/s)^2 */
    FILE *fp;
    fp=fopen("1a.dat","w+");
    double r[4];
    double v[5];
    double E_tot[1];
    fprintf(fp,"%-7s%-14.4s%-14.4s%-14.4s%-14.4s%-14.4s%-14.4s%-14.4s%-14.4s%-14.4s\n","index","r component","x component","y component","z component","speed","v component","w component","u component","Energy");
    for (int i=0;i<kMaxParticles;i++) 
    {
        position(r,v,E_tot);
        fprintf(fp,"%-7d%-14.4f%-14.4f%-14.4f%-14.4f%-14.4f%-14.4f%-14.4f%-14.4f%-14.4f\n",i,r[0],r[1],r[2],r[3],v[1],v[2],v[3],v[4],E_tot[0]); /* this prints out the index of the star, r, x, y, z, V, v, w, u, v component of the star */
    }
    fclose(fp); 
        
    return 0;
}

void position(r,v,E_tot)
double r[];
double v[];
double E_tot[];
{
    double q;
    double X1, X2, X3, X4, X5, X6, X7;
    double U;
    X1=rand_0_1();
    r[0]=pow(pow(X1,-2.0d/3.0d)-1,-0.5); /* radius */
    X2=rand_0_1();
    r[1]=(1-2*X2)*r[0]; /*  z component */
    X3=rand_0_1();
    r[2]=pow(r[0]*r[0]-r[1]*r[1],0.5)*cos(2*PI*X3); /* x component */
    r[3]=pow(r[0]*r[0]-r[1]*r[1],0.5)*sin(2*PI*X3); /* y component */
    v[0]=pow(2,0.5)*pow(1+r[0]*r[0]/(R**2),-0.25)*sqrt(G*M/R); /* escape velocity */
    X4=rand_0_1();
    X5=rand_0_1();
    while (0.1*X5>=g(X4)) {
        X4=rand_0_1();
        X5=rand_0_1();
    }
    q=X4;
    X6=rand_0_1();
    X7=rand_0_1();
    v[1]=q*v[0]; /* velocity component */
    v[2]=(1-2*X6)*v[1]; /* w component */
    v[3]=pow(v[1]*v[1]-v[2]*v[2],0.5)*cos(2*PI*X7); /* u component */
    v[4]=pow(v[1]*v[1]-v[2]*v[2],0.5)*sin(2*PI*X7); /* v component */
    for (int k=0;k<4;k++)
    {
        r[k]=r[k]*3*PI/64.0*pow(M,2)/(-E); /* in kiloparsec */
    }
    U=-4.302*pow(10,5)/1.5*pow(1+pow(r[0]/1.5,2),-0.5); /* in (km/s)^2 */
    E_tot[0]=U+v[1]*v[1]/2;
}

double g(q)
double q;
{
   return q*q*pow(1-q*q, 3.5); 
}

double rand_0_1(void)
{
    return (double) rand()/((double) RAND_MAX);
}
