/*
 * LEAPINT.C: program to integrate hamiltonian system using leapfrog.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAXPNT 400000            /* maximum number of points */
#define PI 3.141592654

double M=pow(10.0,11)*4.302*pow(10,-6); /* in new mass */
double G=4.302*pow(10,-6); /* in kiloparsec/newmass*(km/s)^2 */ 
double R=1.5; /* in kiloparsec */
double m=pow(10.0,11)*4.302*pow(10,-6)/200000; /* mass per star */

void leapstep();                /* routine to take one step */

void accel();                   /* accel. for harmonic osc. */

double rand_0_1();              /* generate random numbers */

void initial();                 /* generate initial conditions */

double g();                     /* generate q distribution */

void printstate();              /* print out system state   */

void main(argc, argv)
int argc;
char *argv[];
{
    int n, mstep, nout, nstep;
    n = 20000;
    double x[n], y[n], z[n], r[n], V[n], v[n], w[n], u[n], E[n], tnow, dt;

    /* store init.dat */
    FILE *fp;
    fp=fopen("init.dat","w+");

    /* first, set up initial conditions */
    initial(r,x,y,z,V,v,w,u,E,n); /* set initial vel and posi of all points */
    fprintf(fp,"%-7s%-14.4s%-14.4s%-14.4s%-14.4s%-14.4s%-14.4s%-14.4s%-14.4s%-14.4s\n","index","r component","x component","y component","z component","speed","v component","w component","u component","Energy");
    for (int i=0;i<n;i++) 
    {
        fprintf(fp,"%-7d%-14.4f%-14.4f%-14.4f%-14.4f%-14.4f%-14.4f%-14.4f%-14.4f%-14.4f\n",i,r[i],x[i],y[i],z[i],V[i],v[i],w[i],u[i],E[i]); /* this prints out the index of the star, r, x, y, z, V, v, w, u, v component of the star */
    }
    fclose(fp); 
// 
//    tnow = 0.0;                 /* set initial time         */
//
//    /* next, set integration parameters */
//    mstep = 300;                /* number of steps to take  */
//    nout = 1;                   /* steps between outputs    */
//    dt = 0.07;            /* timestep for integration, in year */
//
//    /* now, loop performing integration */
//
//    for (nstep = 0; nstep < mstep; nstep++) {   /* loop mstep times in all  */
//    if (nstep % nout == 0)          /* if time to output state  */
//        printstate(x, v, n, tnow);      /* then call output routine */
//    leapstep(r, x, y, z, V, w, v, u, n, dt); /* take integration step    */
//    tnow = tnow + dt;           /* and update value of time */
//    }
//    if (mstep % nout == 0)          /* if last output wanted    */
//    printstate(x, v, n, tnow);      /* then output last step    */
}

/* set up initial conditions 
 * The positions are in kiloparsecs
 * The energies are in (km/s)^2
 * The velocities are in (km/s)
 * */
void initial(r,x,y,z,V,v,w,u,E,n)
double r[];
double x[];
double y[];
double z[];
double V[];
double v[];
double u[];
double w[];
double E[];
int n;
{
    for (int i=0;i<n;i++){
        double q;
        double X1, X2, X3, X4, X5, X6, X7;
        double U;
        double Ve;
        X1=rand_0_1();
        double K=pow(X1*pow(R,3)/M,2.0d/3.0d);
        r[i]=pow(K*R*R/(R*R-K),0.5); /* radius */
        X2=rand_0_1();
        z[i]=(1-2*X2)*r[i]; /*  z component */
        X3=rand_0_1();
        x[i]=pow(r[i]*r[i]-r[i]*r[i],0.5)*cos(2*PI*X3); /* x component */
        y[i]=pow(r[i]*r[i]-r[i]*r[i],0.5)*sin(2*PI*X3); /* y component */
        Ve=pow(2,0.5)*pow(1+r[i]*r[i]/(R*R),-0.25)*sqrt(G*M/R); /* escape velocity */
        X4=rand_0_1();
        X5=rand_0_1();
        while (0.1*X5>=g(X4)) {
            X4=rand_0_1();
            X5=rand_0_1();
        }
        q=X4;
        X6=rand_0_1();
        X7=rand_0_1();
        V[i]=q*Ve; /* velocity component */
        w[i]=(1-2*X6)*V[i]; /* w component */
        u[i]=pow(V[i]*V[i]-w[i]*w[i],0.5)*cos(2*PI*X7); /* u component */
        v[i]=pow(V[i]*V[i]-w[i]*w[i],0.5)*sin(2*PI*X7); /* v component */
        U=-G*M/R*pow(1+(r[i]/R)*(r[i]/R),-0.5);
        E[i]=U+V[i]*V[i];
    }
}

/*
 * Probability distribution functions
 */
double g(q)
double q;
{
   return q*q*pow(1-q*q, 3.5); 
}

double rand_0_1(void)
{
    return (double) rand()/((double) RAND_MAX);
}

/*
 * LEAPSTEP: take one step using the leap-from integrator, formulated
 * as a mapping from t to t + dt.  WARNING: this integrator is not
 * accurate unless the timestep dt is fixed from one call to another.
 */

void leapstep(r, x, y, z, V, w, v, u, n, dt)
double r[];                 /* r-positions of all points  */
double x[];                 /* x-positions of all points  */
double y[];                 /* y-positions of all points  */
double z[];                 /* z-positions of all points  */
double V[];                 /* velocities of all points */
double w[];                 /* z velocities of all points */
double v[];                 /* y velocities of all points */
double u[];                 /* x velocities of all points */
int n;                      /* number of points         */
double dt;                  /* timestep for integration */
{
    int i;
    double a[n],ax[n],ay[n],az[n];
    accel(a, ax, ay, az, r, x, y, z, n); /* call acceleration code   */
    for (i = 0; i < n; i++)         /* loop over all points...  */
    {
        V[i] = V[i] + 0.5 * dt * a[i];      /* advance vel by half-step */
        w[i] = w[i] + 0.5 * dt * az[i];      
        u[i] = u[i] + 0.5 * dt * ax[i];      
        v[i] = v[i] + 0.5 * dt * ay[i];      
    }
    for (i = 0; i < n; i++)         /* loop over points again...*/
    {
        r[i] = r[i] + dt * V[i];        /* advance pos by full-step */
        x[i] = x[i] + dt * u[i];        
        y[i] = y[i] + dt * v[i];        
        z[i] = z[i] + dt * w[i];        
    }
    accel(a, ax, ay, az, r, x, y, z, n);             /* call acceleration code   */
    for (i = 0; i < n; i++)             /* loop over all points...  */
    {
        V[i] = V[i] + 0.5 * dt * a[i];      /* advance vel by half-step */
        w[i] = w[i] + 0.5 * dt * az[i];      
        u[i] = u[i] + 0.5 * dt * ax[i];      
        v[i] = v[i] + 0.5 * dt * ay[i];      
    }
}

/*
 * ACCEL: compute accelerations for harmonic oscillator(s).
 */

void accel(a, ax, ay, az, r, x, y, z, n)
double a[];                 /* accelerations of points  */
double ax[];                 /* x accelerations of points  */
double ay[];                 /* y accelerations of points  */
double az[];                 /* z accelerations of points  */
double r[];                 /* r acceleration of points */
double x[];                 /* x acceleration of points */
double y[];                 /* y acceleration of points */
double z[];                 /* z acceleration of points */
int n;                      /* index of points         */
{
    for (int i=0;i<n;i++){ /* calculate acceleration for every point i */
        for (int j=0;j<n;j++){ /* calculate acceleration exerted on i by every other point */
            if (i!=j){
                double rij=sqrt(pow(x[i]-x[j],2)+pow(y[i]-y[j],2)+pow(z[i]-z[j],2));
                a[i]=a[i]-G*m/pow(rij,2);
                az[i]=az[i]+a[i]*(z[j]-z[i])/rij;
                ax[i]=ax[i]+a[i]*(x[j]-x[i])/rij;
                ay[i]=ay[i]+a[i]*(y[j]-y[i])/rij;
            }
        }
    }
}

/*
 * PRINTSTATE: output system state variables.
 */

void printstate(x, v, n, tnow)
double x[];                 /* positions of all points  */
double v[];                 /* velocities of all points */
int n;                      /* number of points         */
double tnow;                    /* current value of time    */
{

    printf("%8.4f%12.6f%12.6f%12.6f%12.6f\n", tnow, x[0], v[0], x[1], v[1]);
}



