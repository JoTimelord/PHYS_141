/*
 * LEAPINT.C: program to integrate hamiltonian system using leapfrog.
 */

#include <math.h>
#include <stdio.h>
#include <string.h>

#define MAXPNT 100             /* maximum number of points */

char planetname[] = "Mercury";

double m;

void leapstep();                /* routine to take one step */

void accel();                   /* accel. for harmonic osc. */

void printstate();              /* print out system state   */

void startx();

void main(argc, argv)
int argc;
char *argv[];
{
    int n, mstep, nout, nstep;
    double x[MAXPNT], v[MAXPNT], tnow, dt;

    /* first, set up initial conditions */

    n = 2;                      /* set number of points     */
    v[0] = 0.0;                 /* set initial velocity     */
    x[1] = 0;
    startx(planetname, v, x); 
    tnow = 0.0;                 /* set initial time         */

    /* next, set integration parameters */

    mstep = 300;                /* number of steps to take  */
    nout = 1;                   /* steps between outputs    */
    dt = 0.07;            /* timestep for integration, in year */

    /* now, loop performing integration */

    for (nstep = 0; nstep < mstep; nstep++) {   /* loop mstep times in all  */
    if (nstep % nout == 0)          /* if time to output state  */
        printstate(x, v, n, tnow);      /* then call output routine */
    leapstep(x, v, n, dt);          /* take integration step    */
    tnow = tnow + dt;           /* and update value of time */
    }
    if (mstep % nout == 0)          /* if last output wanted    */
    printstate(x, v, n, tnow);      /* then output last step    */
}

/*
 * LEAPSTEP: take one step using the leap-from integrator, formulated
 * as a mapping from t to t + dt.  WARNING: this integrator is not
 * accurate unless the timestep dt is fixed from one call to another.
 */

void leapstep(x, v, n, dt)
double x[];                 /* positions of all points  */
double v[];                 /* velocities of all points */
int n;                      /* number of points         */
double dt;                  /* timestep for integration */
{
    int i;
    double a[MAXPNT];

    accel(a, x, n);             /* call acceleration code   */
    for (i = 0; i < n; i++)         /* loop over all points...  */
    v[i] = v[i] + 0.5 * dt * a[i];      /* advance vel by half-step */
    for (i = 0; i < n; i++)         /* loop over points again...*/
    x[i] = x[i] + dt * v[i];        /* advance pos by full-step */
    accel(a, x, n);             /* call acceleration code   */
    for (i = 0; i < n; i++)         /* loop over all points...  */
    v[i] = v[i] + 0.5 * dt * a[i];      /* and complete vel. step   */
}

/*
 * ACCEL: compute accelerations for harmonic oscillator(s).
 */

void accel(a, x, n)
double a[];                 /* accelerations of points  */
double x[];                 /* positions of points      */
int n;                      /* number of points         */
{
    int i;
    double G = 0.00011859645;   /* in AU^3/earthmass/year^2 */ 
    double M = 332946.05;       /* in earthmass */
    double r = pow(x[0]*x[0]+x[1]*x[1], 0.5);
    for (i = 0; i < n; i++)         /* loop over all points...  */
        a[i] = - G*M*x[i]/pow(r, 1.5);                  /* use linear force law     */
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


/* findout planets info
 *
 */
void startx(planetname, v, x) 
char planetname[];
double v[];
double x[];
{
    double a;  /* semi-major axis */
    double e;  /* eccentricity    */
    if (strcmp(planetname, "Jupiter")==0) {
        e = 0.0485;
        a = 77.8;      /* unit: 1 m */
        v[1] = 12.44*1000;   /* unit: 1 m/s */
        m = 1898.0;      /* unit: 10^24 kg */
    }
    else if (strcmp(planetname, "Venus")==0) {
        e = 0.0068;
        a = 10.8;
        v[1] = 34.79*1000;   /* unit: 1 m/s */
        m = 4.87;
    }
    else if (strcmp(planetname, "Mars")==0) {
        e = 0.0934;
        a = 22.8;
        v[1] = 21.97*1000;   /* unit: 1 m/s */
        m = 0.642;
    }
    else if (strcmp(planetname, "Earth")==0) {
        e = 0.0167;
        a = 15.0;
        v[1] = 29.29*1000;   /* unit: 1 m/s */
        m = 5.97;
    }
    else if (strcmp(planetname, "Saturn")==0) {
        e = 0.0556;
        a = 143;
        v[1] = 9.09*1000;   /* unit: 1 m/s */
        m = 568;
    }
    else if (strcmp(planetname, "Uranus")==0) {
        e = 0.0472;
        a = 287;
        v[1] = 6.49*1000;   /* unit: 1 m/s */
        m = 86.8;
    }
    else if (strcmp(planetname, "Neptune")==0) {
        e = 0.0086;
        a = 450;
        v[1] = 5.37*1000;   /* unit: 1 m/s */
        m = 102;
    }
    else if (strcmp(planetname, "Pluto")==0) {
        e = 0.25;
        a = 590;
        v[1] = 3.71*1000;
        m = 0.0130;
    }
    else if (strcmp(planetname, "Mercury")==0) {
        e = 0.206;
        a = 5.79;
        v[1] = 47.36*1000;
        m = 0.330;
    }
    else {
        printf("I should never get here. \n");
    }
    x[0] = a*(1+e)*pow(10,10)/(1.496*pow(10,11)); /* starting position in AU */
    v[1] = v[1]/(1.496*pow(10,11))*3.156*pow(10,7); 
}


