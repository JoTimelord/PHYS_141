/*
 * LEAPINT.C: program to integrate hamiltonian system using leapfrog.
 */

#include <math.h>
#include <stdio.h>

#define MAXPNT 100             /* maximum number of points */

void leapstep();                /* routine to take one step */

void accel();                   /* accel. for harmonic osc. */

void printstate();              /* print out system state   */

void main(argc, argv)
int argc;
char *argv[];
{
    int n, mstep, nout, nstep;
    double x[MAXPNT], v[MAXPNT], tnow, dt;

    /* first, set up initial conditions */

    n = 1;                      /* set number of points     */
    x[0] = 1.0;                 /* set initial position     */
    v[0] = 3.0;                 /* set initial velocity     */
    tnow = 0.0;                 /* set initial time         */

    /* next, set integration parameters */

    mstep = 256;                /* number of steps to take  */
    nout = 4;                   /* steps between outputs    */
    dt = 1.0 / 32.0;            /* timestep for integration */

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

    for (i = 0; i < n; i++)         /* loop over all points...  */
    a[i] = - x[i];                  /* use linear force law     */
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
        printf("I should never get here. \n");
    }
    x[0] = a*(1+e)*pow(10,10)/(1.496*pow(10,11)); /* starting position in AU */
    vy[0] = vy[0]/(1.496*pow(10,11))*3.156*pow(10,7); 
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
    int i;

    for (i = 0; i < n; i++)         /* loop over all points...  */
    printf("%8.4f%4d%12.6f%12.6f\n", tnow, i, x[i], v[i]);
}
