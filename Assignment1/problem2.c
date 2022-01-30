/*
 * Problem2.c: compute orbit (x,y) of a planet using leapfrog integration and Kepler's law 
 */

/* Referencd: http://www.gnm.cl/fgonzalez/pmwiki.php/Programas/KeplerOrbits#Animation
 * leapint.c
 */

#include <math.h>
#include <stdio.h>
#include <string.h>

#define MAXPNT 500             /* maximum number of points */

double M;                       /* mass of the sun */
double m;                       /* mass of the planet */
float G; 
void leapstep();                /* routine to take one step */
char planetname[] = "Jupiter";

void accel();                   /* accel. for harmonic osc. */

void printstate();              /* print out system state   */

void startx();                /* set relevant parameters for a planet */

int main(argc, argv)
int argc;
char *argv[];
{
    int n, mstep, nout, nstep;
    double x[MAXPNT], vx[MAXPNT], y[MAXPNT], vy[MAXPNT], tnow, dt;

    G = 6.6743*pow(10,-11);
    /* first, set up initial conditions */
    n = 1;                      /* set number of points     */
    tnow = 0.0;                 /* set initial time         */
    startx(planetname, vy, x);                 /* set initial position     */
    y[0] = 0.0;
    vx[0] = 0.0;

    /* next, set integration parameters */

    mstep = 256*20;                /* number of steps to take  */
    nout = 10;                   /* steps between outputs    */
    dt = 0.01;            /* timestep for integration */

    
    /* now, loop performing integration */

    for (nstep = 0; nstep < mstep; nstep++) {   /* loop mstep times in all  */
        if (nstep % nout == 0)          /* if time to output state  */
            printstate(x, vx, y, vy, n, tnow);      /* then call output routine */
        leapstep(x, y, vx, vy, n, dt);   /* take integration step    */
        tnow = tnow + dt;           /* and update value of time */
    }
    if (mstep % nout == 0)          /* if last output wanted    */
        printstate(x, vx, y, vy, n, tnow);      /* then output last step    */
    
    return 0;
}

void startx(planetname, vy, x) 
char planetname[];
double vy[];
double x[];
{
    double a;  /* semi-major axis */
    double e;  /* eccentricity    */
    G = 6.67408*pow(10,-11);
    if (strcmp(planetname, "Jupiter")==0) {
        e = 0.0485;
        a = 77.8;      /* unit: 1 m */
        vy[0] = 12.44*1000;   /* unit: 1 m/s */
        m = 1898.0;      /* unit: 10^24 kg */
    }
    else if (strcmp(planetname, "Venus")==0) {
        e = 0.0068;
        a = 10.8;
        vy[0] = 34.79*1000;   /* unit: 1 m/s */
        m = 4.87;
    }
    else if (strcmp(planetname, "Mars")==0) {
        e = 0.0934;
        a = 22.8;
        vy[0] = 21.97*1000;   /* unit: 1 m/s */
        m = 0.642;
    }
    else if (strcmp(planetname, "Earth")==0) {
        e = 0.0167;
        a = 15.0;
        vy[0] = 29.29*1000;   /* unit: 1 m/s */
        m = 5.97;
    }
    else if (strcmp(planetname, "Saturn")==0) {
        e = 0.0556;
        a = 143;
        vy[0] = 9.09*1000;   /* unit: 1 m/s */
        m = 568;
    }
    else if (strcmp(planetname, "Uranus")==0) {
        e = 0.0472;
        a = 287;
        vy[0] = 6.49*1000;   /* unit: 1 m/s */
        m = 86.8;
    }
    else if (strcmp(planetname, "Neptune")==0) {
        e = 0.0086;
        a = 450;
        vy[0] = 5.37*1000;   /* unit: 1 m/s */
        m = 102;
    }
    else if (strcmp(planetname, "Pluto")==0) {
        e = 0.25;
        a = 590;
        vy[0] = 3.71*1000;
        m = 0.0130;
    }
    else if (strcmp(planetname, "Mercury")==0) {
        e = 0.206;
        a = 5.79;
        vy[0] = 47.36*1000;
        m = 0.330;
    }
    else { 
        printf("I should never get here. \n");
    }
    x[0] = (a*(1+e))*pow(10.0, 10.0)/1000; /* convert scale into km */
    vy[0] = vy[0]*86400/1000; /* convert velocity into meter per day */
    M = 1.9891*pow(10, 30);
}

/*
 * LEAPSTEP: take one step using the leap-from integrator, formulated
 * as a mapping from t to t + dt.  WARNING: this integrator is not
 * accurate unless the timestep dt is fixed from one call to another.
 */

void leapstep(x, y, vx, vy, n, dt)
double x[];                 /* x-positions of all points  */
double y[];                 /* y-positions of all points  */
double vx[];                 /* x-velocities of all points */
double vy[];                 /* y-velocities of all points */
int n;                      /* number of points         */
double dt;                  /* timestep for integration */
{
    int i;
    double ax[MAXPNT];
    double ay[MAXPNT];

    accel(ax, x, ay, y, n);             /* call acceleration code   */
    for (i = 0; i < n; i++)         /* loop over all points...  */
    {
        vx[i] = vx[i] + 0.5 * dt * ax[i];      /* advance vel by half-step */
        vy[i] = vy[i] + 0.5 * dt * ay[i];      /* advance vel by half-step */
    }
    for (i = 0; i < n; i++)         /* loop over points again...*/
    {
        x[i] = x[i] + dt * vx[i];        /* advance pos by full-step */
        y[i] = y[i] + dt * vy[i];        /* advance pos by full-step */
    }
    accel(ax, x, ay, y, n);             /* call acceleration code   */
    for (i = 0; i < n; i++)         /* loop over all points...  */
    {
        vx[i] = vx[i] + 0.5 * dt * ax[i];      /* and complete vel. step   */
        vy[i] = vy[i] + 0.5 * dt * ay[i];      /* and complete vel. step   */
    }
}

/*
 * ACCEL: compute accelerations for a planet based on Kepler orbit.
 */

void accel(ax, x, ay, y, n)
double ax[];                 /* accelerations of points  */
double ay[];                 /* accelerations of points  */
double x[];                 /* positions of points      */
double y[];                 /* positions of points      */
int n;                      /* number of points         */
{
    int i;
    G = 6.6743*pow(10,-11);
    for (i = 0; i < n; i++)
    {
        /* loop over all points...  */
        ax[i] = - G*M*x[i]/pow(x[i]*x[i]+y[i]*y[i], 3.0/2.0)*86400*86400/1000;                  /* acceleration km per day  */
        ay[i] = - G*M*y[i]/pow(x[i]*x[i]+y[i]*y[i], 3.0/2.0)*86400*86400/1000;                  /* acceleration km per day  */
    }
}

/*
 * PRINTSTATE: output system state variables.
 */

void printstate(x, vx, y, vy, n, tnow)
double x[];                 /* positions of all points  */
double y[];                 /* positions of all points  */
double vx[];                 /* velocities of all points */
double vy[];                 /* velocities of all points */
int n;                      /* number of points         */
double tnow;                    /* current value of time    */
{
    int i;
    for (i = 0; i < n; i++)         /* loop over all points...  */
    printf("%8.4f%4d%23.6f%23.6f%23.6f%23.6f\n", tnow, i, x[i], vx[i], y[i], vy[i]);
}
