/*
 * LEAPINT.C: program to integrate hamiltonian system using leapfrog.
 */

#include <math.h>
#include <stdio.h>
#include <string.h>

#define MAXPNT 200000            /* maximum number of points */

double m;

void leapstep();                /* routine to take one step */

void accel();                   /* accel. for harmonic osc. */

double rand_0_1();              /* generate random numbers */

void initial();                /* generate initial conditions */

void printstate();              /* print out system state   */

void main(argc, argv)
int argc;
char *argv[];
{
    int n, mstep, nout, nstep;
    double x[MAXPNT], y[MAXPNT], v[MAXPNT], tnow, dt;

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

/* set up initial conditions */
void initial(r,v,E_tot)
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



