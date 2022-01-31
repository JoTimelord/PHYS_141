//
// This program solves the differential equation
//             m R''(t) = - GMm/r² r^
// where R=(x,y) and r^ = (x,y)/sqrt(x²+y²)
// using Runge-Kutta methods.

#include<iostream>
#include<fstream>
#include<cmath>
using namespace std;

const double G=0.00011859645; // in AU³/(earthmass*year²)
const double m=1;             // earthmass
const double M=332946.05;     // sunmass=332946.05 earthmass
const double tau=0.01, tmax=5; // year

typedef struct vect{
 double x, y, z;
} Vector;


// We separate the 2nd order vector equation in 4 first-order differential equations:
// x1(t) = x(t)
// y1(t) = y(t)
// x2(t) = x'(t)
// y2(t) = y'(t)

double fx1(double  x1, double x2, double y1, double y2, double t){ return x2; } // x1'(t) = x2(t) = fx1(x1,x2,y1,y2,t)
double fy1(double  x1, double x2, double y1, double y2, double t){ return y2; } // y1'(t) = y2(t) = fx2(x1,x2,y1,y2,t)
double fx2(double  x1, double x2, double y1, double y2, double t)
{
 return -G*M*x1/pow(x1*x1+y1*y1,1.5); // x2'(t) = -G*M*x1/R³ = fy1(x1,x2,y1,y2,t)
}
double fy2(double  x1, double x2, double y1, double y2, double t)
{
 return -G*M*y1/pow(x1*x1+y1*y1,1.5); // y2'(t) = -G*M*y1/R³ = fy1(x1,x2,y1,y2,t)
}

Vector RungeLenz(double x, double y, double vx, double vy)
{
 double r=sqrt(x*x+y*y);
 double px=m*vx, py=m*vy;
 double Lz=x*py-y*px;
 Vector A;
 A.x =  Lz*py-m*(G*M*m)*x/r;
 A.y = -Lz*px-m*(G*M*m)*y/r;
 return A;
}

///////////////////////////////
// METHODS                   //
///////////////////////////////

// Runge-Kutta first order
void RK1(double &x1, double &x2, double &y1, double &y2, double &t)
{
 double k1x1, k1x2, k1y1, k1y2;

 k1x1 = tau*fx1(x1,x2,y1,y2,t);
 k1x2 = tau*fx2(x1,x2,y1,y2,t);
 k1y1 = tau*fy1(x1,x2,y1,y2,t);
 k1y2 = tau*fy2(x1,x2,y1,y2,t);

 x1 = x1 + k1x1;
 x2 = x2 + k1x2;
 y1 = y1 + k1y1;
 y2 = y2 + k1y2;
}

// Runge-Kutta de second order
void RK2(double &x1, double &x2, double &y1, double &y2, double &t)
{
 double k1x1, k1x2, k1y1, k1y2;
 double k2x1, k2x2, k2y1, k2y2;

 k1x1 = tau*fx1(x1,x2,y1,y2,t);
 k1x2 = tau*fx2(x1,x2,y1,y2,t);
 k1y1 = tau*fy1(x1,x2,y1,y2,t);
 k1y2 = tau*fy2(x1,x2,y1,y2,t);

 k2x1 = tau*fx1(x1+k1x1,x2+k1x2,y1+k1y1,y2+k1y2,t);
 k2x2 = tau*fx2(x1+k1x1,x2+k1x2,y1+k1y1,y2+k1y2,t);
 k2y1 = tau*fy1(x1+k1x1,x2+k1x2,y1+k1y1,y2+k1y2,t);
 k2y2 = tau*fy2(x1+k1x1,x2+k1x2,y1+k1y1,y2+k1y2,t);

 x1 = x1 + 0.5*(k1x1+k2x1);
 x2 = x2 + 0.5*(k1x2+k2x2);
 y1 = y1 + 0.5*(k1y1+k2y1);
 y2 = y2 + 0.5*(k1y2+k2y2);
}

// Runge-Kutta third order
void RK3(double &x1, double &x2, double &y1, double &y2, double &t)
{
 double k1x1, k1x2, k1y1, k1y2;
 double k2x1, k2x2, k2y1, k2y2;
 double k3x1, k3x2, k3y1, k3y2;

 k1x1 = tau*fx1(x1,x2,y1,y2,t);
 k1x2 = tau*fx2(x1,x2,y1,y2,t);
 k1y1 = tau*fy1(x1,x2,y1,y2,t);
 k1y2 = tau*fy2(x1,x2,y1,y2,t);

 k2x1 = tau*fx1(x1+0.5*k1x1, x2+0.5*k1x2, y1+0.5*k1y1, y2+0.5*k1y2, t+0.5*tau);
 k2x2 = tau*fx2(x1+0.5*k1x1, x2+0.5*k1x2, y1+0.5*k1y1, y2+0.5*k1y2, t+0.5*tau);
 k2y1 = tau*fy1(x1+0.5*k1x1, x2+0.5*k1x2, y1+0.5*k1y1, y2+0.5*k1y2, t+0.5*tau);
 k2y2 = tau*fy2(x1+0.5*k1x1, x2+0.5*k1x2, y1+0.5*k1y1, y2+0.5*k1y2, t+0.5*tau);

 k3x1 = tau*fx1(x1-k1x1+2*k2x1, x2-k1x2+2*k2x2, y1-k1y1+2*k2y1, y2-k1y2+2*k2y2, t+tau);
 k3x2 = tau*fx2(x1-k1x1+2*k2x1, x2-k1x2+2*k2x2, y1-k1y1+2*k2y1, y2-k1y2+2*k2y2, t+tau);
 k3y1 = tau*fy1(x1-k1x1+2*k2x1, x2-k1x2+2*k2x2, y1-k1y1+2*k2y1, y2-k1y2+2*k2y2, t+tau);
 k3y2 = tau*fy2(x1-k1x1+2*k2x1, x2-k1x2+2*k2x2, y1-k1y1+2*k2y1, y2-k1y2+2*k2y2, t+tau);

 x1 = x1 +(k1x1+4*k2x1+k3x1)/6.0;
 x2 = x2 +(k1x2+4*k2x2+k3x2)/6.0;
 y1 = y1 +(k1y1+4*k2y1+k3y1)/6.0;
 y2 = y2 +(k1y2+4*k2y2+k3y2)/6.0;
}

double KineticEnergy(double x1, double x2, double y1, double y2, double t)
{
 return 0.5*m*(x2*x2+y2*y2);       // in earthmass*(AU/year)²
}

double PotentialEnergy(double x1, double x2, double y1, double y2, double t)
{
 return -G*M*m/sqrt(x1*x1+y1*y1);  // in earthmass*(AU/year)²
}


int main()
{
 double x1, x2, y1, y2, t; // x1= x, x2=vx, y1=y, y2=vy
 ofstream archivo("sol-tierraRK1.dat");
 for(x1=5, x2=0, y1=0, y2=5, t=0; t<tmax; t+=tau)
 {
  archivo << t << " " << x1 << " " << x2 << " " << y1 << " " << y2 << endl;
  RK1(x1,x2,y1,y2,t);
 }
 archivo.close();

 ofstream archivo2("sol-tierraRK2.dat");
 for(x1=5, x2=0, y1=0, y2=5, t=0; t<tmax; t+=tau)
 {
  archivo2 << t << " " << x1 << " " << x2 << " " << y1 << " " << y2 << endl;
  RK2(x1,x2,y1,y2,t);
 }
 archivo2.close();

 ofstream archivo3("sol-tierraRK3.dat");
 ofstream prop("props.dat");
//cout << "set terminal " << endl; //wxt size 600,600" << endl; // FIX WINDOW SIZE FOR GNUPLOT
 double A1, A2;
 for(x1=1, x2=0, y1=0, y2=7.28, t=0; t<tmax; t+=tau)
 {
  archivo3 << t << " " << x1 << " " << x2 << " " << y1 << " " << y2 << endl;
  A1 = RungeLenz(x1,y1,x2,y2).x;
  A2 = RungeLenz(x1,y1,x2,y2).y;
  prop << t << " " << KineticEnergy(x1,x2,y1,y2,t) << " " << PotentialEnergy(x1,x2,y1,y2,t)  << " " << A1 << " " << A2 << endl;

 // cout<<"plot [-2:2][-2:2] sqrt(1-x**2)  lc 1 t 'Circulo', -sqrt(1-x**2) lc 1 notitle, '-' ls 7 lc 6 ps 8 t 'Sol', '-' ls 7 lc 3 ps 2 t 'Tierra'" << endl;
 // cout << 0 << " "<< 0 << endl <<'e'<< endl ;
 // cout << x1 << " "<< y1 << endl <<'e'<< endl << flush;

  RK3(x1,x2,y1,y2,t);
 }
 archivo3.close();
 prop.close();

 return 0;
}
