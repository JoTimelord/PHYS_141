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
#define kNumDims 3


// global variables
int M=pow(10,11); /* in solar mass */
int G=4.302*pow(10,-6); /* in kiloparsec/solarmass*(km/s)^2 */ 
double R=1.5; /* in kiloparsec */


void position();

int main(argc, argv)
int argc;
char *argv[];
{
    
        
    return 0;
}

void position(r)
double r[];
{
    
}

