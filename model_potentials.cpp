#include <iostream>
#include <math.h>

using namespace std;

/* three-hole potential, shifted so that the two deep minima are at x[1] = -2./3. x should be 2-dimensional.
   Suggested domain: -2 < x[0] < 2; -2 < x[1] < 2 */
double three_hole_func (double *x) {

    return 3.*exp(-pow(x[0],2.)-pow(x[1]+(1./3.),2.)) - 3.*exp(-pow(x[0],2.)-pow(x[1]-1.,2.)) \
           - 5.*exp(-pow(x[0]-1.,2.)-pow(x[1]+(2./3.),2.)) - 5.*exp(-pow(x[0]+1.,2.)-pow(x[1]+(2./3.),2.)) \
           + 0.2*pow(x[0],4.) + 0.2*pow(x[1],4.);
}

/* bird function. x should be 2-dimensional.
   two global minima f=-106.764537 at (4.70104,3.15294) and (-1.58214,-3.13024).
   suggested domain: any, depending on desired number of metastable states */
double bird_func (double *x) {

    return (sin(x[0])*exp(pow(1.-cos(x[1]),2.)))+(cos(x[1])*exp(pow(1.-sin(x[0]),2.)))+pow(x[0]-x[1],2.);
}

/* egg crate function. x should be 2-dimensional.
   global minimum f=0. at (0.,0.)
   suggested domain: any, depending on desired number of metastable states */
double egg_crate_func (double *x) {

    return pow(x[0],2.)+pow(x[1],2.)+(25.*(pow(sin(x[0]),2.)+pow(sin(x[1]),2.)));
}

/* Rastrigin function. x should be 2-dimensional.
   global minimum f=0. at (0.,0.)
   suggested domain: any, depending on desired number of metastable states */
double rastrigin_func (double *x) {

    return 20.+(pow(x[0],2.)-(10.*cos(2.*4.*atan(1.)*x[0])))+(pow(x[1],2.)-(10.*cos(2.*4.*atan(1.)*x[1])));
}

int main() {

    double x[] = {4.70104,3.15294};
    cout << bird_func(x) << endl;
    x[0]=-1.58214;x[1]=-3.13024;
    cout << bird_func(x) << endl;
    x[0]=0.;x[1]=0.;
    cout << bird_func(x) << endl;
    x[0]=-0.554887;x[1]=0.463602;
    cout << bird_func(x) << endl;
    x[0]=0.,x[1]=0.;
    cout << egg_crate_func(x) << endl;
    cout << rastrigin_func(x) << endl;

    return 0;
}
