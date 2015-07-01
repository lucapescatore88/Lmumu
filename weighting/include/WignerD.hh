#ifndef WIGNERD_HH
#define WIGNERD_HH

#include "TComplex.h"
#include "DPWignerFunctionJ0.hh"
#include "DPWignerFunctionJ1.hh"
#include "DPWignerFunctionJ1over2.hh"

class WignerD
{

  public:
 
    static TComplex WignerDfunction(double J, double m, double n, double cosTheta, double phi);

  private:

    static DPWignerFunctionJ0 wignerJ0;
    static DPWignerFunctionJ1 wignerJ1;
    static DPWignerFunctionJ1over2 wignerJh;

};

#endif
