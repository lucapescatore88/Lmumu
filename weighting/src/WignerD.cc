
#include "WignerD.hh"

#include "TComplex.h"
#include "TMath.h"
#include <iostream>

DPWignerFunctionJ0 WignerD::wignerJ0;
DPWignerFunctionJ1 WignerD::wignerJ1;
DPWignerFunctionJ1over2 WignerD::wignerJh;

TComplex WignerD::WignerDfunction(double J, double m, double n, double cosTheta,
                                  double phi)
{

  TComplex result = TComplex::Exp(-(m - n) * TComplex::I() * phi);
  if (J == 0)
  {
    result *= wignerJ0.function(cosTheta, m, n);
  }
  else if (TMath::Abs(J - 1) < 1e-4)
  {
    result *= wignerJ1.function(cosTheta, m, n);
  }
  else
  {
    result *= wignerJh.function(cosTheta, m, n);
  }

  return result;
}

