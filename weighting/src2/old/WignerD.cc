
#include "WignerD.hh"

#include "TComplex.h"
#include "TMath.h"
#include <iostream>

DPWignerFunctionJ0 WignerD::wignerJ0;
DPWignerFunctionJ1 WignerD::wignerJ1;
DPWignerFunctionJ1over2 WignerD::wignerJh;

TComplex WignerD::WignerDfunction(double J, double m, double n, double cosTheta, double phi)
{

//  std::cout<<"Entering wigner D with J="<<J<<std::endl;

  TComplex result=TComplex::Exp(-(m-n)*TComplex::I()*phi);
  if ( J==0 )
  {
    result*=wignerJ0.function(cosTheta,m,n);
//    std::cout<<"Calling wigner function for J=0\n";
  } 
  else if ( TMath::Abs(J-1)<1e-4 )
  {
    result*=wignerJ1.function(cosTheta,m,n);
//    std::cout<<"Calling wigner function for J=1\n";
  } 
  else 
  {
    result*=wignerJh.function(cosTheta,m,n);
//    std::cout<<"Calling wigner function for J=1/2\n";
  } 

  return result;
}

