#ifndef LBJPSIL2ANG_HH
#define LBJPSIL2ANG_HH

#include "TComplex.h"
#include "TLorentzVector.h"

/*
 * Code has many constants we might in principle want to update hardcoded,
 * but for now I worry about having all formulas correct rather than
 * changing numerical values somewhere.
 */

class LbJpsiL2
{

public:
  LbJpsiL2();

  ~LbJpsiL2() {};

  double physicsRate(double cosTheta, double cosThetaL, double cosThetaB,
                     double phiL, double phiB);
  double physicsRate(double cosTheta, double cosThetaL, double phiL);

  void setPolarization(double polarization)
  {
    Pb = polarization;
  };
  void setLambdaAsymmetry(double as)
  {
    alphaLambda = as;
  };

  void setAlphaLambda(double al)
  {
    alphaLambda = al;
  }

  void setAlphabR0R1(double alphaB, double R0, double R1, double apPhase = 0,
                     double amPhase = 0, double bpPhase = 0,
                     double bmPhase = 0);

  void setVariation(int variation);

  static void LbPsiRAngles(TLorentzVector initialProton, TLorentzVector pB,
                              TLorentzVector pJpsi, TLorentzVector pLambda,
                              TLorentzVector pmp, TLorentzVector pmm,
                              TLorentzVector pp, TLorentzVector ppi,
                              int pCharge, double &cosTheta, double &cosThetaL,
                              double &cosThetaB, double &phiL, double &phiB,
                              double &dphi);
private:
  double Pb, alphaLambda;

  double hLepton(int J, double ll);
  double spinDensityMatrix(double lLb, double lLbp);
  TComplex hLb(double lL, int lj);
  double hLambda(double lp);

  double spinDensity[2][2];
  void updateSpinDensity(double cosTheta);

  TComplex aplus, aminus, bplus, bminus;
};

#endif
