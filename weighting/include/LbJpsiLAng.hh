#ifndef LBJPSILANG_HH
#define LBJPSILANG_HH

#include "TComplex.h"
#include "TLorentzVector.h"

class LbJpsiLAng
{

public:
  LbJpsiLAng();
  LbJpsiLAng(double polarization, double as, double alphaB, double R0,
             double R1, double apPhase = 0, double amPhase = 0,
             double bpPhase = 0, double bmPhase = 0);

  ~LbJpsiLAng() {};

  double physicsRate(double cosTheta, double cosThetaL, double cosThetaB,
                     double phiL, double phiB);

  void setAlphabR0R1(double alphaB, double R0, double R1, double apPhase = 0,
                     double amPhase = 0, double bpPhase = 0,
                     double bmPhase = 0);

  void setPolarization(double polarization)
  {
    Pb = polarization;
  };
  void setLambdaAsymmetry(double as)
  {
    alphaLambda = as;
  };

  static void LbPsiRAngles(TLorentzVector initialProton, TLorentzVector pB,
                           TLorentzVector pJpsi, TLorentzVector pLambda,
                           TLorentzVector pmp, TLorentzVector pmm,
                           TLorentzVector pp, TLorentzVector ppi, int pCharge,
                           double &cosTheta, double &cosThetaL,
                           double &cosThetaB, double &phiL, double &phiB,
                           double &dphi);

  void setVariation(int variation);

private:
  double term0(TComplex ap, TComplex am, TComplex bp, TComplex bm, double Pb,
               double alphaLambda, double cosTheta1, double cosTheta2,
               double phi1, double phi2, double cosTheta);
  double term1(TComplex ap, TComplex am, TComplex bp, TComplex bm, double Pb,
               double alphaLambda, double cosTheta1, double cosTheta2,
               double phi1, double phi2, double cosTheta);
  double term2(TComplex ap, TComplex am, TComplex bp, TComplex bm, double Pb,
               double alphaLambda, double cosTheta1, double cosTheta2,
               double phi1, double phi2, double cosTheta);
  double term3(TComplex ap, TComplex am, TComplex bp, TComplex bm, double Pb,
               double alphaLambda, double cosTheta1, double cosTheta2,
               double phi1, double phi2, double cosTheta);
  double term4(TComplex ap, TComplex am, TComplex bp, TComplex bm, double Pb,
               double alphaLambda, double cosTheta1, double cosTheta2,
               double phi1, double phi2, double cosTheta);
  double term5(TComplex ap, TComplex am, TComplex bp, TComplex bm, double Pb,
               double alphaLambda, double cosTheta1, double cosTheta2,
               double phi1, double phi2, double cosTheta);
  double term6(TComplex ap, TComplex am, TComplex bp, TComplex bm, double Pb,
               double alphaLambda, double cosTheta1, double cosTheta2,
               double phi1, double phi2, double cosTheta);
  double term7(TComplex ap, TComplex am, TComplex bp, TComplex bm, double Pb,
               double alphaLambda, double cosTheta1, double cosTheta2,
               double phi1, double phi2, double cosTheta);
  double term8(TComplex ap, TComplex am, TComplex bp, TComplex bm, double Pb,
               double alphaLambda, double cosTheta1, double cosTheta2,
               double phi1, double phi2, double cosTheta);
  double term9(TComplex ap, TComplex am, TComplex bp, TComplex bm, double Pb,
               double alphaLambda, double cosTheta1, double cosTheta2,
               double phi1, double phi2, double cosTheta);
  double term10(TComplex ap, TComplex am, TComplex bp, TComplex bm, double Pb,
                double alphaLambda, double cosTheta1, double cosTheta2,
                double phi1, double phi2, double cosTheta);
  double term11(TComplex ap, TComplex am, TComplex bp, TComplex bm, double Pb,
                double alphaLambda, double cosTheta1, double cosTheta2,
                double phi1, double phi2, double cosTheta);
  double term12(TComplex ap, TComplex am, TComplex bp, TComplex bm, double Pb,
                double alphaLambda, double cosTheta1, double cosTheta2,
                double phi1, double phi2, double cosTheta);
  double term13(TComplex ap, TComplex am, TComplex bp, TComplex bm, double Pb,
                double alphaLambda, double cosTheta1, double cosTheta2,
                double phi1, double phi2, double cosTheta);
  double term14(TComplex ap, TComplex am, TComplex bp, TComplex bm, double Pb,
                double alphaLambda, double cosTheta1, double cosTheta2,
                double phi1, double phi2, double cosTheta);
  double term15(TComplex ap, TComplex am, TComplex bp, TComplex bm, double Pb,
                double alphaLambda, double cosTheta1, double cosTheta2,
                double phi1, double phi2, double cosTheta);
  double term16(TComplex ap, TComplex am, TComplex bp, TComplex bm, double Pb,
                double alphaLambda, double cosTheta1, double cosTheta2,
                double phi1, double phi2, double cosTheta);
  double term17(TComplex ap, TComplex am, TComplex bp, TComplex bm, double Pb,
                double alphaLambda, double cosTheta1, double cosTheta2,
                double phi1, double phi2, double cosTheta);
  double term18(TComplex ap, TComplex am, TComplex bp, TComplex bm, double Pb,
                double alphaLambda, double cosTheta1, double cosTheta2,
                double phi1, double phi2, double cosTheta);
  double term19(TComplex ap, TComplex am, TComplex bp, TComplex bm, double Pb,
                double alphaLambda, double cosTheta1, double cosTheta2,
                double phi1, double phi2, double cosTheta);

  TComplex aplus, aminus, bplus, bminus;
  double Pb, alphaLambda;
};

#endif
