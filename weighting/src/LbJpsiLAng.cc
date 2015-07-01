#include "LbJpsiLAng.hh"

#include <iostream>
#include "TMath.h"

double LbJpsiLAng::term0(TComplex ap, TComplex am, TComplex bp, TComplex bm,
                         double Pb, double alphaLambda, double cosThetaB,
                         double cosThetaL, double phiB, double phiL,
                         double cosTheta)
{
  return ap.Rho2() + am.Rho2() + bp.Rho2() + bm.Rho2();
}

double LbJpsiLAng::term1(TComplex ap, TComplex am, TComplex bp, TComplex bm,
                         double Pb, double alphaLambda, double cosThetaB,
                         double cosThetaL, double phiB, double phiL,
                         double cosTheta)
{
  double f1 = ap.Rho2() - am.Rho2() + bp.Rho2() - bm.Rho2();
  double f2 = Pb;
  return f1 * f2 * cosTheta;
}

double LbJpsiLAng::term2(TComplex ap, TComplex am, TComplex bp, TComplex bm,
                         double Pb, double alphaLambda, double cosThetaB,
                         double cosThetaL, double phiB, double phiL,
                         double cosTheta)
{
  double f1 = ap.Rho2() - am.Rho2() - bp.Rho2() + bm.Rho2();
  double f2 = alphaLambda;
  return f1 * f2 * cosThetaB;
}

double LbJpsiLAng::term3(TComplex ap, TComplex am, TComplex bp, TComplex bm,
                         double Pb, double alphaLambda, double cosThetaB,
                         double cosThetaL, double phiB, double phiL,
                         double cosTheta)
{
  double f1 = ap.Rho2() + am.Rho2() - bp.Rho2() - bm.Rho2();
  double f2 = alphaLambda * Pb;
  return f1 * f2 * cosTheta * cosThetaB;
}

double LbJpsiLAng::term4(TComplex ap, TComplex am, TComplex bp, TComplex bm,
                         double Pb, double alphaLambda, double cosThetaB,
                         double cosThetaL, double phiB, double phiL,
                         double cosTheta)
{
  double f1 = -ap.Rho2() - am.Rho2() + 0.5 * (bp.Rho2() + bm.Rho2());
  double f2 = 1;
  double wigner = 1.5 * cosThetaL * cosThetaL - 0.5;
  return f1 * f2 * wigner;
}

double LbJpsiLAng::term5(TComplex ap, TComplex am, TComplex bp, TComplex bm,
                         double Pb, double alphaLambda, double cosThetaB,
                         double cosThetaL, double phiB, double phiL,
                         double cosTheta)
{
  double f1 = -ap.Rho2() + am.Rho2() + 0.5 * (bp.Rho2() - bm.Rho2());
  double f2 = Pb;
  double wigner = 1.5 * cosThetaL * cosThetaL - 0.5;
  return f1 * f2 * wigner * cosTheta;
}

double LbJpsiLAng::term6(TComplex ap, TComplex am, TComplex bp, TComplex bm,
                         double Pb, double alphaLambda, double cosThetaB,
                         double cosThetaL, double phiB, double phiL,
                         double cosTheta)
{
  double f1 = -ap.Rho2() + am.Rho2() + 0.5 * (-bp.Rho2() + bm.Rho2());
  double f2 = alphaLambda;
  double wigner = 1.5 * cosThetaL * cosThetaL - 0.5;
  return f1 * f2 * wigner * cosThetaB;
}

double LbJpsiLAng::term7(TComplex ap, TComplex am, TComplex bp, TComplex bm,
                         double Pb, double alphaLambda, double cosThetaB,
                         double cosThetaL, double phiB, double phiL,
                         double cosTheta)
{
  double f1 = -ap.Rho2() - am.Rho2() + 0.5 * (-bp.Rho2() - bm.Rho2());
  double f2 = alphaLambda * Pb;
  double wigner = 1.5 * cosThetaL * cosThetaL - 0.5;
  return f1 * f2 * wigner * cosThetaB * cosTheta;
}

double LbJpsiLAng::term8(TComplex ap, TComplex am, TComplex bp, TComplex bm,
                         double Pb, double alphaLambda, double cosThetaB,
                         double cosThetaL, double phiB, double phiL,
                         double cosTheta)
{
  double f1 = -3.0 * (ap * TComplex::Conjugate(am)).Re();
  double f2 = alphaLambda * Pb;
  double angular = TMath::Sqrt(1 - cosTheta * cosTheta) *
                   TMath::Sqrt(1 - cosThetaB * cosThetaB) *
                   (1 - cosThetaL * cosThetaL) * TMath::Cos(phiB);
  return f1 * f2 * angular;
}

double LbJpsiLAng::term9(TComplex ap, TComplex am, TComplex bp, TComplex bm,
                         double Pb, double alphaLambda, double cosThetaB,
                         double cosThetaL, double phiB, double phiL,
                         double cosTheta)
{
  double f1 = 3.0 * (ap * TComplex::Conjugate(am)).Im();
  double f2 = alphaLambda * Pb;
  double angular = TMath::Sqrt(1 - cosTheta * cosTheta) *
                   TMath::Sqrt(1 - cosThetaB * cosThetaB) *
                   (1 - cosThetaL * cosThetaL) * TMath::Sin(phiB);
  return f1 * f2 * angular;
}

double LbJpsiLAng::term10(TComplex ap, TComplex am, TComplex bp, TComplex bm,
                          double Pb, double alphaLambda, double cosThetaB,
                          double cosThetaL, double phiB, double phiL,
                          double cosTheta)
{
  double f1 = -1.5 * (bm * TComplex::Conjugate(bp)).Re();
  double f2 = alphaLambda * Pb;
  double angular = TMath::Sqrt(1 - cosTheta * cosTheta) *
                   TMath::Sqrt(1 - cosThetaB * cosThetaB) *
                   (1 - cosThetaL * cosThetaL) * TMath::Cos(phiB + 2 * phiL);
  return f1 * f2 * angular;
}

double LbJpsiLAng::term11(TComplex ap, TComplex am, TComplex bp, TComplex bm,
                          double Pb, double alphaLambda, double cosThetaB,
                          double cosThetaL, double phiB, double phiL,
                          double cosTheta)
{
  double f1 = 1.5 * (bm * TComplex::Conjugate(bp)).Im();
  double f2 = alphaLambda * Pb;
  double angular = TMath::Sqrt(1 - cosTheta * cosTheta) *
                   TMath::Sqrt(1 - cosThetaB * cosThetaB) *
                   (1 - cosThetaL * cosThetaL) * TMath::Sin(phiB + 2 * phiL);
  return f1 * f2 * angular;
}

double LbJpsiLAng::term12(TComplex ap, TComplex am, TComplex bp, TComplex bm,
                          double Pb, double alphaLambda, double cosThetaB,
                          double cosThetaL, double phiB, double phiL,
                          double cosTheta)
{
  double f1 = -3.0 / TMath::Sqrt(2) * (bm * TComplex::Conjugate(ap) +
                                       am * TComplex::Conjugate(bp)).Re();
  double f2 = alphaLambda * Pb;
  double angular = TMath::Sqrt(1 - cosTheta * cosTheta) * cosThetaB *
                   TMath::Sqrt(1 - cosThetaL * cosThetaL) * cosThetaL *
                   TMath::Cos(phiL);
  return f1 * f2 * angular;
}

double LbJpsiLAng::term13(TComplex ap, TComplex am, TComplex bp, TComplex bm,
                          double Pb, double alphaLambda, double cosThetaB,
                          double cosThetaL, double phiB, double phiL,
                          double cosTheta)
{
  double f1 = 3.0 / TMath::Sqrt(2) * (bm * TComplex::Conjugate(ap) +
                                      am * TComplex::Conjugate(bp)).Im();
  double f2 = alphaLambda * Pb;
  double angular = TMath::Sqrt(1 - cosTheta * cosTheta) * cosThetaB *
                   TMath::Sqrt(1 - cosThetaL * cosThetaL) * cosThetaL *
                   TMath::Sin(phiL);
  return f1 * f2 * angular;
}

double LbJpsiLAng::term14(TComplex ap, TComplex am, TComplex bp, TComplex bm,
                          double Pb, double alphaLambda, double cosThetaB,
                          double cosThetaL, double phiB, double phiL,
                          double cosTheta)
{
  double f1 = -3.0 / TMath::Sqrt(2) * (bm * TComplex::Conjugate(am) +
                                       ap * TComplex::Conjugate(bp)).Re();
  double f2 = alphaLambda * Pb;
  double angular = cosTheta * TMath::Sqrt(1 - cosThetaB * cosThetaB) *
                   TMath::Sqrt(1 - cosThetaL * cosThetaL) * cosThetaL *
                   TMath::Cos(phiB + phiL);
  return f1 * f2 * angular;
}

double LbJpsiLAng::term15(TComplex ap, TComplex am, TComplex bp, TComplex bm,
                          double Pb, double alphaLambda, double cosThetaB,
                          double cosThetaL, double phiB, double phiL,
                          double cosTheta)
{
  double f1 = 3.0 / TMath::Sqrt(2) * (bm * TComplex::Conjugate(am) +
                                      ap * TComplex::Conjugate(bp)).Im();
  double f2 = alphaLambda * Pb;
  double angular = cosTheta * TMath::Sqrt(1 - cosThetaB * cosThetaB) *
                   TMath::Sqrt(1 - cosThetaL * cosThetaL) * cosThetaL *
                   TMath::Sin(phiB + phiL);
  return f1 * f2 * angular;
}

double LbJpsiLAng::term16(TComplex ap, TComplex am, TComplex bp, TComplex bm,
                          double Pb, double alphaLambda, double cosThetaB,
                          double cosThetaL, double phiB, double phiL,
                          double cosTheta)
{
  double f1 = +3.0 / TMath::Sqrt(2) * (am * TComplex::Conjugate(bp) -
                                       bm * TComplex::Conjugate(ap)).Re();
  double f2 = Pb;
  double angular = TMath::Sqrt(1 - cosTheta * cosTheta) *
                   TMath::Sqrt(1 - cosThetaL * cosThetaL) * cosThetaL *
                   TMath::Cos(phiL);
  return f1 * f2 * angular;
}

double LbJpsiLAng::term17(TComplex ap, TComplex am, TComplex bp, TComplex bm,
                          double Pb, double alphaLambda, double cosThetaB,
                          double cosThetaL, double phiB, double phiL,
                          double cosTheta)
{
  double f1 = -3.0 / TMath::Sqrt(2) * (am * TComplex::Conjugate(bp) -
                                       bm * TComplex::Conjugate(ap)).Im();
  double f2 = Pb;
  double angular = TMath::Sqrt(1 - cosTheta * cosTheta) *
                   TMath::Sqrt(1 - cosThetaL * cosThetaL) * cosThetaL *
                   TMath::Sin(phiL);
  return f1 * f2 * angular;
}

double LbJpsiLAng::term18(TComplex ap, TComplex am, TComplex bp, TComplex bm,
                          double Pb, double alphaLambda, double cosThetaB,
                          double cosThetaL, double phiB, double phiL,
                          double cosTheta)
{
  double f1 = 3.0 / TMath::Sqrt(2) * (bm * TComplex::Conjugate(am) -
                                      ap * TComplex::Conjugate(bp)).Re();
  double f2 = alphaLambda;
  double angular = TMath::Sqrt(1 - cosThetaB * cosThetaB) *
                   TMath::Sqrt(1 - cosThetaL * cosThetaL) * cosThetaL *
                   TMath::Cos(phiB + phiL);
  return f1 * f2 * angular;
}

double LbJpsiLAng::term19(TComplex ap, TComplex am, TComplex bp, TComplex bm,
                          double Pb, double alphaLambda, double cosThetaB,
                          double cosThetaL, double phiB, double phiL,
                          double cosTheta)
{
  double f1 = -3.0 / TMath::Sqrt(2) * (bm * TComplex::Conjugate(am) -
                                       ap * TComplex::Conjugate(bp)).Im();
  double f2 = alphaLambda;
  double angular = TMath::Sqrt(1 - cosThetaB * cosThetaB) *
                   TMath::Sqrt(1 - cosThetaL * cosThetaL) * cosThetaL *
                   TMath::Sin(phiB + phiL);
  return f1 * f2 * angular;
}

double LbJpsiLAng::physicsRate(double cosTheta, double cosThetaL,
                               double cosThetaB, double phiL, double phiB)
{
  double result = 0;

  double termResults[20];
  termResults[0] = term0(aplus, aminus, bplus, bminus, Pb, alphaLambda,
                         cosThetaB, cosThetaL, phiB, phiL, cosTheta);
  termResults[1] = term1(aplus, aminus, bplus, bminus, Pb, alphaLambda,
                         cosThetaB, cosThetaL, phiB, phiL, cosTheta);
  termResults[2] = term2(aplus, aminus, bplus, bminus, Pb, alphaLambda,
                         cosThetaB, cosThetaL, phiB, phiL, cosTheta);
  termResults[3] = term3(aplus, aminus, bplus, bminus, Pb, alphaLambda,
                         cosThetaB, cosThetaL, phiB, phiL, cosTheta);
  termResults[4] = term4(aplus, aminus, bplus, bminus, Pb, alphaLambda,
                         cosThetaB, cosThetaL, phiB, phiL, cosTheta);
  termResults[5] = term5(aplus, aminus, bplus, bminus, Pb, alphaLambda,
                         cosThetaB, cosThetaL, phiB, phiL, cosTheta);
  termResults[6] = term6(aplus, aminus, bplus, bminus, Pb, alphaLambda,
                         cosThetaB, cosThetaL, phiB, phiL, cosTheta);
  termResults[7] = term7(aplus, aminus, bplus, bminus, Pb, alphaLambda,
                         cosThetaB, cosThetaL, phiB, phiL, cosTheta);
  termResults[8] = term8(aplus, aminus, bplus, bminus, Pb, alphaLambda,
                         cosThetaB, cosThetaL, phiB, phiL, cosTheta);
  termResults[9] = term9(aplus, aminus, bplus, bminus, Pb, alphaLambda,
                         cosThetaB, cosThetaL, phiB, phiL, cosTheta);
  termResults[10] = term10(aplus, aminus, bplus, bminus, Pb, alphaLambda,
                           cosThetaB, cosThetaL, phiB, phiL, cosTheta);
  termResults[11] = term11(aplus, aminus, bplus, bminus, Pb, alphaLambda,
                           cosThetaB, cosThetaL, phiB, phiL, cosTheta);
  termResults[12] = term12(aplus, aminus, bplus, bminus, Pb, alphaLambda,
                           cosThetaB, cosThetaL, phiB, phiL, cosTheta);
  termResults[13] = term13(aplus, aminus, bplus, bminus, Pb, alphaLambda,
                           cosThetaB, cosThetaL, phiB, phiL, cosTheta);
  termResults[14] = term14(aplus, aminus, bplus, bminus, Pb, alphaLambda,
                           cosThetaB, cosThetaL, phiB, phiL, cosTheta);
  termResults[15] = term15(aplus, aminus, bplus, bminus, Pb, alphaLambda,
                           cosThetaB, cosThetaL, phiB, phiL, cosTheta);
  termResults[16] = term16(aplus, aminus, bplus, bminus, Pb, alphaLambda,
                           cosThetaB, cosThetaL, phiB, phiL, cosTheta);
  termResults[17] = term17(aplus, aminus, bplus, bminus, Pb, alphaLambda,
                           cosThetaB, cosThetaL, phiB, phiL, cosTheta);
  termResults[18] = term18(aplus, aminus, bplus, bminus, Pb, alphaLambda,
                           cosThetaB, cosThetaL, phiB, phiL, cosTheta);
  termResults[19] = term19(aplus, aminus, bplus, bminus, Pb, alphaLambda,
                           cosThetaB, cosThetaL, phiB, phiL, cosTheta);

  for (int i = 0; i <= 19; ++i)
  {
    result += termResults[i];
  }

  if (result < -1e-5 || TMath::IsNaN(result))
  {
    std::cout << "Negative rate: \n";
    for (int i = 0; i <= 19; ++i)
      std::cout << i << " " << termResults[i] << std::endl;

    std::cout << "-----------------------------\n";
    std::cout << aplus << " " << aminus << " " << bplus << " " << bminus << " "
              << Pb << " " << alphaLambda << std::endl;
    ;
    std::cout << "=============================\n";
  }

  return result;
}

void LbJpsiLAng::setAlphabR0R1(double alphaB, double R0, double R1,
                               double apPhase, double amPhase, double bpPhase,
                               double bmPhase)
{
  double apSq = 0.5 * (R0 + R1);
  double amSq = 0.5 * (R0 - R1);
  double bpSq = 0.5 * (1 + alphaB - R0 - R1);
  double bmSq = 0.5 * (1 - alphaB - R0 + R1);

  if (apSq < 0 || amSq < 0 || bpSq < 0 || bmSq < 0)
  {
    std::cout << "Amplitudes are unphysical: " << apSq << " " << amSq << " "
              << bpSq << " " << bmSq << std::endl;
  }

  if (apSq > -0.2 &&
      apSq < 0) // Put it to zero if we are negative and close to zero
  {
    apSq = 0;
  }
  if (amSq > -0.2 &&
      amSq < 0) // Put it to zero if we are negative and close to zero
  {
    amSq = 0;
  }
  if (bpSq > -0.2 &&
      bpSq < 0) // Put it to zero if we are negative and close to zero
  {
    bpSq = 0;
  }
  if (bmSq > -0.2 &&
      bmSq < 0) // Put it to zero if we are negative and close to zero
  {
    bmSq = 0;
  }

  // Want to renormalizae amplitudes to sum=1 in case we moved some of them
  // to zero
  // This is not really needed as things are normalized away by weighting
  // denominator for efficiency

  aplus = TComplex(TMath::Sqrt(apSq) * TMath::Cos(apPhase),
                   TMath::Sqrt(apSq) * TMath::Sin(apPhase));
  aminus = TComplex(TMath::Sqrt(amSq) * TMath::Cos(amPhase),
                    TMath::Sqrt(amSq) * TMath::Sin(amPhase));
  bplus = TComplex(TMath::Sqrt(bpSq) * TMath::Cos(bpPhase),
                   TMath::Sqrt(bpSq) * TMath::Sin(bpPhase));
  bminus = TComplex(TMath::Sqrt(bmSq) * TMath::Cos(bmPhase),
                    TMath::Sqrt(bmSq) * TMath::Sin(bmPhase));

  std::cout << "Parameters setup: \n";
  std::cout << aplus << " " << aminus << " " << bplus << " " << bminus
            << std::endl;
}

LbJpsiLAng::LbJpsiLAng()
{
  setPolarization(0.06);
  setLambdaAsymmetry(0.643);
  setAlphabR0R1(0.04, 0.58, -0.56);
}

LbJpsiLAng::LbJpsiLAng(double polarization, double as, double alphaB, double R0,
                       double R1, double apPhase, double amPhase,
                       double bpPhase, double bmPhase)
{
  setPolarization(polarization);
  setLambdaAsymmetry(as);
  setAlphabR0R1(alphaB, R0, R1, apPhase, amPhase, bpPhase, bmPhase);
}

/*
 * Calculation of all relevant angles for Lb decays
 */
void LbJpsiLAng::LbPsiRAngles(TLorentzVector initialProton, TLorentzVector pB,
                              TLorentzVector pJpsi, TLorentzVector pLambda,
                              TLorentzVector pmp, TLorentzVector pmm,
                              TLorentzVector pp, TLorentzVector ppi,
                              int pCharge, double &cosTheta, double &cosThetaL,
                              double &cosThetaB, double &phiL, double &phiB,
                              double &dphi)
{

  bool decayLambda = true;
  if (pp.M() <
      1e-5) // Should be zero exactly, but allow for numerical treatment
  {
    decayLambda = false;
  }
  TVector3 analyzer;
  analyzer = initialProton.Vect().Cross(pB.Vect());
  analyzer *= 1.0 / analyzer.Mag();

  // Move everything to Lambda_b rest frame
  if (pB.Vect().Mag() > 1e-6)
  {
    TVector3 boost = -1.0 * pB.BoostVector();
    pB.Boost(boost);
    pJpsi.Boost(boost);
    pLambda.Boost(boost);
    pmp.Boost(boost);
    pmm.Boost(boost);
    pp.Boost(boost);
    ppi.Boost(boost);
  }
  // First level
  cosTheta =
      analyzer.Dot(pLambda.Vect()) / analyzer.Mag() / pLambda.Vect().Mag();

  // Frame 1
  TVector3 z1 = pLambda.Vect() * (1. / pLambda.Vect().Mag());
  TVector3 y1 = analyzer.Cross(z1);
  y1 = y1 * (1. / y1.Mag());
  TVector3 x1 = z1.Cross(y1);
  x1 = x1 * (1. / x1.Mag());
  // Proton in Lambda rest frame
  TLorentzVector ppLRest = pp;
  TVector3 h1;
  if (decayLambda)
  {
    ppLRest.Boost(-1.0 * pLambda.BoostVector());
    cosThetaB = z1.Dot(ppLRest.Vect()) / ppLRest.Vect().Mag();
    TVector3 pperp =
        ppLRest.Vect() - 1.0 * ppLRest.Vect().Mag() * z1 * cosThetaB;
    double cosPhi = pperp.Dot(x1) / pperp.Mag();
    double sinPhi = pperp.Dot(y1) / pperp.Mag();
    phiB = TMath::ACos(cosPhi);
    if (sinPhi < 0)
      phiB = -phiB;
    //   phiB=2*TMath::Pi()-phiB;

    h1=((ppLRest.Vect().Cross(pLambda.Vect())));
  }
  // Frame 2
  TVector3 z2 = pJpsi.Vect() * (1. / pJpsi.Vect().Mag());
  TVector3 y2 = analyzer.Cross(z2);
  y2 = y2 * (1. / y2.Mag());
  TVector3 x2 = z2.Cross(y2);
  x2 = x2 * (1. / x2.Mag());
  // Proton in Lambda rest frame
  ppLRest = pmp;
  ppLRest.Boost(-1.0 * pJpsi.BoostVector());
  cosThetaL = z2.Dot(ppLRest.Vect()) / ppLRest.Vect().Mag();
  TVector3 pperp = ppLRest.Vect() - 1.0 * ppLRest.Vect().Mag() * z2 * cosThetaL;
  double cosPhi = pperp.Dot(x2) / pperp.Mag();
  double sinPhi = pperp.Dot(y2) / pperp.Mag();
  phiL = TMath::ACos(cosPhi);
  if (sinPhi < 0)
    phiL = -phiL;

  if (decayLambda)
  {
    TVector3 h2((ppLRest.Vect().Cross(pJpsi.Vect())));

    // Calculate delta phi
    //  TVector3 h1((pmp.Vect()).Cross(pmm.Vect()));
    //  TVector3 h2((pp.Vect()).Cross(ppi.Vect()));
    TVector3 dir(h1.Cross(h2));
    float a1 = dir.Dot(pJpsi.Vect()) / dir.Mag() / pJpsi.Vect().Mag();
    a1 = TMath::ACos(a1);
    dphi = h1.Dot(h2) / h1.Mag() / h2.Mag();
    //  infoToStore.dphi=TMath::Pi()+TMath::ACos(infoToStore.dphi);
    dphi = TMath::ACos(dphi);
    if (a1 > TMath::Pi() - 0.001)
    {
      dphi = -dphi;
    }
  }
  // Now we handle antiparticles (charge of the proton being negative)

  if (pCharge < 0)
  {
    cosThetaL *= -1;
    if (phiB > 0)
    {
      phiB = TMath::Pi() - phiB;
    }
    else
    {
      phiB = -TMath::Pi() - phiB;
    }
    if (phiL > 0)
    {
      phiL = TMath::Pi() - phiL;
    }
    else
    {
      phiL = -TMath::Pi() - phiL;
    }
    if (dphi > 0)
    {
      dphi = TMath::Pi() - dphi;
    }
    else
    {
      dphi = -TMath::Pi() - dphi;
    }
  }
}

void LbJpsiLAng::setVariation(int variation)
{

  switch (variation)
  {
  case 0:
    setPolarization(0.06);
    setAlphabR0R1(0.04, 0.58, -0.56);
    break;
  case 1:
    setPolarization(0.068185);
    setAlphabR0R1(0.221468, 0.568174, -0.462248);
    break;
  case 2:
    setPolarization(0.051815);
    setAlphabR0R1(-0.141468, 0.591827, -0.657751);
    break;
  case 3:
    setPolarization(0.131683);
    setAlphabR0R1(0.0320301, 0.580041, -0.551201);
    break;
  case 4:
    setPolarization(-0.0116826);
    setAlphabR0R1(0.0479705, 0.57996, -0.568798);
    break;
  case 5:
    setPolarization(0.0502729);
    setAlphabR0R1(0.0116308, 0.580536, -0.506455);
    break;
  case 6:
    setPolarization(0.069727);
    setAlphabR0R1(0.0683698, 0.579465, -0.613545);
    break;
  case 7:
    setPolarization(0.0600591);
    setAlphabR0R1(0.041035, 0.59897, -0.559631);
    break;
  case 8:
    setPolarization(0.0599409);
    setAlphabR0R1(0.0389655, 0.561031, -0.560369);
    break;
  // Hand crafted shift without taking correlation into account, which moves
  // all magnitudes to be physical. It does 1sigma up/down for R0/R1 and
  // small fraction of uncertainty on alphaB
  case 9:
    setPolarization(0.06);
    setAlphabR0R1(-0.06, 0.6, -0.46);
    break;
  default:
    std::cout << "Using default settings\n";
    setPolarization(0.06);
    setAlphabR0R1(0.04, 0.58, -0.56);
    break;
  }
}
