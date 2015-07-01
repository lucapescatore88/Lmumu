#include "LbJpsiL2.hh"
#include "WignerD.hh"

#include "TMath.h"
#include "TComplex.h"
#include <iostream>

LbJpsiL2::LbJpsiL2()
{
  Pb = 0.06;
  alphaLambda =
      0.643; // This will really need check that all comes right for Lbbar.

  setAlphabR0R1(0.04, 0.58, -0.56);
}

/*
 * Function which calculates rate. Version which decays Lambda --> p pi
 */
double LbJpsiL2::physicsRate(double cosTheta, double cosThetaL,
                             double cosThetaB, double phiL, double phiB)
{
  double result = 0;

  updateSpinDensity(cosTheta);

  // While result is real, we do sum with complex numbers to have control that
  // we really get real number
  TComplex tmp = 0;
  int J = 1;
  int Jp = 1;
  // Nice nasty nested sum
  for (double ll = -1; ll <= 1.0; ll += 2)
  {
    for (int lj = -J; lj <= J; ++lj)
    {
      for (int ljp = -Jp; ljp <= Jp; ++ljp)
      {
        for (double lL = -0.5; lL <= 0.5; ++lL)
        {
          for (double lLp = -0.5; lLp <= 0.5; ++lLp)
          {
            for (double lp = -0.5; lp <= 0.5; ++lp)
            {
              tmp +=
                  hLepton(J, ll) * hLepton(Jp, ll) *
                  spinDensityMatrix(lj - lL, ljp - lLp) *
                  TMath::Power(-1, J + Jp) * hLb(lL, lj) *
                  TComplex::Conjugate(hLb(lLp, ljp)) * hLambda(lp) *
                  hLambda(lp) *
                  WignerD::WignerDfunction(J, lj, ll, cosThetaL, phiL) *
                  TComplex::Conjugate(
                      WignerD::WignerDfunction(Jp, ljp, ll, cosThetaL, phiL)) *
                  WignerD::WignerDfunction(0.5, lL, lp, cosThetaB, phiB) *
                  TComplex::Conjugate(
                      WignerD::WignerDfunction(0.5, lLp, lp, cosThetaB, phiB));
            }
          }
        }
      }
    }
  }

  if (TMath::Abs(tmp.Im()) > 1e-5)
  {
    std::cout << "Non-zero imaginary part, " << tmp.Im() << std::endl;
  }
  result = 3 * tmp.Re();

  if (TMath::IsNaN(result))
  {
	result = -999;
    //std::cout << "NAN detected: " << cosTheta << " " << cosThetaL << " "
    //          << cosThetaB << " " << phiL << " " << phiB << std::endl;
  }

  return result;
}

/*
 * Function which calculates rate. Version which treats Lambda as stable
 * final state particle.
 */
double LbJpsiL2::physicsRate(double cosTheta, double cosThetaL, double phiL)
{
  double result = 0;

  updateSpinDensity(cosTheta);

  // While result is real, we do sum with complex numbers to have control that
  // we really get real number
  TComplex tmp = 0;
  int J = 1;
  int Jp = 1;
  // Nice nasty nested sum
  for (double ll = -1; ll <= 1.0; ll += 2)
  {
    for (int lj = -J; lj <= J; ++lj)
    {
      for (int ljp = -Jp; ljp <= Jp; ++ljp)
      {
        for (double lL = -0.5; lL <= 0.5; ++lL)
        {
          tmp += hLepton(J, ll) * hLepton(Jp, ll) *
                 spinDensityMatrix(-lj - lL, -ljp - lL) *
                 TMath::Power(-1, J + Jp) * hLb(lL, -lj) *
                 TComplex::Conjugate(hLb(lL, -ljp)) *
                 WignerD::WignerDfunction(J, -lj, ll, cosThetaL, phiL) *
                 TComplex::Conjugate(
                     WignerD::WignerDfunction(Jp, -ljp, ll, cosThetaL, phiL));
        }
      }
    }
  }

  if (TMath::Abs(tmp.Im()) > 1e-5)
  {
    std::cout << "Non-zero imaginary part, " << tmp.Im() << std::endl;
  }
  result = 3. * tmp.Re();

  return result;
}

double LbJpsiL2::hLepton(int J, double ll)
{
  if (TMath::Abs(ll) > J)
  {
    return 0;
  }

  if (TMath::Abs(ll) == 1)
  {
    return 1;
  }
  return 0;
}

double LbJpsiL2::spinDensityMatrix(double lLb, double lLbp)
{
  int i = 0;
  int j = 0;
  if (lLb < 0)
    i = 1;
  if (lLbp < 0)
    j = 1;

  return spinDensity[i][j];
}

TComplex LbJpsiL2::hLb(double lL, int lj)
{
  // Here we need to return one of the four complex amplitudes

  if (lL > 0 && lj > 0)
    return bminus;
  else if (lL > 0 && lj == 0)
    return aplus;
  else if (lL < 0 && lj == 0)
    return aminus;
  else if (lL < 0 && lj < 0)
    return bplus;

  return 0;
}

void LbJpsiL2::setAlphabR0R1(double alphaB, double R0, double R1,
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

double LbJpsiL2::hLambda(double lp)
{
  if (lp > 0)
  {
    return TMath::Sqrt(1 + alphaLambda);
  }
  else
  {
    return TMath::Sqrt(1 - alphaLambda);
  }
}

void LbJpsiL2::updateSpinDensity(double cosTheta)
{
  spinDensity[0][0] = (1. - Pb * cosTheta) / 2.;
  spinDensity[1][1] = (1. + Pb * cosTheta) / 2.;
  spinDensity[0][1] = 0.5 * Pb * TMath::Sqrt(1 - cosTheta * cosTheta);
  spinDensity[1][0] = 0.5 * Pb * TMath::Sqrt(1 - cosTheta * cosTheta);
}

void LbJpsiL2::setVariation(int variation)
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

/*
 * Calculation of all relevant angles for Lb decays
 */
void LbJpsiL2::LbPsiRAngles(TLorentzVector initialProton, TLorentzVector pB,
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
