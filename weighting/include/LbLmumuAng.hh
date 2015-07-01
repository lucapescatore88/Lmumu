#ifndef LBLMUMUANG_HH
#define LBLMUMUANG_HH

#include "TComplex.h"
#include "TLorentzVector.h"

/* 
 * Code has many constants we might in principle want to update hardcoded,
 * but for now I worry about having all formulas correct rather than
 * changing numerical values somewhere.
 */

class LbLmumuAng
{

  public:

    LbLmumuAng( double pol = 0, int wilson = 0 );

    ~LbLmumuAng() {};

    double physicsRate(double qsq, double cosTheta, double cosTheta1, double cosThetaB, double phiL, double phiB);
    double physicsRate(double qsq, double cosTheta, double cosThetaL, double phiL);

    void setPolarization(double polarization) {Pb=polarization;};
    void setLambdaAsymmetry(double as) {alphaLambda=as;};

    void setLeptonMass(double ml) { mlepton=ml; }
    void setAlphaLambda(double al) { alphaLambda=al;}

    void setq2(double ff) {fq2=ff;}

    typedef enum
    {
      CovQuarkModel, LQCD 
    } FFtype;

    void setFFtype(int ff)
    {
      formFactors = ff;
    }
//  private:

    double Pb, alphaLambda;
    double mlepton;
    double v;
    double fq2;
    double Qplus, Qminus;
    double Mplus, Mminus;
    double mLb, mL;  // M1 and M2 in Gutsche paper
    TComplex C[11], C9eff, CNP[11];
    int formFactors;

    double hLepton(int J, int m, double l1, double l2);
    double spinDensityMatrix(double lLb, double lLbp);
    TComplex hLb(int J, int m, double lL, int lj);
    TComplex hLbV(int J, int m, double lL, int lj);
    TComplex hLbA(int J, int m, double lL, int lj);
    double hLambda(double lp);
    TComplex FV(int i, int m);
    TComplex FA(int i, int m);
    double fV(int i);
    double fTV(int i);
    double fA(int i);
    double fTA(int i);
    double formFactorParametrization(double s, double f0, double a, double b);
    double formFactorParametrizationLQCD(double var, double N, double X);
    double fVconsts[3][3]; 
    double fAconsts[3][3]; 
    double fTVconsts[3][3]; 
    double fTAconsts[3][3]; 
    double FiLQCDconsts[3][2];
    double alphaSmZ;
    double alphaS;
    double mZ;
    double cgamma, cv, csigma;
    void calculateCurrentMatchingConstants();

    double spinDensity[2][2];
    void updateSpinDensity(double cosTheta);
    void updateKinematics(double q2);
    void updateC9eff(double q2);
    TComplex h(double m, double s);
    double mc, mb, mu, mJpsi, mpsi, gJpsi, gpsi, gJpsill, gpsill;
    double alpha;

};


#endif
