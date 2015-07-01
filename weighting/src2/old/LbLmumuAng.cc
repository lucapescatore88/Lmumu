#include "LbLmumuAng.hh"
#include "WignerD.hh"

#include "TMath.h"
#include "TComplex.h"
#include <iostream>

LbLmumuAng::LbLmumuAng(double polarization)
{
  mb=4.19;
  mu=mb; // This is what is done in Gutsche paper, mu is renormalization scale for Wilson coeeficients.
  mc=1.27;
  mLb=5.6194;
  mlepton=0.1056583715; // muons
  Pb=polarization;  // Do unpolarized for now, once we are sure about implementation, change this
  alphaLambda=0.642;  // This will really need check that all comes right for Lbbar.
  mL=1.115683;
  mJpsi=3.096916;
  mpsi=3.68609;
  gJpsi=0.0929;
  gpsi=0.304;
  gJpsill=0.00555;  // Putting gJpsill and gpsill to zero should remove those big spikes.
  gpsill=0.00235;
  gJpsill=0;
  gpsill=0;
  alpha=1/129.; // I guess this is fine structure constant at who knows which scale
  // Willson coeeficients
  C[1]=-0.248;
  C[2]= 1.107;
  C[3]= 0.011;
  C[4]= -0.026;
  C[5]=  0.007;
  C[6]= -0.031;
  C[7]= -0.313;
  C[8]= 0;
  C[9]= 4.344; 
  C[10]= -4.669;
  C[0]=3.0*C[1]+C[2]+3.0*C[3]+C[4]+3.0*C[5]+C[6];
  // Form factors  
  fVconsts[0][0]=0.107;
  fVconsts[0][1]=2.27;
  fVconsts[0][2]=1.367;
  fVconsts[1][0]=0.043;
  fVconsts[1][1]=2.411;
  fVconsts[1][2]=1.531;
  fVconsts[2][0]=0.003;
  fVconsts[2][1]=2.815;
  fVconsts[2][2]=2.041;

  fAconsts[0][0]=0.104;
  fAconsts[0][1]=2.232;
  fAconsts[0][2]=1.328;
  fAconsts[1][0]=0.003;
  fAconsts[1][1]=2.955;
  fAconsts[1][2]=3.620;
  fAconsts[2][0]=-0.052;
  fAconsts[2][1]=2.437;
  fAconsts[2][2]=1.559;

  fTVconsts[0][0]=-0.043;
  fTVconsts[0][1]=2.411;
  fTVconsts[0][2]=1.531;
  fTVconsts[1][0]=-0.105;
  fTVconsts[1][1]=0.072;
  fTVconsts[1][2]=0.001;
  fTVconsts[2][0]=0;  // Not used anywhere
  fTVconsts[2][1]=0;
  fTVconsts[2][2]=0;

  fTAconsts[0][0]=0.003;
  fTAconsts[0][1]=2.955;
  fTAconsts[0][2]=3.620;
  fTAconsts[1][0]=-0.105;
  fTAconsts[1][1]=2.233;
  fTAconsts[1][2]=1.328;
  fTAconsts[2][0]=0;  // Not used anywhere
  fTAconsts[2][1]=0;
  fTAconsts[2][2]=0;
}

double LbLmumuAng::physicsRate(double q2, double cosTheta, double cosThetaL, double cosThetaB, double phiL, double phiB)
{
  double result=0;

  updateSpinDensity(cosTheta);
  updateKinematics(q2);
  updateC9eff(q2);

  // While result is real, we do sum with complex numbers to have control that we really get real number
  TComplex tmp=0;
  // Nice nasty nested sum
  for (int J=0; J<=1; ++J)
  {
    for (int Jp=0; Jp<=1; ++Jp)
    {
      for (double l1=-0.5; l1<=0.5; ++l1)
      {
        for (double l2=-0.5; l2<=0.5; ++l2)
        {
          for (int lj=-J;lj<=J;++lj)
          {
            for (int ljp=-Jp;ljp<=Jp;++ljp)
            {
              for (int m=1; m<=2; ++m)
              {
                for (int mp=1; mp<=2; ++mp)
                {
                  for (double lL=-0.5; lL<=0.5; ++lL)
                  {
                    for (double lLp=-0.5; lLp<=0.5; ++lLp)
                    {
                      for (double lp=-0.5; lp<=0.5; ++lp)
                      {
//   std::cout<<"Adding "<<J<<" "<<Jp<<std::endl;
                        tmp+=hLepton(J,m,l1,l2)*hLepton(Jp,mp,l1,l2)*
                             spinDensityMatrix(lj-lL,ljp-lLp)*
                             TMath::Power(-1,J+Jp)*
                             hLb(J,m,lL,lj)*TComplex::Conjugate(hLb(Jp,mp,lLp,ljp))*
                             hLambda(lp)*hLambda(lp)*  
                             WignerD::WignerDfunction(J,lj,l1-l2,cosThetaL,phiL)*
                             TComplex::Conjugate(WignerD::WignerDfunction(Jp,ljp,l1-l2,cosThetaL,phiL))*
                             WignerD::WignerDfunction(0.5,lL,lp,cosThetaB,phiB)*
                             TComplex::Conjugate(WignerD::WignerDfunction(0.5,lLp,lp,cosThetaB,phiB));
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  if ( TMath::Abs(tmp.Im())>1e-5 )
  {
    std::cout<<"Non-zero imaginary part, "<<tmp.Im()<<std::endl;
  }
  result=tmp.Re();

  if ( TMath::IsNaN(result) )
  {
    //std::cout<<"NAN detected: "<<q2<<" "<<cosTheta<<" "<<cosThetaL<<" "<<cosThetaB<<" "<<phiL<<" "<<phiB<<std::endl;
    result = -999;
  }

  return result;
}

double LbLmumuAng::hLepton(int J, int m, double l1, double l2)
{ 
  if ( TMath::Abs(l1-l2)>J )
  {
    return 0;
  }
  if ( l1<0 )
  {
    if ( m==1 )
    {
      return hLepton(J,m,-l1,-l2);
    }
    else
    {
      return -hLepton(J,m,-l1,-l2);
    }
  }

  if ( J==0 && m==1 && l1>0 && l2>0 )
  {
    return 0;
  }
  if ( J==0 && m==2 && l1>0 && l2>0 )
  {
    return 2*mlepton;
  }
  if ( J==1 )
  {
    if ( m==1 )
    {
      if ( l1>0 && l2>0 )
      {
        return 2*mlepton;
      }
      if ( l1>0 && l2<0 )
      {
        return -TMath::Sqrt(2*fq2);
      }
    }
    if ( m==2 )
    {
      if ( l1>0 && l2>0 )
      {
        return 0;
      }
      if ( l1>0 && l2<0 )
      {
        return TMath::Sqrt(2*fq2)*v;
      }
    }
  }

  std::cout<<"Impossible inputs in hLepton: "<<
     J<<" "<<m<<" "<<l1<<" "<<l2<<std::endl;

  return -999;
}
double LbLmumuAng::spinDensityMatrix(double lLb, double lLbp)
{ 
  int i=0;
  int j=0;
  if ( lLb<0 )
    i=1;
  if ( lLbp<0 )
    j=1;

  return spinDensity[i][j];
}

TComplex LbLmumuAng::hLb(int J, int m, double lL, int lj)
{ 
  return hLbV(J, m, lL, lj)-hLbA(J, m, lL, lj);
}

TComplex LbLmumuAng::hLbV(int J, int m, double lL, int lj)
{ 
  if ( J==0 && lj!=0 )
    return 0;
  if ( J==0 && lj==0 )
  {
    // For this one, lL = +1/2 and -1/2 are same
    return TMath::Sqrt(Qplus/fq2)*(Mminus*FV(1,m)+fq2/mLb*FV(3,m));
  }
  if ( J==1 )
  {
    if ( lL>0 )
    {
      switch (lj)
      {
        case 1: return TMath::Sqrt(2*Qminus)*(FV(1,m)+Mplus/mLb*FV(2,m));
                break;
        case 0: return TMath::Sqrt(Qminus/fq2)*(Mplus*FV(1,m)+fq2/mLb*FV(2,m));
                break;
        case -1: return 0;
                break;
      }
    }
    else
    {
      return hLbV(J,m,-lL,-lj);
    }
  }

  return -999;
}

TComplex LbLmumuAng::hLbA(int J, int m, double lL, int lj)
{ 
//  std::cout<<"Entering hLbA with "<<J<<" "<<m<<" "<<lL<<" "<<lj<<std::endl;
  if ( J==0 && lj!=0 )
    return 0;
  if ( J==0 && lj==0 )
  {
    if ( lL>0 )
    {
      return TMath::Sqrt(Qminus/fq2)*(Mplus*FA(1,m)-fq2/mLb*FA(3,m));
    }
    else
    {
      return -hLbA(J, m, -lL, -lj);
    }
  }
  if ( J==1 )
  {
    if ( lL>0 )
    {
      switch (lj)
      {
        case 1: return TMath::Sqrt(2*Qplus)*(FA(1,m)-Mminus/mLb*FA(2,m));
                break;
        case 0: return TMath::Sqrt(Qplus/fq2)*(Mminus*FA(1,m)-fq2/mLb*FA(2,m));
                break;
        case -1: return 0;
                 break;
      }
    }
    else
    {
//      std::cout<<"hLbA with "<<lL<<" "<<lj<<" transformed to "<<-lL<<" "<<-lj<<std::endl;
      return -hLbA(J,m,-lL,-lj);
    }
  }

  return -999;
}

TComplex LbLmumuAng::FV(int i, int m)
{
  if ( m==2 )
  {
    return C[10]*fV(m);
  }
  if ( m==1 )
  {
    switch (i)
    {
      case 1: return C9eff*fV(1)-2*mb/mLb*C[7]*fTV(1);
              break;
      case 2: return C9eff*fV(2)-2*mb*mLb/fq2*C[7]*fTV(2);
              break;
      case 3: return C9eff*fV(3)+2*mb*Mminus/fq2*C[7]*fTV(1);
              break;
    }
  }
  return -999;
}

TComplex LbLmumuAng::FA(int i, int m)
{
  if ( m==2 )
  {
    return C[10]*fA(m);
  }
  if ( m==1 )
  {
    switch (i)
    {
      case 1: return C9eff*fA(1)+2*mb/mLb*C[7]*fTA(1);
              break;
      case 2: return C9eff*fA(2)+2*mb*mLb/fq2*C[7]*fTA(2);
              break;
      case 3: return C9eff*fA(3)+2*mb*Mplus/fq2*C[7]*fTA(1);
              break;
    }
  }
  return -999;
}

double LbLmumuAng::formFactorParametrization(double s, double f0, double a, double b)
{
  return f0/(1-a*s+b*s*s);
}

double LbLmumuAng::fV(int i)
{
  return formFactorParametrization(fq2/mLb/mLb, fVconsts[i-1][0], fVconsts[i-1][1], fVconsts[i-1][2]);
}

double LbLmumuAng::fTV(int i)
{
  return formFactorParametrization(fq2/mLb/mLb, fTVconsts[i-1][0], fTVconsts[i-1][1], fTVconsts[i-1][2]);
}

double LbLmumuAng::fA(int i)
{
  return formFactorParametrization(fq2/mLb/mLb, fAconsts[i-1][0], fAconsts[i-1][1], fAconsts[i-1][2]);
}

double LbLmumuAng::fTA(int i)
{
  return formFactorParametrization(fq2/mLb/mLb, fTAconsts[i-1][0], fTAconsts[i-1][1], fTAconsts[i-1][2]);
}

void LbLmumuAng::updateC9eff(double q2)
{
  C9eff=C[9]+C[0]*(h(mc/mLb,q2/mLb/mLb)+
         3*TMath::Pi()/alpha/alpha*1/C[0]*
           (gJpsill*mJpsi/(mJpsi*mJpsi-q2-TComplex::I()*mJpsi*gJpsi)+
            gpsill*mpsi/(mpsi*mpsi-fq2-TComplex::I()*mpsi*gpsi)))-
         h(1,q2/mLb/mLb)*(4.0*C[3]+4.0*C[4]+3.0*C[5]+C[6])*0.5-
         h(0,q2/mLb/mLb)*(C[3]+3.0*C[4])*0.5+
         2./9*(3.0*C[3]+C[4]+3.0*C[5]+C[6]);
}

TComplex LbLmumuAng::h(double m, double s)
{
  if ( m==0 )
  {
    return (8./27-8./9*TMath::Log(mb/mu)-4./9*TMath::Log(s)+4./9*TComplex::I()*TMath::Pi());
  }

  TComplex A;
  double x=4*mc*mc/s;
  if ( x > 1 )
  {
    A=2*TMath::ATan(1/TMath::Sqrt(x-1));
  }
  else
  {
    A=-1.0*TComplex::I()*TMath::Pi()+
      TMath::Log(TMath::Abs((TMath::Sqrt(1-x)+1)/(TMath::Sqrt(1-x)-1)));
  }

  return -8./9*TMath::Log(mb/mu)-8./9*TMath::Log(mc)+8./27+4./9*x-
         2./9*(2+x)*TMath::Sqrt(TMath::Abs(1-x))*A;
}

double LbLmumuAng::hLambda(double lp)
{
  if ( lp>0 )
  {
    return 1+alphaLambda;
  }
  else 
  {
    return 1-alphaLambda;
  }
}

void LbLmumuAng::updateSpinDensity(double cosTheta)
{
  spinDensity[0][0]=(1.+Pb*cosTheta)/2.;
  spinDensity[1][1]=(1.-Pb*cosTheta)/2.;
  spinDensity[0][1]=Pb*TMath::Sqrt(1-cosTheta*cosTheta);
  spinDensity[1][0]=Pb*TMath::Sqrt(1-cosTheta*cosTheta);
}

void LbLmumuAng::updateKinematics(double q2)
{
  fq2=q2;
  v=TMath::Sqrt(1-4*mlepton*mlepton/q2); 
  Mplus=mLb+mL;
  Mminus=mLb-mL;
  Qplus=Mplus*Mplus-q2;
  Qminus=Mminus*Mminus-q2;
}
