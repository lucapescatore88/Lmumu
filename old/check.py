from ROOT import *


rare_yields = {'Lb' : { 'LL' :  3.13545, 'DD' : 11.8138 }, 'Lbbar': { 'LL' :  2.8047, 'DD' : 10.9496 }}
jpsi_yields = {'Lb' : { 'LL' :  2.2132, 'DD' : 7.1118 }, 'Lbbar': { 'LL' :  2.0126, 'DD' : 6.34853 }}

#jpsi_yields = {'Lb' : { 'LL' :  5919, 'DD' : 20860 }, 'Lbbar': { 'LL' :  5382, 'DD' : 18742 }}
#rare_yields = {'Lb' : { 'LL' :  2092, 'DD' : 8038 }, 'Lbbar': { 'LL' :  1837, 'DD' : 7460 }}

## for Lb
ratio_LL_Lb = float(rare_yields['Lb']['LL']) / jpsi_yields['Lb']['LL']
ratio_DD_Lb = float(rare_yields['Lb']['DD']) / jpsi_yields['Lb']['DD']

## for Lbbar
ratio_LL_Lbbar = float(rare_yields['Lbbar']['LL']) / jpsi_yields['Lbbar']['LL']
ratio_DD_Lbbar = float(rare_yields['Lbbar']['DD']) / jpsi_yields['Lbbar']['DD']

double_ratio_LL = ratio_LL_Lb / ratio_LL_Lbbar
double_ratio_DD = ratio_DD_Lb / ratio_DD_Lbbar

print 'Ratios Lb'
print 'LL', ratio_LL_Lb, 'DD', ratio_DD_Lb

print 'Ratios Lb'
print 'LL', ratio_LL_Lbbar, 'DD', ratio_DD_Lbbar

print 'Double ratios'
print 'LL', double_ratio_LL, 'DD', double_ratio_DD




##Data

Lb_LL = {'val' : 0.0118, 'err' : 0.0027 }
Lb_DD = {'val' : 0.0231, 'err' : 0.0023 }

Lbbar_LL = {'val' : 0.0127, 'err' : 0.0029 }
Lbbar_DD = {'val' : 0.0157, 'err' : 0.0021 }

ratio_LL = Lb_LL['val']/Lbbar_LL['val'];
ratio_DD = Lb_DD['val']/Lbbar_DD['val'];

ratio_LL_err = ratio_LL*(TMath.Sqrt(TMath.Power(Lb_LL['err']/Lb_LL['val'],2) + TMath.Power(Lbbar_LL['err']/Lbbar_LL['val'],2)));
ratio_DD_err = ratio_DD*(TMath.Sqrt(TMath.Power(Lb_DD['err']/Lb_DD['val'],2) + TMath.Power(Lbbar_DD['err']/Lbbar_DD['val'],2)));

print "Data double ratios"
print ratio_LL, "+-", ratio_LL_err
print ratio_DD, "+-", ratio_DD_err



