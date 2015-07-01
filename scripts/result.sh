
#./createTuple.sh

###################  BR Analysis

#~/python/submit.py -q 1nd -d yield/comb Lb_yield.out
#~/python/submit.py -q 1nd -d yield/comb -n Lb_yield_norm_plus 'Lb_yield.out -Snorm_plus'
#~/python/submit.py -q 1nd -d yield/comb -n Lb_yield_norm_minus 'Lb_yield.out -Snorm_minus'
#~/python/submit.py -q 1nd -d yield/comb -n Lb_yield_eff_plus 'Lb_yield.out -Seff_plus -W'
#~/python/submit.py -q 1nd -d yield/comb -n Lb_yield_eff_minus 'Lb_yield.out -Seff_minus -W'
#~/python/submit.py -q 1nd -d yield/LL 'Lb_yield.out -tLL'
#~/python/submit.py -q 1nd -d yield/LL -n Lb_yield_norm_plus 'Lb_yield.out -Snorm_plus -tLL'
#~/python/submit.py -q 1nd -d yield/LL -n Lb_yield_norm_minus 'Lb_yield.out -Snorm_minus -tLL'
#~/python/submit.py -q 1nd -d yield/LL -n Lb_yield_eff_plus 'Lb_yield.out -Seff_plus -W -tLL'
#~/python/submit.py -q 1nd -d yield/LL -n Lb_yield_eff_minus 'Lb_yield.out -Seff_minus -W -tLL'
#~/python/submit.py -q 1nd -d yield/DD 'Lb_yield.out -tDD'
#~/python/submit.py -q 1nd -d yield/DD -n Lb_yield_norm_plus 'Lb_yield.out -Snorm_plus -tDD'
#~/python/submit.py -q 1nd -d yield/DD -n Lb_yield_norm_minus 'Lb_yield.out -Snorm_minus -tDD'
#~/python/submit.py -q 1nd -d yield/DD -n Lb_yield_eff_plus 'Lb_yield.out -Seff_plus -W -tDD'
#~/python/submit.py -q 1nd -d yield/DD -n Lb_yield_eff_minus 'Lb_yield.out -Seff_minus -W -tDD'

#~/python/submit.py -q 1nd -d yield/LL -n Lb_yield_eff_plus_noWC 'Lb_yield.out -Seff_plus -tLL'
#~/python/submit.py -q 1nd -d yield/LL -n Lb_yield_eff_minus_noWC 'Lb_yield.out -Seff_minus -tLL'
#~/python/submit.py -q 1nd -d yield/DD -n Lb_yield_eff_plus_noWC 'Lb_yield.out -Seff_plus -tDD'
#~/python/submit.py -q 1nd -d yield/DD -n Lb_yield_eff_minus_noWC 'Lb_yield.out -Seff_minus -tDD'

~/python/submit.py -q 1nd -d yield/comb -n Lb_yield_eff_plus_noWC 'Lb_yield.out -Seff_plus'
~/python/submit.py -q 1nd -d yield/comb -n Lb_yield_eff_minus_noWC 'Lb_yield.out -Seff_minus'

echo 'DDsys 2bin'
#./Lb_sys.out -sys -percent -tDD -r -c[1.1,6.0,15.0,20.0] '-Cpplus_ID>0' -oLbreleffAndSysvsq2_DD_2bins.root > log_Rel_DD_2bins
echo 'LLsys 2bin'
#./Lb_sys.out -sys -percent -tLL -r -c[1.1,6.0,15.0,20.0] '-Cpplus_ID>0' -oLbreleffAndSysvsq2_LL_2bins.root > log_Rel_LL_2bins
echo 'DDsys 7bin'
#./Lb_sys.out -sys -percent -tDD -r '-Cpplus_ID>0' -oLbreleffAndSysvsq2_DD_7bins.root > log_Rel_DD_7bins
echo 'LLsys 7bin'
#./Lb_sys.out -sys -percent -tLL -r '-Cpplus_ID>0' -oLbreleffAndSysvsq2_LL_7bins.root > log_Rel_LL_7bins
#mv log* $WORK/results/
#mv Lb*.root $WORK/results/
#./extractLbresult.out > log_BR_results
#./extractLbresult.out -peff >> log_BR_results

#######################################

echo 'jpsi'
#./Lb_sys.out -jpsi -tLL > log_Abs_LL_jpsi
#./Lb_sys.out -jpsi -tDD > log_Abs_DD_jpsi
echo 'DDsys 2bin'
#./Lb_sys.out -sys -percent -tDD -c[1.1,6.0,15.0,20.0] -oLbeffAndSysvsq2_DD_2bins.root > log_Abs_DD_2bins
echo 'LLsys 2bin'
#./Lb_sys.out -sys -percent -tLL -c[1.1,6.0,15.0,20.0] -oLbeffAndSysvsq2_LL_2bins.root > log_Abs_LL_2bins
echo 'DDsys 7bin'
#./Lb_sys.out -sys -percent -tDD -oLbeffAndSysvsq2_DD_7bins.root > log_Abs_DD_7bins
echo 'LLsys 7bin'
#./Lb_sys.out -sys -percent -tLL -oLbeffAndSysvsq2_LL_7bins.root > log_Abs_LL_7bins
#mv Lb*.root $WORK/results/
#mv log* $WORK/results/


########## Angular Analysis

#./Lb_2Deff.out -XcosThetaB -x[10,-1,1] -YcosTheta -y[10,-1,1] -Dcolztext
#./Lb_2Deff.out -XcosThetaL -x[10,-1,1] -YcosTheta -y[10,-1,1] -Dcolztext
#./Lb_2Deff.out -XcosThetaL -x[10,-1,1] -YcosThetaB -y[10,-1,1] -Dcolztext
#./Lb_2Deff.out -XcosThetaB -x[10,-1,1] -YcosTheta -y[10,-1,1] -Dcolztext -tDD
#./Lb_2Deff.out -XcosThetaL -x[10,-1,1] -YcosTheta -y[10,-1,1] -Dcolztext -tDD
#./Lb_2Deff.out -XcosThetaL -x[10,-1,1] -YcosThetaB -y[10,-1,1] -Dcolztext -tDD
#./Lb_2Deff.out -XcosThetaB -x[10,-1,1] -YcosTheta -y[10,-1,1] -Dcolztext -tLL
#./Lb_2Deff.out -XcosThetaL -x[10,-1,1] -YcosTheta -y[10,-1,1] -Dcolztext -tLL
#./Lb_2Deff.out -XcosThetaL -x[10,-1,1] -YcosThetaB -y[10,-1,1] -Dcolztext -tLL

#./Lb_2Deff.out -YcosThetaL -y[30,-1,1] -x2bin -tDD -Dcolztext -OLbeff2D_cosThetaL_vs_q2_DD_2bins.root
#./Lb_2Deff.out -YcosThetaL -y[30,-1,1] -x2bin -tLL -Dcolztext -OLbeff2D_cosThetaL_vs_q2_LL_2bins.root
#./Lb_2Deff.out -YcosThetaB -y[30,-1,1] -x2bin -tDD -Dcolztext -OLbeff2D_cosThetaB_vs_q2_DD_2bins.root
#./Lb_2Deff.out -YcosThetaB -y[30,-1,1] -x2bin -tLL -Dcolztext -OLbeff2D_cosThetaB_vs_q2_LL_2bins.root

#./Lb_2Deff.out -YcosThetaL -y[30,-1,1] -tDD -Dcolztext
#./Lb_2Deff.out -YcosThetaL -y[30,-1,1] -tLL -Dcolztext
#./Lb_2Deff.out -YcosThetaB -y[30,-1,1] -tDD -Dcolztext
#./Lb_2Deff.out -YcosThetaB -y[30,-1,1] -tLL -Dcolztext

#./Lb_sys.out -XcosThetaL -tDD -jpsi -b[30,-1,1] -Dcolztext
#./Lb_sys.out -XcosThetaL -tLL -jpsi -b[30,-1,1] -Dcolztext
#./Lb_sys.out -XcosThetaB -tDD -jpsi -b[30,-1,1] -Dcolztext
#./Lb_sys.out -XcosThetaB -tLL -jpsi -b[30,-1,1] -Dcolztext

#mv Lb*.root $WORK/results/
#mv log* $WORK/results/

#./AnafitAfb_nosW.out > log_Afb

#mv *eff*.pdf plots/effs
#mv *fit*.pdf plots/fits
#mv *.pdf plots/results

