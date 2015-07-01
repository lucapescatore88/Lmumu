echo Abs geometry

grep eff $WORK/results/log_Abs_DD_7bins | awk '{ print $2, "\t& ", $4, $5, $6, " \\\\" }'
echo \\hline
grep eff $WORK/results/log_Abs_DD_2bins | awk '{ print $2, "\t& ", $4, $5, $6, " \\\\" }'


echo ''
echo Abs detection and reco

grep eff $WORK/results/log_Abs_DD_7bins | awk '{ print $2, "\t& ", $8, $9, $10 }' > logDet
grep eff $WORK/results/log_Abs_DD_2bins | awk '{ print $2, "\t& ", $8, $9, $10 }' >> logDet

grep eff $WORK/results/log_Abs_DD_7bins | awk '{ print "\t& ", $12, $13, $14 }' > logDD
grep eff $WORK/results/log_Abs_DD_2bins | awk '{ print "\t& ", $12, $13, $14 }' >> logDD

grep eff $WORK/results/log_Abs_LL_7bins | awk '{ print "\t& ", $12, $13, $14, " \\\\" }' > logLL
grep eff $WORK/results/log_Abs_LL_2bins | awk '{ print "\t& ", $12, $13, $14, " \\\\" }' >> logLL

python joinFiles.py logDet logDD >> logTmp
python joinFiles.py logTmp logLL




echo ''
echo Abs MVA

grep eff $WORK/results/log_Abs_DD_7bins | awk '{ print $2, "\t& ", $16, $17, $18 }' > logDD
grep eff $WORK/results/log_Abs_DD_2bins | awk '{ print $2, "\t& ", $16, $17, $18 }' >> logDD

grep eff $WORK/results/log_Abs_LL_7bins | awk '{ print "\t& ", $16, $17, $18, " \\\\" }' > logLL
grep eff $WORK/results/log_Abs_LL_2bins | awk '{ print "\t& ", $16, $17, $18, " \\\\" }' >> logLL

python joinFiles.py logDD logLL



echo ''
echo Abs trig

grep eff $WORK/results/log_Abs_DD_7bins | awk '{ print $2, "\t& ", $20, $21, $22 }' > logDD
grep eff $WORK/results/log_Abs_DD_2bins | awk '{ print $2, "\t& ", $20, $21, $22 }' >> logDD

grep eff $WORK/results/log_Abs_LL_7bins | awk '{ print "\t& ", $20, $21, $22, " \\\\" }' > logLL
grep eff $WORK/results/log_Abs_LL_2bins | awk '{ print "\t& ", $20, $21, $22, " \\\\" }' >> logLL

python joinFiles.py logDD logLL



echo ''
echo Abs total

grep eff $WORK/results/log_Abs_DD_7bins | awk '{ print $2, "\t& ", $28, $29, $30 }' > logDD
grep eff $WORK/results/log_Abs_DD_2bins | awk '{ print $2, "\t& ", $28, $29, $30 }' >> logDD

grep eff $WORK/results/log_Abs_LL_7bins | awk '{ print "\t& ", $28, $29, $30, " \\\\" }' > logLL
grep eff $WORK/results/log_Abs_LL_2bins | awk '{ print "\t& ", $28, $29, $30, " \\\\" }' >> logLL

python joinFiles.py logDD logLL


########################################################################################################


echo ''
echo Jpsi abs 

grep eff $WORK/results/log_Abs_DD_jpsi | awk '{ print $4, $5, $6, "\n", $8, $9, $10, "\n",$12, $13, $14, "\n",$16, $17, $18, "\n", $20, $21, $22, "\n", $24, $25, $26, "\n", $28, $29, $30 }' > logDD
grep eff $WORK/results/log_Abs_LL_jpsi | awk '{ print "\t& ", $4, $5, $6, " \\\\", "\n&", $8, $9, $10, " \\\\", "\n&",$12, $13, $14, " \\\\", "\n&",$16, $17, $18, " \\\\", "\n&", $20, $21, $22, " \\\\",  "\n&", $24, $25, $26, " \\\\",  "\n&", $28, $29, $30, " \\\\" }' > logLL

python joinFiles.py logDD logLL






###################################################################################################################

echo ''
echo Rel geometry and detection

grep eff $WORK/results/log_Rel_DD_7bins | awk '{ print $2, "\t& ", $4, $5, $6, "\t& ", $8, $9, $10, " \\\\" }'
echo \\hline
grep eff $WORK/results/log_Rel_DD_2bins | awk '{ print $2, "\t& ", $4, $5, $6, "\t& ", $8, $9, $10, " \\\\" }'

echo ''
echo Rel selection DD

grep eff $WORK/results/log_Rel_DD_7bins | awk '{ print $2, "\t& ", $12, $13, $14, "\t& ", $16, $17, $18, "\t& ", $20, $21, $22, "\t& ", $24, $25, $26, " \\\\" }'
echo \\hline
grep eff $WORK/results/log_Rel_DD_2bins | awk '{ print $2, "\t& ", $12, $13, $14, "\t& ", $16, $17, $18, "\t& ", $20, $21, $22, "\t& ", $24, $25, $26, " \\\\" }'

echo ''
echo Rel selection LL

grep eff $WORK/results/log_Rel_LL_7bins | awk '{ print $2, "\t& ", $12, $13, $14, "\t& ", $16, $17, $18, "\t& ", $20, $21, $22, "\t& ", $24, $25, $26, " \\\\" }'
echo \\hline
grep eff $WORK/results/log_Rel_LL_2bins | awk '{ print $2, "\t& ", $12, $13, $14, "\t& ", $16, $17, $18, "\t& ", $20, $21, $22, "\t& ", $24, $25, $26, " \\\\" }'

echo ''
echo Rel total

grep eff $WORK/results/log_Rel_DD_7bins | awk '{ print $2, "\t& ", $28, $29, $30 }' > logDD
grep eff $WORK/results/log_Rel_DD_2bins | awk '{ print $2, "\t& ", $28, $29, $30 }' >> logDD

grep eff $WORK/results/log_Rel_LL_7bins | awk '{ print "\t& ", $28, $29, $30, " \\\\" }' > logLL
grep eff $WORK/results/log_Rel_LL_2bins | awk '{ print "\t& ", $28, $29, $30, " \\\\" }' >> logLL

python joinFiles.py logDD logLL

rm logTmp logDD logLL logDet




