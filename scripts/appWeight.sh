BASE=/afs/cern.ch/work/p/pluca/Lmumu
CASTORBASE=root://castorlhcb.cern.ch//castor/cern.ch/user/p/pluca/Lmumu
OUTDIR=/afs/cern.ch/work/p/pluca/Lmumu/weighted
MICHAL2012=/afs/cern.ch/work/k/kreps/public/LbLMuMuAna/early2012HLT

./trainLb_NeuroBayes.out > log_train


#B -> KS and Kst MC

#./appNB_Lb.out data $OUTDIR/Bd2JpsiKS_MC12_NBweighted.root $BASE/Bd2JpsiKS_MC12_Strip20_Pythia6_MagDown.root $BASE/Bd2JpsiKS_MC12_Strip20_Pythia6_MagUp.root

#./appNB_Lb.out data $OUTDIR/Bd2KSmumu_MC12_NBweighted.root $BASE/Bd2KSmumu_MC12_Strip20_Pythia6_MagDown.root $BASE/Bd2KSmumu_MC12_Strip20_Pythia6_MagUp.root

#./appNB_Lb.out data $OUTDIR/Bu2Kstmumu_MC12_NBweighted.root $BASE/Bu2Kstmumu_MC12_Pythia8_MagDown.root $BASE/Bu2Kstmumu_MC12_Pythia8_MagUp.root


#$OUTDIR/add.sh

#./optimizeMVAcut.out > log_optimize

#./plotTrain.out



#./appNB_Lb.out Data orig $OUTDIR/Lb2Lmumu_CL_NBweighted.root $CASTORBASE/Lb2Lmumu_CL12_MagDown.root $CASTORBASE/Lb2Lmumu_CL12_MagUp.root $CASTORBASE/Lb2Lmumu_CL11_MagDown.root $CASTORBASE/Lb2Lmumu_CL11_MagUp.root

#./appNB_Lb.out MC orig $OUTDIR/Lb2JpsiL_MC_Pythia8_NBweighted.root $BASE/Lb2JpsiL_MC12_Pythia8_MagDown.root $BASE/Lb2JpsiL_MC12_Pythia8_MagUp.root $MICHAL2012/jpsiLMagUpP8.root $MICHAL2012/jpsiLMagDownP8.root $BASE/Lb2JpsiL_MC11_Pythia8_MagDown.root $BASE/Lb2JpsiL_MC11_Pythia8_MagUp.root



#./appNB_Lb.out MC orig MCtree $OUTDIR/Lb2Lmumu_MC_Pythia8_NBweighted_MCtree_12.root $BASE/Lb2Lmumu_MC12_Pythia8_MagDown.root $BASE/Lb2Lmumu_MC12_Pythia8_MagUp.root $MICHAL2012/LmumuMagUpP8.root $MICHAL2012/LmumuMagDownP8.root
#./appNB_Lb.out MC orig MCtree $OUTDIR/Lb2Lmumu_MC_Pythia8_NBweighted_MCtree_11.root $BASE/Lb2Lmumu_MC11_Pythia8_MagDown.root $BASE/Lb2Lmumu_MC11_Pythia8_MagUp.root

#./appNB_Lb.out MC orig tree $OUTDIR/Lb2Lmumu_MC_Pythia8_NBweighted_tree.root $BASE/Lb2Lmumu_MC12_Pythia8_MagDown.root $BASE/Lb2Lmumu_MC12_Pythia8_MagUp.root $MICHAL2012/LmumuMagUpP8.root $MICHAL2012/LmumuMagDownP8.root  $BASE/Lb2Lmumu_MC11_Pythia8_MagDown.root $BASE/Lb2Lmumu_MC11_Pythia8_MagUp.root

#./appNB_Lb.out MC orig MCDecaytree $OUTDIR/Lb2Lmumu_MC_Pythia8_NBweighted_MC12_MD.root $BASE/Lb2Lmumu_MC12_Pythia8_MagDown.root
#./appNB_Lb.out MC orig MCDecaytree $OUTDIR/Lb2Lmumu_MC_Pythia8_NBweighted_MC12_MU.root $BASE/Lb2Lmumu_MC12_Pythia8_MagUp.root
#./appNB_Lb.out MC orig MCDecaytree $OUTDIR/Lb2Lmumu_MC_Pythia8_NBweighted_MC11.root $BASE/Lb2Lmumu_MC11_Pythia8_MagDown.root $BASE/Lb2Lmumu_MC11_Pythia8_MagUp.root
#./appNB_Lb.out MC orig MCDecaytree $OUTDIR/Lb2Lmumu_MC_Pythia8_NBweighted_early12.root $MICHAL2012/LmumuMagUpP8.root $MICHAL2012/LmumuMagDownP8.root







#Last to use!!!


#./appNBW_Lb.out sigTestSample $OUTDIR/sigTestSample.root $OUTDIR/samplesMVA_12.root
#./appNBW_Lb.out sigTrainSample $OUTDIR/sigTrainSample.root $OUTDIR/samplesMVA_12.root
#./appNBW_Lb.out bkgTestSample $OUTDIR/bkgTestSample.root $OUTDIR/samplesMVA_12.root
#./appNBW_Lb.out bkgTrainSample $OUTDIR/bkgTrainSample.root $OUTDIR/samplesMVA_12.root

#./addAngles.out $OUTDIR/sigTestSample_angles.root $OUTDIR/sigTestSample.root sigTestSample
#./addAngles.out $OUTDIR/bkgTestSample_angles.root $OUTDIR/bkgTestSample.root bkgTestSample

#rm $OUTDIR/trainingSamples.root
#hadd $OUTDIR/trainingSamples.root $OUTDIR/bkgTrainSample.root $OUTDIR/bkgTestSample.root $OUTDIR/sigTestSample.root $OUTDIR/sigTrainSample.root
#rm $OUTDIR/sigT* $OUTDIR/bkgT*



#./appNB_Lb.out MC ww $OUTDIR/Lb2JpsiL_MC_Pythia8_NBweighted.root $WORK/weighted/Lmumu/old/Lb2JpsiL_MC_Pythia8_NBweighted.root
#./appNB_Lb.out MC ww $OUTDIR/Lb2Lmumu_MC_Pythia8_NBweighted.root $WORK/weighted/Lmumu/old/Lb2Lmumu_MC_Pythia8_NBweighted.root

#./appNB_Lb.out data orig $OUTDIR/Bu2Kstmumu_MC12_NBweighted.root $BASE/Bu2Kstmumu_MC12_Pythia8_MagDown.root $BASE/Bu2Kstmumu_MC12_Pythia8_MagUp.root
#./appNB_Lb.out data orig $OUTDIR/Bd2JpsiKS_MC12_NBweighted.root $BASE/Bd2JpsiKS_MC12_Pythia8_MagDown.root $BASE/Bd2JpsiKS_MC12_Pythia8_MagUp.root
#./appNB_Lb.out data orig $OUTDIR/Bd2KSmumu_MC12_NBweighted.root $BASE/Bd2KSmumu_MC12_Pythia8_MagDown.root $BASE/Bd2KSmumu_MC12_Pythia8_MagUp.root

#./appNB_Lb.out MCGeom ww all $OUTDIR/Lb2JpsiL_geomMC_Pythia8_NBweighted_new.root $WORK/Lmumu/weighted/Lb2JpsiL_geomMC_Pythia8_NBweighted.root
./appNB_Lb.out MCGeom ww all $OUTDIR/Lb2Lmumu_geomMC_Pythia8_NBweighted_new.root $WORK/Lmumu/weighted/Lb2Lmumu_geomMC_Pythia8_NBweighted.root


