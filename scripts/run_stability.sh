#./Lb_yield.out -highQ2only -mLb -useNcand
#./Lb_yield.out -highQ2only -mLbbar -useNcand

#./Lb_yield.out -highQ2only -y2011 -useNcand
#./Lb_yield.out -highQ2only -y2012 -useNcand

#./Lb_yield.out -highQ2only '-p1' -useNcand
#./Lb_yield.out -highQ2only '-p-1' -useNcand

#./Lb_yield.out -highQ2only '-p1' -mLb -useNcand
#./Lb_yield.out -highQ2only '-p-1' -mLb -useNcand
#./Lb_yield.out -highQ2only '-p1' -mLbbar -useNcand
#./Lb_yield.out -highQ2only '-p-1' -mLbbar -useNcand

#./Lb_yield.out -highQ2only '-p1' -y2011 -useNcand
#./Lb_yield.out -highQ2only '-p-1' -y2011 -useNcand
#./Lb_yield.out -highQ2only '-p1' -y2012 -useNcand
#./Lb_yield.out -highQ2only '-p-1' -y2012 -useNcand

#./Lb_yield.out -highQ2only -y2011 -mLb -useNcand
#./Lb_yield.out -highQ2only -y2012 -mLb -useNcand
#./Lb_yield.out -highQ2only -y2011 -mLbbar -useNcand
#./Lb_yield.out -highQ2only -y2012 -mLbbar -useNcand

./Lb_yield.out -highQ2only -y2011 -c'' -useNcand


grep -n LL log* | sed -e 's/log//g' | sed -e 's/_//g' | sed -e 's/\.txt/ /g' | sed -e 's/:1://g' | sed -e 's/:2://g' | sed -e 's/rare\/jpsi/\t/g'

grep -n DD log* | sed -e 's/log//g' | sed -e 's/_//g' | sed -e 's/\.txt/ /g' | sed -e 's/:1://g' | sed -e 's/:2://g' | sed -e 's/rare\/jpsi/\t/g'


