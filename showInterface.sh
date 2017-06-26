#for i in "::BasicInfo("   "::Traj(" "::Reader(" "::XmlReader(" "::BinReader(" "::XmlModify(" "::AnalFunc(" "::SegRgDis(" "::Analysis(" "::BoxLen(" "::Traj2CMBin(" 
#for i in $(grep export_ tackle.cc | awk -F '[_)]' '{printf "\::""$2"\" "}END{"\n"}' )
for i in $(grep export_ tackle.cc | awk -F '[_)]' '{printf "::"$2" "}END{"\n"}'  )
do
#echo $i
  grep $i  *.cc  *.h particles/*.cc particles/*.h
  echo " ------------------------------- "
done
