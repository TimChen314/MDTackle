
dirname="$(basename $(pwd))"
release_name="${dirname}_release"
if [  -d $release_name ];then
    rm -r $release_name
fi
if [  -d ${release_name}.tar.gz ];then
    rm -r ${release_name}.tar.gz
fi
mkdir  -p ../$release_name/particles
mkdir  ../$release_name/compile
mkdir  ../$release_name/test

cp configure *.h *.cc ../$release_name
cp compile/compilefile ../$release_name/compile
cp particles/*.h particles/*.cc ../$release_name/particles
cp -r ~/PyTackleTest/standard_* ../$release_name/test
cp ~/PyTackleTest/coor.xml ~/PyTackleTest/coor_type3.xml ~/PyTackleTest/Test.py ~/PyTackleTest/input*so ~/PyTackleTest/particle.00100[12]0000.bin ../$release_name/test

#tar -zcf ../${release_name}.tar.gz ../${release_name} 
#rm -r ../$release_name

