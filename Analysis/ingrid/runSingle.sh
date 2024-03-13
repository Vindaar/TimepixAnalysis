#!/usr/bin/zsh
dbPath=/home/basti/CastData/ExternCode/TimepixAnalysis/InGridDatabase/src/resources
cp $dbPath/ingridDatabase2018.h5 $dbPath/ingridDatabase.h5
runName=`basename $1`
dirName=`dirname $1`
outName=$dirName/$runName.h5
echo "Raw data manipulation of" $outName
./raw_data_manipulation /mnt/1TB/CAST/CDL_2019/$runName --runType=calib --out=$outName
echo "Reconstruction"
./reconstruction $outName
echo "Applying TOT calibration"
./reconstruction $outName --only_charge
echo "Fitting Fe spectra"
./reconstruction $outName --only_fe_spec
echo "Calculating FADC properties"
./reconstruction $outName --only_fadc
echo "Calculating gas gain"
./reconstruction $outName --only_gas_gain
echo "Calculating energy from charge"
./reconstruction $outName --only_energy_from_e
