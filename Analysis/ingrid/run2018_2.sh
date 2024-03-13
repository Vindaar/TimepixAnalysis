#!/usr/bin/zsh
dbPath=/home/basti/CastData/ExternCode/TimepixAnalysis/InGridDatabase/src/resources
cp $dbPath/ingridDatabase2018.h5 $dbPath/ingridDatabase.h5
echo "Raw data manipulation"
#./raw_data_manipulation /mnt/1TB/CAST/2018_2/CalibrationRuns/ --runType=calib --out=/mnt/1TB/CAST/2018_2/CalibrationRuns.h5
#./raw_data_manipulation /mnt/1TB/CAST/2018_2/DataRuns/ --runType=back --out=/mnt/1TB/CAST/2018_2/DataRuns.h5
echo "Reconstruction"
#./reconstruction /mnt/1TB/CAST/2018_2/CalibrationRuns.h5
#./reconstruction /mnt/1TB/CAST/2018_2/DataRuns.h5
echo "Calculating FADC properties"
#./reconstruction /mnt/1TB/CAST/2018_2/CalibrationRuns.h5 --only_fadc
#./reconstruction /mnt/1TB/CAST/2018_2/DataRuns.h5 --only_fadc
#echo "Applying TOT calibration"
#./reconstruction /mnt/1TB/CAST/2018_2/DataRuns.h5 --only_charge
## already done by default in reconstruction
## ./reconstruction /mnt/1TB/CAST/2018_2/CalibrationRuns.h5 --only_charge
#echo "Calculating gas gain"
#./reconstruction /mnt/1TB/CAST/2018_2/CalibrationRuns.h5 --only_gas_gain
#./reconstruction /mnt/1TB/CAST/2018_2/DataRuns.h5 --only_gas_gain
#echo "Fitting gas gain on calibration run"
./reconstruction /mnt/1TB/CAST/2018_2/CalibrationRuns.h5 --only_gain_fit
#echo "Calculating energy from charge"
./reconstruction /mnt/1TB/CAST/2018_2/CalibrationRuns.h5 --only_energy_from_e
./reconstruction /mnt/1TB/CAST/2018_2/DataRuns.h5 --only_energy_from_e


echo "Performing Likelihood"
# perform likelihood cuts
./likelihood /mnt/Daten/CAST/2014_15/DataRuns.h5 --reference /mnt/Daten/CAST/CDL-reference/XrayReferenceDataSet.h5 --h5out 2014_15_jan4_2018.h5
./likelihood /mnt/1TB/CAST/2018_2/DataRuns.h5 --reference /mnt/Daten/CAST/CDL-reference/XrayReferenceDataSet.h5 --h5out 2018_2_jan4_2018.h5


# creating plots
cd /home/basti/CastData/ExternCode/TimepixAnalysis/PlottingPython/Plotting/
python3 ./PyS_createBackgroundRate.py ../../Analysis/ingrid/2014_15_jan4_2018.h5 --file2 ../../Analysis/ingrid/2018_2_jan4_2018.h5 --chip 3 --fancy
