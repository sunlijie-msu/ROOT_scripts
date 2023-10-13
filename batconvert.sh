 #!/bin/bash
 #rootdir=/media/Si22/rootfile
 num1=200
 num2=5400
 #num1<num2
 for((i=${num1};i<=${num2};i+=200));
 do
 echo Procesing file_lifetime ${i} ....
 /mnt/analysis/s1582/Simulation/source/sim_dsam /mnt/analysis/s1582/Simulation/InputCards/S31_Eg1248_tau${i}_EMG.txt
 done                                        

 #try this command first: chmod a+x batconvert.sh
