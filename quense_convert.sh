#!/bin/bash
rootdir=/media/Si23/rootfile
echo please input the start num:
read start_num
echo please input the stop num:
read stop_num
echo start_num is ${start_num}, stop_num is ${stop_num}
for((i=${start_num};i<=${stop_num};i++));
do
echo Procesing file_num ${i} ....
/home/daq/Ribll_401/r2root/r2root/convert ${i} >>${rootdir}/data0${i}.log
done
echo congratulation all the file from num ${start_num} to num ${stop_num} have been converted



 #!/bin/bash
 rootdir=/media/Si22/rootfile
 num1=516
 num2=518
 #num1<num2
 for((i=${num1};i<=${num2};i++));
 do /home/daq/Ribll_401/r2root/r2root/convert ${i} >>${rootdir}/data0${i}.log
 done                                        
