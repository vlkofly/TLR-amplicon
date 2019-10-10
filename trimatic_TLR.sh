#! /bin/bash

for f in `cat forwardsamplelist.txt`
do
base=${f/'000000000-B4M7B_TLRmockers_17s003050-1-1_Vlcek_lane1'/''}
base=${base/'_sequence.txt.gz'}
logname='../outtrim/'$base'log.txt'
basename='../outtrim/'$base'fastq.gz'



java -jar ~/programs/Trimmomatic-0.36/trimmomatic-0.36.jar PE \
-threads 7 -trimlog $logname -basein $f -baseout $basename \
ILLUMINACLIP:../adapters.txt:5:20:8 TRAILING:25 SLIDINGWINDOW:3:25 MINLEN:120
echo $base
done >> ../outtrim/trimlog.txt 2>&1
