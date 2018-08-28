#!/bin/bash -e
# wrapper script for launching 3 Illumina NSv4 workflows:
# 1) the "Resequencing" workflow on the normal sample
# 2) the "Resequencing" workflow on the tumour sample
# 3) the "TumourNormal" workflow using the results of 1) and 2), pending their completion

while getopts "b:c:n:NTs:t:f:d:" opt; do
    case $opt in
	b) basedir=$OPTARG ;;	# /scratch
	n) normalbampath=$OPTARG ;; # /path/to/LP1000130-DNA_D11
	c) cancerbampath=$OPTARG ;; # /path/to/LP1000130-DNA_D12
	N) skipN=1 ;; # skip germline resequencing, use existing mapping
	T) skipT=1 ;; # skip tumour resequencing, use existing mapping
	s) runNS=$OPTARG ;;   # runNS4.sh path
	t) TemplateDir=$OPTARG ;;# eg. /path/containing/SampleSheet-*.csv
	f) fit2tpl=$OPTARG ;;   # fit2V4template.sh path
	d) deliverydir="$OPTARG" ;;	# final delivery dir base
	\?) echo "Invalid option: -"$OPTARG"" >&2
	    exit 1;;
	: ) echo "Option -"$OPTARG" requires an arg." >&2
	    exit 1;;
    esac
done




/home/mzarowiecki/scratch/TRACERx/2_Recalling/run_cancer_analysis.20171026.sh -n /home/mzarowiecki/scratch/TRACERx/2_Recalling/Reruns/2153n -c /home/mzarowiecki/scratch/TRACERx/2_Recalling/Reruns/2153t -s /home/mzarowiecki/scratch/TRACERx/2_Recalling/runNS4.sh -t  /genomes/software/apps/runNSv4/templates -f /genomes/software/apps/runNSv4/fit2V4template.sh -d /home/pipeline/bertha/temp_magz 




#basedir=/scratch
normalbampath=/home/mzarowiecki/scratch/TRACERx/2_Recalling/Reruns/2153n
cancerbampath=/home/mzarowiecki/scratch/TRACERx/2_Recalling/Reruns/2153t
skipN=1
skipT=1 
runNS=/home/mzarowiecki/scratch/TRACERx/2_Recalling/runNS4.sh
TemplateDir=/genomes/software/apps/runNSv4/templates 
fit2tpl=/genomes/software/apps/runNSv4/fit2V4template.sh  
deliverydir=/home/pipeline/bertha/temp_magz


if [ $skipN ]; then sed -i s/"$normalid"_S1/"$normalid"/ SampleSheet.csv; fi
if [ $skipT ]; then sed -i s/"$cancerid"_S1/"$cancerid"/ SampleSheet.csv; fi






if [ -z $basedir ]; then basedir=/scratch; fi
if [ ! -d $basedir ]||[ -z $cancerbampath ]||[ -z $normalbampath ]; then
    echo "Cannot read Base dir $basedir or no sample name given"; exit 1
fi
date; echo "starting $0"; uname -a	# debug
echo "CHECK 1"
normalbamfull=$(find $normalbampath -maxdepth 3 -type f -name "*.bam")
normalbam=${normalbamfull##*/}
normalid=${normalbam%.bam}
cancerbamfull=$(find $cancerbampath -maxdepth 3 -type f -name "*.bam")
cancerbam=${cancerbamfull##*/}
cancerid=${cancerbam%.bam}

# create base directory for the 3 analysis
mkdir -p ${basedir}/Cancer${cancerid}_Normal${normalid}

echo "CHECK 2"
if [ -z $skipT ]; then
# 1. run Resequencing workflow on tumour
  echo "WENT HERE 1 $skipT"
  if [ $(wc -c < $cancerbamfull) -lt 210000000000 ]; then
    echo "WENT HERE 2"
    mkdir -p $basedir/Old${cancerid}
    rsync -a $cancerbamfull $basedir/Old${cancerid}
    rsync -a ${cancerbamfull}.bai $basedir/Old${cancerid}
    cancerbamfull="$basedir/Old${cancerid}/$cancerbam"
  fi
  mkdir -p ${basedir}/Cancer${cancerid}_Normal${normalid}/tumour
  cd ${basedir}/Cancer${cancerid}_Normal${normalid}/tumour
  command="$runNS -n $normalid -c $cancerid -b $basedir -m T"
  cat $TemplateDir/SampleSheet-T.csv|sed "s:#CANCER_ID#:$cancerid:g;s:#CANCER_BAM_PATH#:$cancerbamfull:g" > SampleSheet.csv
  echo $command; $command
  status=$?
  if [ $status -ne 0 ]; then
    echo "WENT HERE 3 $skipT"
    date; echo "rerun $command"; $command
    status=$?
    if [ $status -ne 0 ]; then echo "failed again: $status"
	rsync -av $PWD/analysis/ $deliverydir/../Cancer${cancerid}_Normal${normalid}/T
    fi
  fi
  rm -rf $basedir/Old$cancerid	# save space
else	# skip tumour resequencing, use existing files
  aPATH="${basedir}/Cancer${cancerid}_Normal${normalid}/tumour/analysis"
  mkdir -p $aPATH
  cancerbamdir=$(dirname $cancerbamfull)
  rsync -av ${cancerbamfull}* $cancerbamdir/../LP*.* $aPATH/
  grep ".bam$" $cancerbampath/md5sum.txt|cut -d ' ' -f 1 > $aPATH/${cancerbam}.md5
fi


if [ -z $skipN ]; then
echo "CHECK 3a"
# 2. run Resequencing workflow on normal
  mkdir -p $basedir/Old${normalid}
  rsync -a $normalbamfull $basedir/Old${normalid}
  rsync -a ${normalbamfull}.bai $basedir/Old${normalid}
  normalbamfull="$basedir/Old${normalid}/$normalbam"
  mkdir -p ${basedir}/Cancer${cancerid}_Normal${normalid}/normal
  cd ${basedir}/Cancer${cancerid}_Normal${normalid}/normal
  command="$runNS -n $normalid -c $cancerid -b $basedir -m N"
  cat $TemplateDir/SampleSheet-N.csv|sed "s:#NORMAL_ID#:$normalid:g;s:#NORMAL_BAM_PATH#:$normalbamfull:g" > SampleSheet.csv
  echo $command; $command
  status=$?
  if [ $status -ne 0 ]; then
    date; echo "rerun $command"; $command
    status=$?
    if [ $status -ne 0 ]; then echo "failed again: $status"
	rsync -av $PWD/analysis/ $deliverydir/../Cancer${cancerid}_Normal${normalid}/N
    fi
  fi
  rm -rf $basedir/Old${normalid}
else	# skip germline resequencing, use existing files
  echo "CHECK3b"
  aPATH="${basedir}/Cancer${cancerid}_Normal${normalid}/N/analysis"
  mkdir -p $aPATH
  rsync -av ${normalbamfull} ${normalbamfull}.bai $aPATH/
  grep ".bam$" $normalbampath/md5sum.txt|cut -d ' ' -f 1 > $aPATH/${normalbam}.md5
  rsync -av $normalbampath/V*/* $normalbampath/M*/* $normalbampath/LP*.* $aPATH/
  cd $aPATH/../../; ln -s N normal
fi

# 3. run TumourNormal workflow on normal + tumour
echo "CHECK4"
mkdir -p ${basedir}/Cancer${cancerid}_Normal${normalid}/tumour_normal
cd ${basedir}/Cancer${cancerid}_Normal${normalid}/tumour_normal
# ltime=$(date +%Y%m%d.%H%M)
command="$runNS -n $normalid -c $cancerid -b $basedir -m TN"
cat $TemplateDir/SampleSheet-TN.csv|sed "s:#BASE_DIR#:$basedir:g;s:#NORMAL_ID#:$normalid:g;s:#CANCER_ID#:$cancerid:g" > SampleSheet.csv
if [ $skipN ]; then sed -i s/"$normalid"_S1/"$normalid"/ SampleSheet.csv; fi
if [ $skipT ]; then sed -i s/"$cancerid"_S1/"$cancerid"/ SampleSheet.csv; fi
echo $command; $command
  status=$?
if [ $status -ne 0 ]; then
    date; echo "rerun $command"; $command
    status=$?
    if [ $status -ne 0 ]; then echo "failed again: $status"
	rsync -av $PWD/analysis/ $deliverydir/../Cancer${cancerid}_Normal${normalid}/TN
    fi
fi

# 4. fit to delivery template and md5
if [ ! -d ${basedir}/upload ]; then mkdir ${basedir}/upload; fi
sPATH=$(dirname $fit2tpl)
docker run -t -v /etc/localtime:/etc/localtime:ro -v $sPATH:$sPATH:ro -v $basedir:$basedir:rw illumina/isis:2.6.53.23 bash -c "$fit2tpl ${basedir}/Cancer${cancerid}_Normal${normalid}/ TN ${basedir}/upload"
status=$?
if [ $status -eq 0 ]; then
  if [ -z $skipN ]; then
    rsync -av ${basedir}/upload/${normalid} $deliverydir
  fi
    rsync -av ${basedir}/upload/Cancer${cancerid}_Normal${normalid} $deliverydir
    if [ $status -eq 0 ]; then
	echo -n 'finished rsync ';date
	docker run -t -v /etc/localtime:/etc/localtime:ro -v /scratch:/scratch:rw illumina/isis:2.6.53.23 bash -c "rm -rf ${basedir}/Cancer${cancerid}_Normal${normalid}/ $basedir/upload/Cancer${cancerid}_Normal${normalid} $basedir/upload/$normalid"
    fi
fi


\[\033[01;32m\]\u: \[\033[01;34m\]\W \[\033[01;34m\] \$ \[\033[0m\]



