#!/bin/bash
#keys=`find work/* -type f -name "*sorted*" -print`
#keys=`find work/* -type f -name "*merged*" -print`
#keys=`find work/* -type f -name "*deduped*" -print`
keys=`find work/* -type f -name "*indels*" -print`
#keys=`find work/* -type f -name "*bqsr*" -print`
#keys=`find work/* -type f -name "*clean*" -print`
#keys=`find work/* -type f -name "*.vcf*" -print`


#dir=/wittelab/data2/carioc/wgs-gatk/tumor-normal/aligned;
#dir=/wittelab/data2/carioc/wgs-gatk/tumor-normal/merged;
dir=/wittelab/data2/carioc/wgs-gatk/tumor-normal/cleaned;
#dir=/wittelab/data2/carioc/wgs-gatk/tumor-normal/called;

mkdir -p $dir;
for element in ${keys[@]};
do
    echo ${element};
    mv ${element} $dir;
done
