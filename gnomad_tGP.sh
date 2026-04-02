#!/bin/bash

# Query
cd /mnt/nas_share/data_sb/tools/mutect2_ref/tGP_data_preparation
/mnt/nas_share/data_sb/tools/code/query_gnomAD.sh --gl ../tPanelGene.txt --data-type joint

# Download 
for i in *json ; do /mnt/nas_share/data_sb/tools/code/joint_flatten_gnomad.sh $i --data-type joint ; done

# Catenate all genes
awk -F"\t" 'FNR==1 && NR != 1 {next} {print} ' *txt > all_tGP.txt

# Prepare VCF for Mutect2 
awk -F'\t' 'BEGIN{
OFS="\t";
print "##fileformat=VCFv4.2"

for(i=1;i<=22;i++) print "##contig=<ID=chr"i">"
print "##contig=<ID=chrX>"
print "##contig=<ID=chrY>"

print "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele Count\">"
print "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Allele Number\">"
print "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">"
print "##INFO=<ID=AC_Global,Number=A,Type=Integer,Description=\"Global Allele Count\">"
print "##INFO=<ID=AN_Global,Number=1,Type=Integer,Description=\"Global Allele Number\">"
print "##INFO=<ID=AF_Global,Number=A,Type=Float,Description=\"Global Allele Frequency\">"
print "##INFO=<ID=AC_SAS,Number=A,Type=Integer,Description=\"South Asian Allele Count\">"
print "##INFO=<ID=AN_SAS,Number=1,Type=Integer,Description=\"South Asian Allele Number\">"
print "##INFO=<ID=AF_SAS,Number=A,Type=Float,Description=\"South Asian Allele Frequency\">"

print "#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO"
}
FNR==1{next}
{
split($1,a,"-")
if(length(a)==4){

ID=($4=="" || $4=="NA")?".":$4

AC=$11; AN=$14
AC_SAS=$41; AN_SAS=$42

AF_Global=(AN?AC/AN:0)
AF_SAS=(AN_SAS?AC_SAS/AN_SAS:0)

if(AF_SAS>AF_Global){ACF=AC_SAS; ANF=AN_SAS; AFF=AF_SAS}
else{ACF=AC; ANF=AN; AFF=AF_Global}

print "chr"a[1],a[2],ID,a[3],a[4],".",".", \
"AC="ACF";AN="ANF";AF="sprintf("%.8f",AFF) \
";AC_Global="AC";AN_Global="AN";AF_Global="sprintf("%.8f",AF_Global) \
";AC_SAS="AC_SAS";AN_SAS="AN_SAS";AF_SAS="sprintf("%.8f",AF_SAS)
}
}' all_tGP.txt | bcftools sort -Oz -o ../af_gnomad41_tGP.vcf.gz

