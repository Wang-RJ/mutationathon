# Grab all "true" heterozygotes

sample=$(head -n 1 ../triotable.txt | tail -n 1)
bcftools view -Ou -M2 -m2 -s $sample ../vcf/output.filtered.snps.removed.CalculateGenotypePosteriors.vcf.gz | bcftools annotate -Oz -x INFO/CSQ,INFO/InbreedingCoeff,INFO/BaseQRankSum,INFO/ExcessHet,INFO/ClippingRankSum,INFO/FS,INFO/MLEAC,INFO/MLEAF,INFO/MQRankSum,INFO/ReadPosRankSum,INFO/SOR \
-i'(GT[0]="RR" && GT[1]="AA") || (GT[0]="AA" && GT[1]="RR")' -o hetcandidates.vcf.gz &

# Test whether inclusive filters on bcftools works correctly
# Depth first
bcftools view -i'FORMAT/DP[0] > 20 && FORMAT/DP[1] > 20 && FORMAT/DP[0] < 60 && FORMAT/DP[1] < 80 && FORMAT/GQ[0] > 50 && FORMAT/GQ[1] > 50' hetcandidates.vcf.gz | bcftools query -f '[%DP\t]\n' | \
awk '{if(min1==""){min1=max1=$1}; if($1>max1) {max1=$1}; if($1<min1) {min1=$1}; if(min2==""){min2=max2=$2}; if($2>max2) {max2=$2}; if($2<min2) {min2=$2}; if(min3==""){min3=max3=$3}; if($3>max3) {max3=$3}; if($3<min3) {min3=$3}} END {print max1, min1, max2, min2, max3, min3}'
# Now GQ
bcftools view -i'FORMAT/DP[0] > 20 && FORMAT/DP[1] > 20 && FORMAT/DP[0] < 60 && FORMAT/DP[1] < 80 && FORMAT/GQ[0] > 50 && FORMAT/GQ[1] > 50' hetcandidates.vcf.gz | bcftools query -f '[%GQ\t]\n' | \
awk '{if(min1==""){min1=max1=$1}; if($1>max1) {max1=$1}; if($1<min1) {min1=$1}; if(min2==""){min2=max2=$2}; if($2>max2) {max2=$2}; if($2<min2) {min2=$2}; if(min3==""){min3=max3=$3}; if($3>max3) {max3=$3}; if($3<min3) {min3=$3}} END {print max1, min1, max2, min2, max3, min3}'

# Create denominator by filtering to high quality sites and sampling roughly 250k
bcftools view -Oz -i'FORMAT/DP[0] > 20 && FORMAT/DP[1] > 20 && FORMAT/DP[0] < 60 && FORMAT/DP[1] < 80 && FORMAT/GQ[0] > 50 && FORMAT/GQ[1] > 50' hetcandidates.vcf.gz -o hetcandidates.filtered.vcf.gz

lines=$(zcat hetcandidates.filtered.vcf.gz | grep -v '#' | wc -l)
zcat hetcandidates.filtered.vcf.gz | awk -v lines=$lines 'BEGIN {srand()} (rand() * lines < 250000 || /^#/ )' | bcftools view -Oz -o hetdenom.vcf.gz

# Create numerator by filtering children from denominator
bcftools view -Oz -i'GT[2]="het"' hetdenom.vcf.gz -o hetnum_childhet.vcf.gz

# Filter for child GQ and DP (use same filters as for selecting candidates)
bcftools view -Oz -i'FORMAT/DP[2] > 20 && FORMAT/DP[2] < 80 && FORMAT/GQ[2] > 20' hetnum_childhet.vcf.gz -o hetnum_childhet_DPGQ20.vcf.gz

# In numerator folder
childbam=$(cut -d',' -f3 ../triotable.txt | head -n 1 | tail -n 1).dups.bam
bcftools mpileup -Oz -A -f /N/dcwan/projects/hahnlab-phi/macaque/ref/ref-relabeled.fna -Q 10 \
-R <(bcftools query hetnum_childhet_DPGQ20.vcf.gz -f'%CHROM\t%POS\n') \
-a ADF,ADR,AD,DP ../bam_files/$childbam > hetnum_mpileupGQ20.vcf.gz &

bcftools query hetnum_mpileupGQ20.vcf.gz -i'TYPE="snp" && FORMAT/ADF[:1] > 0 && FORMAT/ADR[:1] > 0' -f'%CHROM\t%POS\n' > pileup_childhetGQ20.txt
bcftools index hetnum_childhet_DPGQ20.vcf.gz
bcftools index hetnum_mpileupGQ20.vcf.gz


for GQ in {30,40,50}
do
	bcftools query hetnum_mpileupGQ20.vcf.gz -R <(bcftools query -i'FORMAT/GQ[2] > '$GQ hetnum_childhet_DPGQ20.vcf.gz -f '%CHROM\t%POS\n') -i'TYPE="snp" && FORMAT/ADF[:1] > 0 && FORMAT/ADR[:1] > 0' -f'%CHROM\t%POS\n' > pileup_childhetGQ$GQ.txt &
done

for GQ in {20,30,40,50}
do
  header="trio\tdenominator\ttransmit_filter\tdpgq_filter\tbam_filter\tab_filtered"
  echo -e $header > het_callability_tableGQ$GQ.txt
  i=1
  echo -e "$i\t\c" >> het_callability_tableGQ$GQ.txt
  echo -e "$(zcat hetdenom.vcf.gz | grep -v '#' | wc -l)\t\c" >> het_callability_tableGQ$GQ.txt
  echo -e "$(zcat hetnum_childhet.vcf.gz | grep -v '#' | wc -l)\t\c" >> het_callability_tableGQ$GQ.txt
  echo -e "$(bcftools view -i'FORMAT/GQ[2] > '$GQ hetnum_childhet_DPGQ20.vcf.gz | grep -v '#' | wc -l)\t\c" >> het_callability_tableGQ$GQ.txt
  echo -e "$(cat pileup_childhetGQ$GQ.txt | wc -l)\t\c" >> het_callability_tableGQ$GQ.txt
  echo -e "$(bcftools query -R pileup_childhetGQ$GQ.txt -f'[%AD\t]\n' hetnum_childhet_DPGQ20.vcf.gz | \
cut -f3 | sed 's/,/\t/g' | awk '$2/($1 + $2) > 0.35' | wc -l)" >> het_callability_tableGQ$GQ.txt
done
