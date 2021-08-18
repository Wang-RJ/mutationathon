# We don't care if it's really homozygote. Just if it's called homozygote.

sample=$(head -n 1 ../triotable.txt | tail -n 1)
bcftools view -Ou -M2 -m2 -i'INFO/AC < 2' -s $sample ../vcf/output.filtered.snps.removed.CalculateGenotypePosteriors.vcf.gz | \
bcftools annotate -x INFO/CSQ,INFO/InbreedingCoeff,INFO/BaseQRankSum,INFO/ExcessHet,INFO/ClippingRankSum,INFO/FS,INFO/MLEAC,INFO/MLEAF,INFO/MQRankSum,INFO/ReadPosRankSum,INFO/SOR \
-i'(GT[0]="RR" && GT[1]="RR")' | awk 'BEGIN {srand()} (rand() <= 0.22 || /^#/ )' | bcftools view -Oz -o homocandidates.vcf.gz &

bcftools view -Oz -i'FORMAT/DP[0] > 20 && FORMAT/DP[1] > 20 && FORMAT/DP[0] < 80 && FORMAT/DP[1] < 80' homocandidates.vcf.gz -o homocandidates.filtered.vcf.gz &

mombam=$(cut -d',' -f1 ../triotable.txt | head -n 1 | tail -n 1).dups.bam
dadbam=$(cut -d',' -f2 ../triotable.txt | head -n 1 | tail -n 1).dups.bam
bcftools mpileup -Oz -A -f /N/dcwan/projects/hahnlab-phi/macaque/ref/ref-relabeled.fna -Q 10 \
-R <(bcftools query homocandidates.filtered.vcf.gz -f'%CHROM\t%POS\n') \
-a ADF,ADR,AD,DP ../bam_files/$mombam ../bam_files/$dadbam > homo_pileup.vcf.gz &

bcftools index homo_pileup.vcf.gz
bcftools index homocandidates.filtered.vcf.gz

for GQ in {20,30,40,50}
do
  bcftools query homo_pileup.vcf.gz -R <(bcftools query -i'FORMAT/GQ[0] > '$GQ' && FORMAT/GQ[1] > '$GQ homocandidates.filtered.vcf.gz -f '%CHROM\t%POS\n') -e'TYPE="snp" && FORMAT/AD[*:1] > 0' -f'%CHROM\t%POS\n' > pileup_homoGQ$GQ.txt &
done

for GQ in {20,30,40,50}
do
  header="trio\tdenominator\tGQ_filtered\tbam_filtered"
  echo -e $header > homo_callabilityraw_tableGQ$GQ.txt
  i=1
  echo -e "$i\t\c" >> homo_callabilityraw_tableGQ$GQ.txt
  echo -e "$(zcat homocandidates.filtered.vcf.gz | grep -v '#' | wc -l)\t\c" >> homo_callabilityraw_tableGQ$GQ.txt
  echo -e "$(bcftools query -i'FORMAT/GQ[0] > '$GQ' && FORMAT/GQ[1] > '$GQ homocandidates.filtered.vcf.gz -f '%CHROM\t%POS\n' | wc -l)\t\c" >> homo_callabilityraw_tableGQ$GQ.txt
  echo -e "$(cat pileup_homoGQ$GQ.txt | wc -l)" >> homo_callability_tableGQ$GQ.txt
done
