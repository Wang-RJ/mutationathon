LC_ALL=C sort -V M_depth.txt > M_depth.sorted.txt &
LC_ALL=C sort -V Heineken_depth.txt > Heineken_depth.sorted.txt &
LC_ALL=C sort -V Noot_depth.txt > Noot_depth.sorted.txt &

join -j 2 -o 1.1,1.2,1.3,2.3 M_depth.sorted.txt Heineken_depth.sorted.txt > MH_depth.sorted.txt &

# join -j 2 -o 1.1,1.2,1.3,1.4,2.3 MN_depth.sorted.txt Heineken_depth.sorted.txt > joint_pileup.sorted.txt &

join -j 1 <(awk 'BEGIN {OFS = "\t"} {print $1":"$2,$3}' M_depth.sorted.txt) <(awk 'BEGIN {OFS = "\t"} {print $1":"$2,$3}' Heineken_depth.sorted.txt) > MH_depth.sorted.txt &

awk 'NR==FNR {h[$1,$2] = $3; next} {print $1, $2, $3, h[$1,$2]}' Heineken_depth.txt M_depth.txt > MH_depth.txt

samtools depth ../bam_files/M.dups.bam ../bam_files/Noot.dups.bam ../bam_files/Heineken.dups.bam > MNH_depthB.txt

awk '{
	if($3 > 20 && $3 < 80 &&
       $4 > 20 && $4 < 80 &&
       $5 > 20 && $5 < 80) count[$1]++
	   }
  END {
    for(key in count) { print key,count[key] }
  }
' joint_pileup.txt > pile_bychr.txt

# grep "chr[1-9]" M_depth_sample.txt | wc -l
# 8343
LC_ALL=C grep -Ff <(grep "chr[1-9]" M_depth_sample.txt | cut -f1,2) Noot_depth.txt > MN_samplematch.txt

bedtools genomecov -ibam bam_files/M.dups.bam -bga > pile_depth/M_gcov_out.txt 


grep ^COV M.dups.stats | awk '{if($3 > 20 && $3 < 80) print $4}' | paste -sd+ - | bc
# 2771898385
grep ^COV M.dups.stats | cut -f4 | paste -sd+ - | bc
# 3083354638

grep ^COV Heineken.dups.stats | awk '{if($3 > 20 && $3 < 80) print $4}' | paste -sd+ - | bc
# 2097338957
grep ^COV Heineken.dups.stats | cut -f4 | paste -sd+ - | bc
# 3078454313

grep ^COV Noot.dups.stats | awk '{if($3 > 20 && $3 < 80) print $4}' | paste -sd+ - | bc
# 2192899065
grep ^COV Noot.dups.stats | cut -f4 | paste -sd+ - | bc
# 3079842150
