#cat ~/MultiPlex/aa_COSMIC.txt | sort | uniq  > variant_info.txt

##cat ~/custom_3/MSI.txt | sort | uniq > variant_info.txt

#cat all_mutations_1.txt | sort | uniq > variant_info.txt

#cp ATCC_variants.txt variant_info.txt

#cp ~/MultiPlex/par_no_location.txt variant_info.txt

#sed 's/fs//g'  /Seagate/CRUKCI_vcfs/parental_5a.txt > variant_info.txt

#cp ~/HD753_variants/aa_COSMIC.txt variant_info.txt

#cp ~/BRCA_HD795/COSMIC_aa.txt variant_info.txt

rm out_variants.txt
rm COSMIC_hg19.vcf

while read p;do

array=($p)
echo ${array[0]}

##grep ${array[0]} ~/COSMIC/CosmicMutantExport.tsv | grep -F ${array[1]} | sort | uniq | awk -F '\t' 'BEGIN{OFS="\t"}{print $1,$17,$18,$19,$24,$25}' >> out_variants.txt
#awk -F '\t' -v x=${array[0]} '{if($1==x){print $0}}' ~/COSMIC/CosmicMutantExport.tsv | grep -F ${array[1]} | sort | uniq | awk -F '\t' 'BEGIN{OFS="\t"}{print $1,$17,$18,$19,$24,$25,$5}' >> out_variants.txt

grep -F ${array[0]} ~/COSMIC/CosmicMutantExport.tsv | grep -F ${array[1]} | sort | uniq | awk -F '\t' 'BEGIN{OFS="\t"}{print $1,$17,$18,$19,$24,$25,$5}' | sed 's/_ENS.*//g' >> out_variants.txt

###| awk -v x=$p '{print $0,x}'   >> out_variants.txt

##done < ~/HD701_end*/variant_info.txt
##done < ~/BRCA_variants/variant_list_1.txt
##done < to_find_info.txt
done < variant_info.txt

awk -v x="NA" 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$5,$6,x}' out_variants.txt | sort | uniq > out_variants_1.txt

cat  out_variants_1.txt | sort | uniq | python parse_loc.py

awk 'BEGIN{OFS="\t"}{print $1,$2-1,$2,$3,$4,$5,$6,$7,$8}' COSMIC.vcf > COSMIC.bed

rm alt_1.bed
while read p;do

array=($p)
echo -e "chr"${array[0]}"\t"${array[1]}"\t"${array[2]}  > line.bed
echo -e ${array[3]}"\t"${array[4]}"\t"${array[5]}"\t"${array[6]}"\t"${array[7]}"\t"${array[8]}  > info.bed
~/bin/liftOver line.bed ~/bin/hg38ToHg19.over.chain.gz output.bed unlifted.bed
if [[ $(cat output.bed | wc -l) > 0  ]];then
paste -d '\t'  output.bed info.bed  | awk 'BEGIN{OFS="\t"}{print $1,$3,$4,$5,$6,$7,$8,$9}' >> COSMIC_hg19.vcf
fi

done < COSMIC.bed

