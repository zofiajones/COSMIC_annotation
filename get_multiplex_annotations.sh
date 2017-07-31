#get mutations from mysql database
python browse_database_3.py

#prep COSMIC vcf as bed file
grep -v '#' CosmicCodingMuts.vcf | python parse_COSMIC.py > CosmicCoding_1.bed
~/hb/tools/bedtools sort -i CosmicCoding_1.bed > CosmicCoding_sort.bed
bgzip -c CosmicCoding_sort.bed > CosmicCoding_sort.bed.gz
~/hb/tools/tabix CosmicCoding_sort.bed.gz



cp ~/Mul*/full.vcf ~/COSMIC/
#get exact matches where genomic location is known
cat full.vcf | python get_COSMIC_1.py | sort | uniq > cell_line_COSMIC.txt
#get matches where genomic location is know - location correct, amino acid change, base change not
cat exceptions.txt | grep -v Intron | grep -v UTR | grep -v splice | python get_COSMIC_2.py > mix_match.txt
cat ~/Mul*/mutations_1.csv | grep -v protein > ~/COSMIC/variant_info.txt
#get possible matches where genomic location is not known
bash search_COSMIC_names.sh
awk 'BEGIN{OFS="\t"}{print $1,$2,$4,$5,$8}' COSMIC_hg19.vcf | sed 's/p.//g' > COSMIC_hg19_1.vcf

#collect annotation data
cat cell_line_COSMIC.txt > multiplex_annotations_1.txt
cat COSMIC_hg19_1.vcf >> multiplex_annotations_1.txt
cat mix_match.txt >> multiplex_annotations_1.txt

#annotate mutations data
sort multiplex_annotations_1.txt | uniq > multiplex_annotations.txt
cat multiplex_annotations.txt | wc -l
python match_mutations.py
grep -v UTR mutations_3.csv | grep -v Intron | grep -v splice | grep -v protein | sed 's/fs//g' > variant_info.txt
#manually check unannotated mutations
bash search_COSMIC_names.sh


#check and query results
cat mix_match.txt | awk '{print $5}' | awk -F ';' '{print $1,$3}' > info.txt

while read p;do
echo $p
array=($p)
grep ${array[0]} COSMIC_hg19.vcf | grep ${array[1]}
done < info.txt


while read p;do

q=$(echo $p | awk '{print $8}' | awk -F ';' '{print $1,$3}' | sed 's/p.//g')
array=($q)

r=$(grep ${array[0]} info.txt | grep ${array[1]} | wc -l)

if [[ $r == 0  ]];then
echo $p
fi

done < COSMIC_hg19.vcf

cat mutations_2.csv | awk -F ',' '{print $6}' | grep -v eng | grep -v 'NA' | awk -F ';' '{print NF-1}' | awk '{total = total + $1}END{print total}'
cat ~/Mul*/cell_lines.csv | grep -v id | grep -v '100%' | awk -F ',' '{print $5}' > 2.txt
cat mutations_5.csv | grep -v cell | awk -F ',' '{print $2}' > 1.txt
grep -vf 1.txt 2.txt
