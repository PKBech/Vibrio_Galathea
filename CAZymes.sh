##Run dbCAN CAZymes profiling


mkdir amino_acid_seqs
AAs='/Users/Pernille/Documents/Vibrio_Galathea/anvio/amino_acid_seqs' 

#Extract aa seqs from anvio dbs 
for stuff in $sample_list 
do 
anvi-get-sequences-for-gene-calls -c $CONTIG_dbs/"$stuff"-contigs.db -o $AAs/"$stuff".faa --get-aa-sequences
done

#Make hmm CAZyme database
hmmpress dbCAN-fam-HMMs.txt.v10.txt # convert file to hmmpress

mkdir hmms_Vibrios
HMMS='/Users/Pernille/Documents/Vibrio_Galathea/anvio/hmms_Vibrios' 

#Run HMM search
for stuff in $sample_list 
do 
echo "converting contig names to fixed names in" "$stuff" 
hmmscan --domtblout $stuff.out.dm ../dbCAN-fam-HMMs.txt.v10.txt $AAs/"$stuff".faa 
done

#Clean up text file and import to R for further summarize

cd hmms_Vibrios

CON='/Users/Pernille/Documents/Vibrio_Galathea/contig-db-skesa'

ls $CON | sed 's/_contigs.fa//g' > sample_list.txt
sample_list="$(cat sample_list.txt)"



echo $sample_list

for stuff in $sample_list 
do
head -n $(($(wc -l < $stuff.out.dm) - 10)) $stuff.out.dm | sed '1d;2d;3d' | tr -s ' ' | sed 's/ /\t/g' | cut -f 1-22 > $stuff.out.dm.txt
done