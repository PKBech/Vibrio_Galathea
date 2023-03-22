FASTQ_DIR='/mnt/raid2/perbec/Vibrio_Galathea'
OUT_DIR='/home/perbec/Vibrios_Galathea/shovill/OUT'
ls $FASTQ_DIR > sample_list.txt
sample_list=$(cat sample_list.txt)

echo $sample_list
ls $FASTQ_DIR

#Run shovill
for stuff in $sample_list 
do 
mkdir $OUT_DIR/"$stuff" 
echo "processing" "$stuff" 
shovill --R1 $FASTQ_DIR/"$stuff"/*1.fq.gz --R2 $FASTQ_DIR/"$stuff"/*2.fq.gz --outdir $OUT_DIR/"$stuff"  --trim --cpus 8 --force; rm -R $FASTQ_DIR/"$stuff" 
done
