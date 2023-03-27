# Quality control reports. Before alignment happens we need to check out our data to make sure it looks okay.

# Drag and drop your sample file data folder into the terminal or copy the file path and use.
#cd /your/file/path/here/

# Make sure scripts are in the same directory and the correct permissions are given. 

for R1_SRC in *R1.fastq.gz*
do
    DIRECTORY=${R1_SRC//R1.fastq.gz/output}
    eval mkdir -v -m 777 $DIRECTORY
done 

# Make sure youre in the trim environment that has fastqc, cutadapt, trim-galore and fq installed.
source activate trim

# Validation of paired end data using fq. Toggle off with # if using single data.
# Exits cleanly if no error. R1 and R2 being first and second or left and right read files of paired data. Must be .gz file.
 #fq lint  R1_SRC R2_SRC 

 # Looped Pair validation. If pairs are named the same except for the ending of either R1.fastq.gz or R2.fastq.gz then it might be faster to use a loop. 
 # For this we will create a parameter R1_SRC that can be tied to a list, *R1.fastq.gz* which is comprised of all files with that ending in the directory.
 # From R1_SRC a second parameter R2_SRC will be created by searching for R1 in the file names and replacing them with R2.
 for R1_SRC in *R1.fastq.gz*
 do
     R2_SRC=${R1//R1.fastq.gz/R2.fastq.gz}
    eval fq lint $R1_SRC $R2_SRC
done

echo "Paired end files validated"

# For a one sample fastq file run (toggle command on an off by removing/adding before running.
# Multithreading maybe applied to your command by adding -t #of threads. Here it is set to 8, if your computer has more threads then you can increase number for faster output.
#fastq yoursample.fastq.gz -t 8

# Creating a loop for fastqc. Realistically, youre probably going to have multiple samples to run. Creating a loop will automate the process. 
# For this we will create a variable FASTQ that can be tied to a list, *fastq.gz* which is comprised of all files with that ending in the directory.
# Multithreading maybe applied to your command by adding -t #of threads. Here it is set to 8, if your computer has more threads then you can increase number for faster output.
for R1_SRC in *R1.fastq.gz*
do
    DIRECTORY=${R1_SRC//R1.fastq.gz/output}
    R2_SRC=${R1_SRC//R1.fastq.gz/R2.fastq.gz}
	eval fastqc $R1_SRC -t 8 --outdir $DIRECTORY
    eval fastqc $R2_SRC -t 8 --outdir $DIRECTORY
done 

# Optional-Trimming your reads. Sometimes, near the end of a long fragment you may get poor quality scores for bases. Most aligners do soft trimmings of these sequences.
# Kallisto is a pseudoaligner and doesnt soft trim. It does not need to trim adaptor sequences. 
#However, if you notice your base quality decreases near the end then trim. No more than 80%.

# trim-galore works with cutadapt and fastqc to sense adaptor sequences and poor quality bases to do trimmings.
# This loop will perform trimmings on all files ending with fastq.gz and then generate: a new trimmed fastq.gz file, a trimming report and will do another fastqc report for the trimmed data.
for R1_SRC in *R1.fastq.gz*
do
    DIRECTORY=${R1_SRC//R1.fastq.gz/output}
    R2_SRC=${R1_SRC//R1.fastq.gz/R2.fastq.gz}
	eval trim_galore --fastqc --paired $R1_SRC $R2_SRC --clip_R1 23 --clip_R2 23  -o $DIRECTORY
    
done 

echo "QCReports have been generated"

conda deactivate
# For Bulk RNA-seq make sure youre in the rnaseq environment that has kallisto, samtools and multiqc installed.
source activate bioinformatics

# Paired Bulk RNA-seq. This loop will pseudoalign samples to the refenence index. It will create output folders for each sample along with a log folder of the run.
# Multithreading maybe applied to your command by adding -t #of threads. Here it is set to 8, if your computer has more threads then you can increase number for faster output.
#for R1_SRC in $(find . -type f -name '*R1_trimmed.fq.gz'); do
#	R2_SRC=${R1_SRC//R1_trimmed.fq.gz/R2_trimmed.fq.gz}
#	OUT="${R1_SRC%R1_trimmed.fq.gz}_mapped"
#	LOG="${R1_SRC%R1_trimmed.fq.gz}_mapped.log"
#	eval kallisto quant -i Homo_sapiens.GRCh38.cdna.all.index -o "$OUT" "$R1_SRC" "$R2_SRC" -t 8 &> "$LOG"
#done

#echo "Pseudoalignment has been completed"

# Optional-For paired end reads (Single Cell and Bulk) you may also want to create a BAM file. 
# kallisto command --pseudobam allows kallisto data to be turned into a pseudobam file. Using "|" the pseudobam file can be piped into samtools to be turned into a BAM file.
for R1_SRC in $(find . -type f -name '*R1_val_1.fq.gz'); do
	R2_SRC=${R1_SRC//R1_val_1.fq.gz/R2_val_2.fq.gz}
	OUT="${R1_SRC%R1_val_1.fq.gz}_mapped"
	LOG="${R1_SRC%R1_val_1.fq.gz}_mapped.log"
	eval kallisto quant -i Homo_sapiens.GRCh38.cdna.all.index -o "$OUT" "$R1_SRC" "$R2_SRC" -t 8 &> "$LOG" 
done

#eval kallisto quant -i Homo_sapiens.GRCh38.cdna.all.index -o "$OUT" "$R1_SRC" "$R2_SRC" -t 8 --rf- stranded--pseudobam &> "$LOG" 
#echo "BAM file has been created"

#kallisto quant -i Homo_sapiens.GRCh38.cdna.all.index -o AMH01_iPSCs_A AMH01_iPSCs_A_S15_L001_R1_val_1.fq.gz AMH01_iPSCs_A_S15_L001_R2_val_2.fq.gz -t 8 --fr-stranded --pseudobam 
#kallisto quant -i Homo_sapiens.GRCh38.cdna.all.index -o AMH01_iPSCs_A AMH01_iPSCs_A_S15_L001_R1_val_1.fq.gz AMH01_iPSCs_A_S15_L001_R2_val_2.fq.gz -t 8 --rf-stranded --pseudobam 
conda deactivate
source activate rnaseq

# Compiled Quality control. Using MultiQC you can compile all the fastqc reports and log folders from quantifications together. 
# If data was trimmed, move old fastqc reports out of directory before running.
# Running multiqc will generate an html file that can be opened in your browser to view your data.
multiqc --dirs .

echo "multiqc report has been generated"

#featureCounts -p -s 1 -a Homo_sapiens.GRCh38.109.gtf.gz -o geneCount.csv *.bam

#gunzip -c AMH04_iPSCs_ECs_B_S28_L001_R1.fastq.gz | wc -l