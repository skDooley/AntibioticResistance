#!/bin/bash -login
: '
to create an environment to run this pipeline using conda you can use the commands below
conda create -n microbe -c bioconda khmer,trimmomatic,samtools,bowtie2,microbecensus
conda activate microbe
This script also assumes that anaconda3 is installed at $USER/bin/anaconda3/bin. If it is not, 
change MappingHeader.sb to the location of your main conda install. If you already of some of these
tools installed, ensure all the packages above are installed and MappingHeader.sb will say "unable to activate microbe"
but everything will still run on your conda.
'

#################################### Paths for storage of files ####################################
SAMPLE_TYPE="treated"
ANTIBIO=/mnt/research/germs/shane/antibioticResistance
RLOGS=$ANTIBIO/data/run_logs
SCRATCH_DIR=$SCRATCH/antibioticResistance

#fastqs dirs
TRIMMED=$SCRATCH_DIR/trimmed/$SAMPLE_TYPE
UNPAIRED=$SCRATCH_DIR/raw/$SAMPLE_TYPE

#Executables dir
SAMPLE_SCRIPTS=$SCRATCH_DIR/scripts
TRIMSTATS=$ANTIBIO/data/stats

#################################### House keeping run Variables ####################################
HEADER=$ANTIBIO/scripts/hpc_scripts/MappingHeader.sb 
cd $UNPAIRED  #Starting Location
THREADS=20    ## of threads for trimming and alignment
files=(*.fastq)   #Get the sample files to process
nsamples=${#files[@]} #Number of samples to process
counter=0 #Counter to keep track of what sample we are processing

#If the dirs for output don't exist create them
mkdir -p $SCRATCH_DIR $TRIMMED $SAMPLE_SCRIPTS $TRIMSTATS
cd $SCRATCH_DIR

#################################### For Each Sample Build HPC Script and Launch ####################################
for fastq in "${files[@]}"; do 
	counter=$((counter + 1))
	sample=${fastq/\.fastq/}

	cat $HEADER >$SAMPLE_SCRIPTS/$sample.sb
	echo -en "cd $SCRATCH_DIR\n\n" >> $SAMPLE_SCRIPTS/$sample.sb
	echo "echo \"$sample\""  >> $SAMPLE_SCRIPTS/$sample.sb


	# #################################### Step 3. Trim adapters and QC reads #################################### 
	echo -en "trimmomatic SE -phred33 -threads $THREADS $UNPAIRED/$sample.fastq $TRIMMED/$sample.fastq.gz ILLUMINACLIP:TruSeq2-SE.fa:2:30:10:8:TRUE LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 2>$TRIMSTATS/$sample.$SAMPLE_TYPE.log \n\n" >>$SAMPLE_SCRIPTS/$sample.sb
	echo "echo \"Completed Trimming Reads\"" >>$SAMPLE_SCRIPTS/$sample.sb
	echo -en "$counter $sample Job: "

	# If it's the last sample to process, then add an email flag so I know when it is done
	if [ $counter = $nsamples ]; then
		#echo "$SAMPLE_SCRIPTS/$sample.sb"
		sbatch --time=04:00:00 --mem=50G --cpus-per-task=20 --job-name=$sample -e "run_logs/$sample.$SAMPLE_TYPE.err" -o "run_logs/$sample.$SAMPLE_TYPE.out" --mail-type=ALL --mail-user=dooley.shanek@gmail.com "$SAMPLE_SCRIPTS/$sample.sb"
		echo "Done"
	else 
		#echo "$SAMPLE_SCRIPTS/$sample.sb"
		sbatch --time=04:00:00 --mem=50G --cpus-per-task=20 --job-name=$sample -e "run_logs/$sample.$SAMPLE_TYPE.err" -o "run_logs/$sample.$SAMPLE_TYPE.out" "$SAMPLE_SCRIPTS/$sample.sb"
	fi
done

SAMPLE_TYPE="untreated"
ANTIBIO=/mnt/research/germs/shane/antibioticResistance
RLOGS=$ANTIBIO/data/run_logs
SCRATCH_DIR=$SCRATCH/antibioticResistance

#fastqs dirs
TRIMMED=$SCRATCH_DIR/trimmed/$SAMPLE_TYPE
UNPAIRED=$SCRATCH_DIR/raw/$SAMPLE_TYPE

#Executables dir
SAMPLE_SCRIPTS=$SCRATCH_DIR/scripts

TRIMSTATS=$ANTIBIO/data/stats

#################################### House keeping run Variables ####################################
HEADER=$ANTIBIO/scripts/hpc_scripts/MappingHeader.sb 
cd $UNPAIRED  #Starting Location

THREADS=20    ## of threads for trimming and alignment


files=(*.fastq)   #Get the sample files to process
nsamples=${#files[@]} #Number of samples to process
counter=0 #Counter to keep track of what sample we are processing

#If the dirs for output don't exist create them
mkdir -p $SCRATCH_DIR $TRIMMED $SAMPLE_SCRIPTS $TRIMSTATS

cd $SCRATCH_DIR

#################################### For Each Sample Build HPC Script and Launch ####################################
for fastq in "${files[@]}"; do 
	counter=$((counter + 1))
	sample=${fastq/\.fastq/}

	cat $HEADER >$SAMPLE_SCRIPTS/$sample.sb
	echo -en "cd $SCRATCH_DIR\n\n" >> $SAMPLE_SCRIPTS/$sample.sb
	echo "echo \"$sample\""  >> $SAMPLE_SCRIPTS/$sample.sb


	# #################################### Step 3. Trim adapters and QC reads #################################### 
	 echo -en "trimmomatic SE -phred33 -threads $THREADS $UNPAIRED/$sample.fastq $TRIMMED/$sample.fastq.gz ILLUMINACLIP:TruSeq2-SE.fa:2:30:10:8:TRUE LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 2>$TRIMSTATS/$sample.E.log \n\n" >>$SAMPLE_SCRIPTS/$sample.sb
	# echo -en "cat $TRIMMED/$sample.fastq.se1.gz $TRIMMED/$sample.fastq.se2.gz > $TRIMMED/$sample.fastq.se12.gz\n #rm $TRIMMED/$sample.fastq.se1.gz $TRIMMED/$sample.fastq.se2.gz\n\n" >>$SAMPLE_SCRIPTS/$sample.sb
	echo "echo \"Completed Trimming Reads\"" >>$SAMPLE_SCRIPTS/$sample.sb
	
	
	# #################################### Step 4. Remove the host related reads #################################### 
	# #4a) bowtie2 mapping against host sequence
	# if [[ $sample == *"G5"* ]]; then
 # 		SAMFILE=$SAMS/$sample.SWGRASS.sam
 # 		BAMFILE=$BAMS/$sample.SWGRASS.bam
 # 		echo -en "bowtie2 --threads $THREADS -x $SWITCHGRASS -1 $TRIMMED/$sample.fastq.pe1.gz -2 $TRIMMED/$sample.fastq.pe2.gz -U $TRIMMED/$sample.fastq.se12.gz -S $SAMS/$sample.SWGRASS.sam >$CLEANING_STATS/$sample.stat 2>&1 \n\n" >> $SAMPLE_SCRIPTS/$sample.sb
 # 	else
 # 		SAMFILE=$SAMS/$sample.MISCANTS.sam
 # 		BAMFILE=$BAMS/$sample.MISCANTS.bam
 # 		echo -en "bowtie2 --threads $THREADS -x $MISCANTHUS -1 $TRIMMED/$sample.fastq.pe1.gz -2 $TRIMMED/$sample.fastq.pe2.gz -U $TRIMMED/$sample.fastq.se12.gz -S $SAMS/$sample.MISCANTS.sam >$CLEANING_STATS/$sample.stat 2>&1 \n\n" >> $SAMPLE_SCRIPTS/$sample.sb
	# fi
	# echo "samtools view -bS $SAMFILE > $BAMFILE" >>$SAMPLE_SCRIPTS/$sample.sb
	# echo "echo \"Completed Alignment\""          >>$SAMPLE_SCRIPTS/$sample.sb

	# #4b) filter required unmapped reads
	# echo "samtools view -b -f 12 -F 256 $BAMFILE > $BAMS/$sample.unmapped.bam"  >>$SAMPLE_SCRIPTS/$sample.sb

	# #4c) split paired-end reads into separated fastq files .._R1 .._R2
	# echo "samtools sort -n $BAMS/$sample.unmapped.bam -o $BAMS/$sample.unmapped_sorted.bam"  >>$SAMPLE_SCRIPTS/$sample.sb
	# echo "bedtools bamtofastq -i $BAMS/$sample.unmapped_sorted.bam -fq $CLEANED_FASTQS/$sample.R1.fastq -fq2 $CLEANED_FASTQS/$sample.R2.fastq"  >>$SAMPLE_SCRIPTS/$sample.sb
	#Finally: Launch the script on the hpc

	echo -en "$counter $sample Job: "

	# If it's the last sample to process, then add an email flag so I know when it is done
	if [ $counter = $nsamples ]; then
		#echo "$SAMPLE_SCRIPTS/$sample.sb"
		sbatch --time=04:00:00 --mem=40G --cpus-per-task=20 --job-name=$sample -e "run_logs/$sample.$SAMPLE_TYPE.err" -o "run_logs/$sample.$SAMPLE_TYPE.out" --mail-type=ALL --mail-user=dooley.shanek@gmail.com "$SAMPLE_SCRIPTS/$sample.sb"
		echo "Done"
	else 
		#echo "$SAMPLE_SCRIPTS/$sample.sb"
		sbatch --time=04:00:00 --mem=40G --cpus-per-task=20 --job-name=$sample -e "run_logs/$sample.$SAMPLE_TYPE.err" -o "run_logs/$sample.$SAMPLE_TYPE.out" "$SAMPLE_SCRIPTS/$sample.sb"
	fi
done



#################################### Alternative way to do the count files ####################################
#module load bedtools
#python coverage-bed-reference.py $CONTIGS_FILE > $CONTIGS_FILE.bed

#for x in mapping-data/*sorted; do echo "bamToBed -i $x > $x.bed"; done > bamtobed.sh
#cat bamtobed.sh | $PAR_PATH/parallel
#for x in mapping-data/*bed; do echo "coverageBed -a $CONTIGS_FILE.bed -b $x -d > $x.bed2"; done > coveragebed.sh
#cat coveragebed.sh | $PAR_PATH/parallel
#for x in mapping-data/*bed2; do echo "python bedcoverage-to-coverage.py $x > $x.counts"; done > bedcoveragefinal.sh
#cat bedcoveragefinal.sh | $PAR_PATH/parallel
