# set the primary input files
# these should already be filtered and quality cleaned
# this was the trimmomatic command to filter the HiC data:
# trimmomatic PE -threads 16 CUTRIO-2_S33_L007_R1_001.fastq.gz CUTRIO-2_S33_L007_R2_001.fastq.gz CUTRIO_F_paired.fq.gz CUTRIO_F_unpaired.fq.gz CUTRIO_R_paired.fq.gz CUTRIO_R_unpaired.fq.gz ILLUMINACLIP:/Volumes/ExtremePro/Cuora_genome_data/Cuora_HiC/RawData/CUTRIO-2/TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:4:5 LEADING:5 TRAILING:5 MINLEN:31
# need get_statistics_from_genome_fasta_v2.py in the working directory
# You need filter_five_end.pl, get_stats.pl, and two_read_bam_combiner.pl from the Arima Genomics pipeline in he working directory

######### these variables need your input
TAXON="Cuora"
HIFIBAM="/Volumes/ExtremePro/Cuora_genome_data/Cuora_HiFi/m84071_240531_080822_s2.hifi_reads.bc2014.bam"
HIFI1="/Volumes/ExtremePro/Cuora_genome_data/Cuora_HiFi/Cuora_HiFi.fastq.gz"
HIFI2=""
HICR1="/Volumes/ExtremePro/Cuora_genome_data/Cuora_HiC/RawData/CUTRIO-2/CUTRIO_F_paired.fq.gz"
HICR2="/Volumes/ExtremePro/Cuora_genome_data/Cuora_HiC/RawData/CUTRIO-2/CUTRIO_R_paired.fq.gz"
PER1=""
PER2=""
DIAMONDDB="/Volumes/T9/dbs/uniprot/reference_proteomes.dmnd"			# database for Diamond

######### these variables do not need your input
HIFIOUT=${TAXON}"_HiFi"
HIFIGFA=${HIFIOUT}".bp.p_ctg.gfa"
HIFIFASTA=${HIFIOUT}".bp.p_ctg.fasta"
HIFIFAI=${HIFIOUT}".bp.p_ctg.fasta.fai"
HICOUT=${TAXON}"_HiC"
HICGFA=${HICOUT}".hic.p_ctg.gfa"
HICFASTA=${HICOUT}".hic.p_ctg.fasta"
HICFASTAFAI=${HICOUT}".hic.p_ctg.fasta.fai"
HIC1BAM=${TAXON}"_HiC_R1.bam"
HIC2BAM=${TAXON}"_HiC_R2.bam"
HIC1BAMFILT=${TAXON}"_HiC_R1_filtered.bam"
HIC2BAMFILT=${TAXON}"_HiC_R2_filtered.bam"
HICCOMBINEDBAM=${TAXON}"_HiC_combined.bam"
HICRG=${TAXON}"_HiC_rgadded.bam"
HICDEDUP=${TAXON}"_HiC_dedup.bam"
HICMETRICS=${TAXON}"_HiC_dedup.metrics.txt"
YAHSFASTA=${TAXON}"_yahs.fasta"
YAHSINDEX=${TAXON}"_yahs.fasta.fai"
YAHSALIGNED=${TAXON}"_yahs_hifi_aligned.bam"
YAHSSORTED=${TAXON}"_yahs_hifi_sorted.bam"
PURGEHAPLOGENCOV=${TAXON}"_yahs_hifi_sorted.bam.200.gencov"
PURGEHAPLOCOV=${TAXON}"_coverage.stats"
PEALIGNED=${TAXON}"_PE.sam"
PESORTED=${TAXON}"_PE_sorted.bam"
PEMAPPED=${TAXON}"_PE_mapped.bam"
PEGENCOV=${PEMAPPED}".200.gencov"
PECOVSTATS=${PEMAPPED}"_coverage.stats"
HIFIPAF=${TAXON}"_mimimap2.paf.gz"
HIFISPLIT=${TAXON}"_yahs.split"
HIFISELF=${TAXON}"_yahs.split.paf.gz"
COMPLEASM=${TAXON}"_compleasm"
COMPLEASMFINAL=${TAXON}"_compleasm_final"
BLOB=${TAXON}".blobdir"
DIAMOND=${TAXON}"_yahs.diamond.blastx.out"
DRAFT=${TAXON}"_V1.full.softmasked.fasta"
source "$(conda info --base)/etc/profile.d/conda.sh"

# prep files if they haven't been processed
#conda activate samtools
#samtools fastq -0 ${HIFI1} ${HIFIBAM}
#pigz ${HIFI1}
#conda deactivate

######## HiFi-only assembly ########
#source "$(conda info --base)/etc/profile.d/conda.sh"	# need to activate conda in this shell
#conda activate hifiasm
#mkdir HiFi_only
#cd HiFi_only
#hifiasm -o ${HIFIOUT} -t 12 ${HIFI1} ${HIFI2}
#awk '/^S/{print ">"$2;print $3}' ${HIFIGFA} > ${HIFIFASTA}
#python3 ../get_statistics_from_genome_fasta_v2.py ${HIFIFASTA}
#cd ..

###### purge duplicates before scaffolding ######
######## purge dups ########
#source ~/.start_miniforge3.sh
#source "$(conda info --base)/etc/profile.d/conda.sh"
#mkdir purge_dups
#cd purge_dups
#conda activate purge_dups
# STEP 1
#minimap2 -xasm20 -I50G ../HiFi_only/${HIFIFASTA} ${HIFI1} | gzip -c - > ${HIFIPAF}
#pbcstat ${HIFIPAF}
#calcuts PB.stat > cutoffs 2>calcults.log
# STEP2 - align genome to self
#split_fa ../HiFi_only/${HIFIFASTA} > ${HIFISPLIT}
#minimap2 -xasm5 -DP ${HIFISPLIT} ${HIFISPLIT} | gzip -c - > ${HIFISELF}
# STEP 3 - purge duplicates and overlaps
#purge_dups -2 -T cutoffs -c PB.base.cov ${HIFISELF} > dups.bed 2> purge_dups.log
# STEP 4 - Get purged primary and haplotig sequences from draft assembly.
#get_seqs -e dups.bed ../HiFi_only/${HIFIFASTA}

######## scaffold using yahs ########
source "$(conda info --base)/etc/profile.d/conda.sh"	# need to activate conda in this shell; if already did in shell script, don't have to do again
conda activate bwa-mem2									# samtools and picard should be installed in this environment, if not, install it
#samtools faidx purge_dups/purged.fa
#bwa-mem2 index purge_dups/purged.fa
#mkdir yahs
cd yahs
bwa-mem2 mem -t 10 ../purge_dups/purged.fa ${HICR1} | samtools view -@ 10 -Sb - > ${HIC1BAM}
bwa-mem2 mem -t 10 ../purge_dups/purged.fa ${HICR2} | samtools view -@ 10 -Sb - > ${HIC2BAM}
samtools view -h ${HIC1BAM} | perl ../filter_five_end.pl | samtools view -Sb - > ${HIC1BAMFILT}
samtools view -h ${HIC2BAM} | perl ../filter_five_end.pl | samtools view -Sb - > ${HIC2BAMFILT}
perl ../two_read_bam_combiner.pl ${HIC1BAMFILT} ${HIC2BAMFILT} samtools 10 | samtools view -bS -t ../purge_dups/purged.fa.fai - | samtools sort -@ 10 -o ${HICCOMBINEDBAM} -

source ~/.start_miniconda3.sh				# osx-64
source /Users/mkfujita/miniconda3/etc/profile.d/conda.sh
conda activate picard
picard AddOrReplaceReadGroups \
    I=${HICCOMBINEDBAM} \
    O=${HICRG} \
    RGID=1 \
    RGLB=lib1 \
    RGPL=illumina \
    RGPU=unit1 \
    RGSM=Heterodon
#picard MarkDuplicates INPUT=${HICRG} OUTPUT=${HICDEDUP} METRICS_FILE=${HICMETRICS} ASSUME_SORTED=TRUE VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=TRUE
# the MarkDuplicates didn't work automatically I think because of memory limitations, so had to do it manuall as below
java -Xmx100g -jar /Users/mkfujita/miniconda3/envs/picard/share/picard-3.4.0-0/picard.jar MarkDuplicates INPUT=${HICRG} OUTPUT=${HICDEDUP} METRICS_FILE=${HICMETRICS} ASSUME_SORTED=TRUE VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=TRUE

source "$(conda info --base)/etc/profile.d/conda.sh"	# need to activate conda in this shell; if already did in shell script, don't have to do again
conda activate bwa-mem2									# samtools and picard should be installed in this environment, if not, install it
samtools index ${HICDEDUP}
perl ../get_stats.pl ${HICDEDUP} > ${HICDEDUP}".stats"
conda deactivate
conda activate yahs
yahs ../purge_dups/purged.fa ${HICDEDUP}
mv yahs.out_scaffolds_final.fa ${YAHSFASTA}
#if [ $(stat -f%z ${YAHSFASTA}) -gt 1000000000 ]; then		# if the scaffolded assembly is at least 1 GB, remove intermediate files
#	rm ${HIC1BAMFILT}		# can always recreate these
#	rm ${HIC2BAMFILT}
#	rm ${HICCOMBINEDBAM
#	rm ${HICRG}
#fi
#if [ $(stat -f%z ${HIC1BAMFILT}) -gt 1000000000 ]; then rm ${HIC1BAMFILT}; fi		# remove if at least 1 gigabyte; can always recreate these
#if [ $(stat -f%z ${HIC2BAMFILT}) -gt 1000000000 ]; then rm ${HIC2BAMFILT}; fi
#if [ $(stat -f%z ${HICCOMBINEDBAM}) -gt 1000000000 ]; then rm ${HICCOMBINEDBAM}; fi
#if [ $(stat -f%z ${HICRG}) -gt 1000000000 ]; then rm ${HICRG}; fi
conda activate bwa-mem2
bwa-mem2 index ${YAHSFASTA}				# for aligning PE reads to it later
python3 ../get_statistics_from_genome_fasta_v2.py ${YAHSFASTA}
cd ..

# run compleasm on scaffolded assembly
source ~/.start_miniconda3.sh				# osx-64
source /Users/mkfujita/miniconda3/etc/profile.d/conda.sh
conda activate compleasm
compleasm run -a yahs/${YAHSFASTA} -o ${COMPLEASM} -l sauropsida_odb12 -t 10 -L /Volumes/ExtremePro/tesselatus_PacBio_Revio_April2025/mb_downloads/squamata_odb12

### Had Sam manual curate using Juicer here
# the output fasta is called out_JBAT_final.FINAL.fa
### used the VGP9 decontaminatin pipeline in Galaxy (release v0.7) AFTER Juicer manual curation
### downloaded the files into decontamination_VGP9 folder from Galaxy
### this is the actual draft genome!
# The genome was downloaded as Galaxy25-[gfastats on data 24 and data 12_ edited sequences].fasta
# rename it Cuora_V1.fasta

# RepeatMasking
mkdir repeatmasking
sh repeat_pipeline.sh			# based on Daren Card's blog post

# final compleasm run
conda activate compleasm
compleasm run -a repeatmasking/05_full_out/${DRAFT} -o ${COMPLEASMFINAL} -l sauropsida_odb12 -t 10 -L /Volumes/ExtremePro/tesselatus_PacBio_Revio_April2025/mb_downloads/squamata_odb12

# align using minimap2
conda activate minimap2
mkdir minimap2
cd minimap2
minimap2 -x asm10 -c /Volumes/ExtremePro/Cuora_genome_data/Mauremys_reevesii/ncbi_dataset/data/GCF_016161935.1/Mauremys_reevesii_simplenames.fasta /Volumes/ExtremePro/Cuora_genome_data/repeatmasking/05_full_out/Cuora_V1.full.softmasked.fasta > Mauremys_Cuora_21Sep25.paf
		# the -c tells it to include CIGAR strings in the PAF output, which is required for blastn2dotplots
# view dotplots with Mauremys reevesii using dgenies
# still in minimap2 folder
conda deactivate
conda activate blastn2dotplots
perl paf2blastn-fmt6.pl Mauremys_Cuora_21Sep25.paf > Mauremys_Cuora_21Sep25_blastn-fmt6-like.txt
blastn2dotplots -i1 Mauremy.txt -i2 Cuora.txt --blastn Mauremys_Cuora_21Sep25_blastn-fmt6-like.txt --figure_size 100 100 --font_size 20 --out Mauremys_Cuora_21Sep25_boxplots

# aligned Mauremys reevesi to the Cuora genome using KegAlign in Galaxy. It took several days but was done! Now I have to convert the maf file to a chain file
# installed maftoaxt and axtChain via conda
cd /Volumes/ExtremePro/Cuora_genome_data/lastz
conda activate ucsc-kent
python3 axt_split_from_maf.py Mauremys_Cuora_KegAlign_22Sept25.maf Mauremys_Cuora_KegAlign_22Sept25.axt		# splits the maf into axt files one for each target scaffold because axtChain might give a memory error on the whole axt file
sh axtChain.sh			# chain file commands for all target scaffolds, e.g.: axtChain -linearGap=loose NC_006082.1.axt /Volumes/ExtremePro/Cuora_genome_data/Mauremys_reevesii/ncbi_dataset/data/GCF_016161935.1/Mauremys_reevesii_simplenames.2bit /Volumes/ExtremePro/Cuora_genome_data/repeatmasking/05_full_out/Cuora_V1.full.softmasked.2bit NC_006082.1.chain
sh chainSort.sh			# sort the chain files, e.g.: chainSort NC_006082.1.chain NC_006082.1.sorted.chain
chainMergeSort *sorted.chain > ../Mauremys_Cuora_KegAlign_22Sept25.sorted.chain
./toga.py /Volumes/ExtremePro/Cuora_genome_data/lastz/Mauremys_Cuora_KegAlign_22Sept25.sorted.chain /Volumes/ExtremePro/Cuora_genome_data/Mauremys_reevesii/ncbi_dataset/data/GCF_016161935.1/Maureev.bed /Volumes/ExtremePro/Cuora_genome_data/Mauremys_reevesii/ncbi_dataset/data/GCF_016161935.1/Mauremys_reevesii_simplenames.2bit /Volumes/ExtremePro/Cuora_genome_data/repeatmasking/05_full_out/Cuora_V1.full.softmasked.2bit --pn TOGA_Cuora_V1_KegAlign --cb 3,5 --cjn 500



# liftoff
cd liftoff
conda activate liftoff
liftoff -g genomic.gff -o Cuora_V1.liftoff.GCF_016161935.1.gff3 Cuora_V1.full.softmasked.fasta GCF_016161935.1_ASM1616193v1_genomic.fna
gffread -g Cuora_V1.full.softmasked.fasta -C -V -x Cuora_V1.liftoff.cds.fasta Cuora_V1.liftoff.GCF_016161935.1.gff3


# orthofinder
mkdir orthofinder
conda activate orthofinder
orthofinder -d -f test				# test is a folder with a few cds files from other turtles


# decided to align with Mauremys reevesii again with lastz_32 since I found out an fairly efficient way to do so.
# Then get annotations for Cuora via TOGA!
cd lastz
script_lastz_Maureev1.sh			# see these files for the commands
script_lastz_Maureev2.sh
script_lastz_Maureev3.sh
cat *axt > Mauree_vs_Cuora_17Oct25.axt
conda activate ucsc-kent
axtChain -linearGap=loose Mauree_vs_Cuora_17Oct25.axt /Volumes/ExtremePro/Cuora_genome_data/Mauremys_reevesii/ncbi_dataset/data/GCF_016161935.1/Mauremys_reevesii_simplenames.2bit /Volumes/ExtremePro/Cuora_genome_data/repeatmasking/05_full_out/Cuora_V1.full.softmasked.2bit Mauree_vs_Cuora_17Oct25.chain
chainSort Mauree_vs_Cuora_17Oct25.chain Mauree_vs_Cuora_17Oct25.sorted.chain
# convert the Mauremys gtf file into a BED12 file using kent tools
cd /Volumes/ExtremePro/Cuora_genome_data/Mauremys_reevesii/ncbi_dataset/data/GCF_016161935.1
gtfToGenePred -ignoreGroupsWithoutExons genomic.gtf tmp.genePred
genePredToBed tmp.genePred Maureev.bed
conda deactivate
cd /Users/mkfujita/TOGA
./toga.py /Volumes/ExtremePro/Cuora_genome_data/lastz/Mauree_vs_Cuora_17Oct25.sorted.chain /Volumes/ExtremePro/Cuora_genome_data/Mauremys_reevesii/ncbi_dataset/data/GCF_016161935.1/Maureev.bed /Volumes/ExtremePro/Cuora_genome_data/Mauremys_reevesii/ncbi_dataset/data/GCF_016161935.1/Mauremys_reevesii_simplenames.2bit /Volumes/ExtremePro/Cuora_genome_data/repeatmasking/05_full_out/Cuora_V1.full.softmasked.2bit --pn Cuora_V1_TOGA --cb 3,5 --cjn 500
conda activate ucsc-kent
bedToGenePred query_annotation.bed query_annotation.genePred
genePredToGtf file query_annotation.genePred Cuora_V1.gtf

# Helixer
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate helixer
cd Helixer
Helixer.py --lineage vertebrate \
  --downloaded-model-path "/Users/mkfujita/Helixer/models" \
  --fasta-path ../Cuora_V1.full.softmasked.fasta \
  --species Cuora \
  --gff-output-path Cuora_V1_helixer.gff3

# tiberius
mkdir Tiberius
cd Tiberius
conda activate tiberius			# had to download Tiberius anyway because didn't know how to use it from conda
python tiberius.py --genome /Volumes/ExtremePro/Cuora_genome_data/Cuora_individual_scaffolds/scaffold_10.fasta --out /Volumes/ExtremePro/Cuora_genome_data/Tiberius/scaffold_10_V1_tiberius.gtf --seq_len 259002 --model_cfg model_cfg/mammalia_softmasking_v2.yaml 

# Mikado on TOGA and Helixer gffs, not Tiberius which was taking too long to run
mkdir Mikado
cd Mikado
conda activate Mikado
mikado configure --list Cuora.mikado.cfg.txt --reference ../Cuora_V1.full.softmasked.fasta --mode permissive --scoring mammalian.yaml Cuora.mikado.cfg.yaml
mikado prepare --json-conf Cuora.mikado.cfg.yaml
prodigal -i mikado_prepared.fasta -g 1 -o mikado.orfs.gff3 -f gff
mikado serialise --json-conf Cuora.mikado.cfg.yaml --orfs mikado.orfs.gff3
mikado pick --configuration Cuora.mikado.cfg.yaml --subloci-out mikado.subloci.gff3
python3 rename_genes_gff.py mikado.loci.gff3 Cuotri			# renames the complicated mikado gene names with something simple with the prefix Cuotri
gffread Cuotri.gff -g ../Cuora_V1.full.softmasked.fasta -w Cuotri1.transcripts.fasta   
TransDecoder.LongOrfs -t Cuotri1.transcripts.fasta -m 50
TransDecoder.Predict -t Cuotri1.transcripts.fasta
awk '/^>/{print $1; next} {print}' Cuotri1.transcripts.fasta.transdecoder.pep > Cuotri1.pep.fasta		# simplifies the transdecoder names into something easier for Orthofinder

# orthofinder between Cuora and Chrysemys; not with Mauremys because Helixer is not good with isoforms
# this is to maximize the annotations
cd /Users/mkfujita/orthofinder				# remember you have to run it on internal drive
awk '/^>/{print $1; next} {print}' Chrysemys_picta_bellii-3.0.3.pep.all.fa > Chrpic3.pep.fasta
orthofinder -f Cuora_Chrysemys -t 12 -a 12 -S diamond -d




### start the decontamination process - old process but still ok
# Note this was done after yahs, not after Juicer as was actually done to get the draft genome
# back in the main folder
source "$(conda info --base)/etc/profile.d/conda.sh"
mkdir decontamination
cd decontamination

# remove mitochondrial genome if it exists in the assembly
# do this in Galaxy using MitoHiFi
# if a mitochondrial genome is present, add it to the "scaffolds_to_remove.txt"

# try kraken2
# on MBP M4 Max
cd decontamination
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate kraken2
# on Mac Studio this took 1.5 hours, copy back the kraken2 outputs to ExtremePro and make a list "scaffolds_to_remove.txt" based on the outputs - these are anything that isn't a vertebrate and/or mitochondrial genome
### ****** DETERMINE WHAT THE KRAKEN2 DATABASE IS ***** ######
kraken2 --db /Users/pantai/Desktop/databases/kraken2 --threads 24 --use-names --report kraken2/Cuora.kraken2.report --output kraken2/Cuora.kraken2.out Cuora_yahs.fasta
python3 remove_sequences_from_fasta.py ${YAHSFASTA} scaffolds_to_remove.txt


































# use Diamond and nr.dmnd
# download and prepare the databases. These are held outside of the folder because they are so big
# the nr database needs taxonomy info, so when making the diamond database, make sure you do this
# this process can take 1-2 days
source "$(conda info --base)/etc/profile.d/conda.sh"
mkdir -p ~/dbs/nr_diamond
cd ~/dbs/nr_diamond
wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz
gunzip nr.gz
wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
tar -xvzf taxdump.tar.gz names.dmp nodes.dmp
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz
gunzip prot.accession2taxid.gz
conda activate diamond
diamond makedb --in nr -d nr --taxonmap prot.accession2taxid --taxonnodes nodes.dmp --taxonnames names.dmp

mkdir -p ~/dbs/kraken2
cd ~/dbs/kraken2
conda activate kraken2
kraken2-build --standard --threads 32 --db .
kraken2-build --download-library vertebrate_other --db .
kraken2-build --build --threads 32 --db .
DIAMONDDB = "/Volumes/T9/dbs/nr_diamond/nr.dmnd"	

# configure the Diamond Uniprot database
wget -q -O uniprot/reference_proteomes.tar.gz \
 ftp.ebi.ac.uk/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/$(curl \
     -vs ftp.ebi.ac.uk/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/ 2>&1 | \
     awk '/tar.gz/ {print $9}')


find reference_proteomes -type f -name '*.fasta.gz' -print0 \
  | xargs -0 pigz -dc \
  > uniprot_reference_proteomes.fasta
touch reference_proteomes.fasta.gz
find . -mindepth 2 | grep "fasta.gz" | grep -v 'DNA' | grep -v 'additional' | xargs cat >> reference_proteomes.fasta.gz
python3 make_taxid_map.py			# this goes into the directories and gets the idmappings and gets information from the NCBI_TaxID lines into a tax_id map  
diamond makedb -p 16 --in reference_proteomes.fasta.gz --taxonmap reference_proteomes.taxid_map --taxonnodes /Volumes/T9/dbs/kraken2/newtaxdump/nodes.dmp -d reference_proteomes.dmnd

# back in the main folder
source "$(conda info --base)/etc/profile.d/conda.sh"
mkdir decontamination
cd decontamination
conda activate diamond
diamond blastx -p 14 -q ../yahs/${YAHSFASTA} -d ${DIAMONDDB} --evalue 1e-25 --max-target-seqs 1 --sensitive --outfmt 6 qseqid staxids bitscore qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore > ${DIAMOND}

# map reads back to assembly
# still in decontamination folder
map reads back to assembly
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate minimap2
minimap2 -ax map-hifi -t 14 ../yahs/${YAHSFASTA} ${HIFI1} | samtools sort -@ 8 -o reads_vs_assembly.bam
samtools index reads_vs_assembly.bam

# run blobtools
conda activate blobtools2		# this is blobtools2, installed as instructed here https://blobtoolkit.genomehubs.org/install/

# Create dataset
blobtools2 create --fasta ../yahs/${YAHSFASTA} .

# Add taxonomy from Kraken2
# downloaded new taxdump from ncbi because kraken2 taxdump was insufficient
mkdir /Volumes/T9/dbs/kraken2/ncbi-taxdump
cd /Volumes/T9/dbs/kraken2/ncbi-taxdump
curl -O https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/newtaxdump/new_taxdump.tar.gz
tar -xzf new_taxdump.tar.gz
cd /Volumes/ExtremePro/Cuora_genome_data/decontamination
blobtools add --taxdump /Volumes/T9/dbs/kraken2/newtaxdump --hits kraken2:Cuora.kraken2.out .

# Add coverage from BAM
blobtools add --bam reads_vs_assembly.bam .

# View interactively
blobtools view --local my_dataset.blobdir











# configure the Diamond database
wget -q -O uniprot/reference_proteomes.tar.gz \
 ftp.ebi.ac.uk/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/$(curl \
     -vs ftp.ebi.ac.uk/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/ 2>&1 | \
     awk '/tar.gz/ {print $9}')


find reference_proteomes -type f -name '*.fasta.gz' -print0 \
  | xargs -0 pigz -dc \
  > uniprot_reference_proteomes.fasta
touch reference_proteomes.fasta.gz
find . -mindepth 2 | grep "fasta.gz" | grep -v 'DNA' | grep -v 'additional' | xargs cat >> reference_proteomes.fasta.gz
python3 make_taxid_map.py			# this goes into the directories and gets the idmappings and gets information from the NCBI_TaxID lines into a tax_id map  
diamond makedb -p 16 --in reference_proteomes.fasta.gz --taxonmap reference_proteomes.taxid_map --taxonnodes /Volumes/T9/dbs/kraken2/newtaxdump/nodes.dmp -d reference_proteomes.dmnd
cd -




find . -mindepth 2 -name "*.idmapping.gz" -print0 | xargs -0 zcat | grep "NCBI_TaxID" | awk '{print $1 "\t" $1 "\t" $3 "\t" 0}' >> reference_proteomes.taxid_map





######## prepare scaffolded genome for Juicer ########			
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate samtools
mkdir Juicer
chdir Juicer

### index assembly ###
samtools faidx ../purge_dups/purged.fa			# use purged.fa because this is what was inputted into yahs

### Make chromosome size file ###
cut -f1,2 ../yahs/${YAHSINDEX} > ../yahs/scaffolds_final.chrom.sizes

### run juicer ###
conda deactivate
conda activate yahs
juicer pre -a -o out_JBAT ../yahs/yahs.out.bin ../yahs/yahs.out_scaffolds_final.agp ../purge_dups/purged.fa.fai > out_JBAT.log 2>&1

cat out_JBAT.log | grep PRE_C_SIZE | awk '{print $2" "$3}'


### make .hic ###
java -Xmx128G -jar /Users/mkfujita/bin/juicer_tools_1.22.01.jar pre out_JBAT.txt out_JBAT.hic <(cat out_JBAT.log | grep PRE_C_SIZE | awk '{print S2" "S3}')
java -Xmx128G -jar /Users/mkfujita/bin/juicer_tools_1.22.01.jar pre out_JBAT.txt out_JBAT.hic <(echo "assembly 1104241180")