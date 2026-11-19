# see Daren Card's great tutorial on how to do this, plus the pipeline is based on his tutorial
# https://darencard.net/blog/2022-07-09-genome-repeat-annotation/

# assumes that RepeatModeler and RepeatMasker are installed with all dependencies

# make sure the dfam library is installed
# for conda via miniconda, they can be found: /Users/mkfujita/miniconda3/envs/repeatmodeler/share/RepeatMasker/Libraries/famdb
# dfam39_full.12.h5 should include Lepidodosauria
# also have Daren Card's script rmOutToGFF3custom in the working directory
# and your script, RepeatModeler_process.py, in the working directory
# have bedtools installed in the environment

# be sure to be working in x86_64 environment ("intel")
# activate the repeatmodeler conda environment (conda activate repeatmodeler)

###### the only input you need to provide
export GENOME="Cuora_V1.fasta"
export GENOME_BLASTDB="Cuora_V1"
export ID="Cuotri1"

###### Create output directories
#mkdir -p logs 01_simple_out 02_tetrapoda_out 03_known_out 04_unknown_out

###### RepeatModeler
export RMFASTA=${GENOME_BLASTDB}"-families.fa"
#BuildDatabase -name ${GENOME_BLASTDB} ${GENOME}
#RepeatModeler -threads 16 -engine ncbi -database ${GENOME_BLASTDB} 2>&1 | tee 00_repeatmodeler.log
#python3 RepeatModeler_process.py ${RMFASTA} ${ID}		# this is a script that I wrote and needs to be in the directory
														# produces two files, "known" and "unknown" repeats

###### Mask simple repeats
# Note RepeatMasker is part of the repeatmodeler conda environment
#RepeatMasker -pa 16 -a -e ncbi -species vertebrata -dir 01_simple_out -noint -xsmall ${GENOME} 2>&1 | tee logs/01_simplemask.log
# rename files to something meaningful
# zsh does not include the rename function so use these commands
# these are the output files
export ALIGN=${GENOME}.align		
export CAT=${GENOME}.cat.gz
export MASKED=${GENOME}.masked
export OUT=${GENOME}.out
export TBL=${GENOME}.tbl
export NEWMASKED=${MASKED/fasta/simple_mask}
echo ${${MASKED}/fasta/simple_mask}
mv 01_simple_out/${ALIGN} 01_simple_out/${ALIGN/fasta/simple_mask}
mv 01_simple_out/${CAT} 01_simple_out/${CAT/fasta/simple_mask}
mv 01_simple_out/${MASKED} 01_simple_out/${NEWMASKED}.fasta
mv 01_simple_out/${OUT} 01_simple_out/${OUT/fasta/simple_mask}
mv 01_simple_out/${TBL} 01_simple_out/${TBL/fasta/simple_mask}

###### Mask vertebrate repeats using Dfam vertebrate partition (partition 12)
RepeatMasker -pa 16 -a -e ncbi -dir 02_tetrapoda_out -nolow -species vertebrata 01_simple_out/${NEWMASKED}.fasta 2>&1 | tee logs/02_tetrapodamask.log
export ALIGN2=${NEWMASKED}.fasta.align
export CAT2=${NEWMASKED}.fasta.cat.gz
export MASKED2=${NEWMASKED}.fasta.masked
export OUT2=${NEWMASKED}.fasta.out
export TBL2=${NEWMASKED}.fasta.tbl
NEWMASKED2=${MASKED2/simple_mask.masked.fasta/vertebrata_mask}
mv 02_tetrapoda_out/${ALIGN2} 02_tetrapoda_out/${ALIGN2/simple_mask.masked.fasta/vertebrata_mask}
mv 02_tetrapoda_out/${CAT2} 02_tetrapoda_out/${CAT2/simple_mask.masked.fasta/vertebrata_mask}
mv 02_tetrapoda_out/${MASKED2} 02_tetrapoda_out/${NEWMASKED2}.fasta
mv 02_tetrapoda_out/${OUT2} 02_tetrapoda_out/${OUT2/simple_mask.masked.fasta/vertebrata_mask}
mv 02_tetrapoda_out/${TBL2} 02_tetrapoda_out/${TBL2/simple_mask.masked.fasta/vertebrata_mask}

###### Mask vertebrate repeats using known repeats from RepeatModeler
RepeatMasker -pa 16 -a -e ncbi -dir 03_known_out -nolow -lib ${RMFASTA}.known 02_tetrapoda_out/${NEWMASKED2}.fasta 2>&1 | tee logs/03_knownmask.log
export ALIGN3=${NEWMASKED2}.fasta.align
export CAT3=${NEWMASKED2}.fasta.cat.gz
export MASKED3=${NEWMASKED2}.fasta.masked
export OUT3=${NEWMASKED2}.fasta.out
export TBL3=${NEWMASKED2}.fasta.tbl
NEWMASKED3=${${MASKED3}/vertebrata_mask.masked.fasta/known_mask}
mv 03_known_out/${ALIGN3} 03_known_out/${ALIGN3}/vertebrata_mask.masked.fasta/known_mask}
mv 03_known_out/${CAT3} 03_known_out/${CAT3}/vertebrata_mask.masked.fasta/known_mask}
mv 03_known_out/${MASKED3} 03_known_out/${NEWMASKED3}.fasta
mv 03_known_out/${OUT3} 03_known_out/${OUT3}/vertebrata_mask.masked.fasta/known_mask}
mv 03_known_out/${TBL3} 03_known_out/${TBL3}/vertebrata_mask.masked.fasta/known_mask}

###### Mask vertebrate repeats using unknown repeats from RepeatModeler
RepeatMasker -pa 16 -a -e ncbi -dir 04_unknown_out -nolow -lib ${RMFASTA}.unknown 03_known_out/${NEWMASKED3}.fasta 2>&1 | tee logs/03_knownmask.log
export ALIGN4=${NEWMASKED3}.fasta.align
export CAT4=${NEWMASKED3}.fasta.cat.gz
export MASKED4=${NEWMASKED3}.fasta.masked
export OUT4=${NEWMASKED3}.fasta.out
export TBL4=${NEWMASKED3}.fasta.tbl
NEWMASKED4=${MASKED4/known_mask.masked.fasta/uknown_mask}
mv 04_unknown_out/${ALIGN4} 04_unknown_out/${ALIGN4}/known_mask.masked.fasta/uknown_mask}
mv 04_unknown_out/${CAT4} 04_unknown_out/${CAT4}/known_mask.masked.fasta/uknown_mask}
mv 04_unknown_out/${MASKED4} 04_unknown_out/${NEWMASKED4}.fasta
mv 04_unknown_out/${OUT4} 04_unknown_out/${OUT4}/known_mask.masked.fasta/uknown_mask}
mv 04_unknown_out/${TBL4} 04_unknown_out/${TBL4}/known_mask.masked.fasta/uknown_mask}

###### Produce combined files
mkdir -p 05_full_out
export FULLCAT=${GENOME/fasta/full_mask.cat.gz}
# combine full RepeatMasker result files - .cat.gz
cat 01_simple_out/${CAT/fasta/simple_mask} \
02_tetrapoda_out/${CAT2/simple_mask.masked.fasta/vertebrata_mask} \
03_known_out/${CAT3/vertebrata_mask.masked.fasta/known_mask} \
04_unknown_out/${CAT4/known_mask.masked.fasta/uknown_mask} \
> 05_full_out/${FULLCAT}

# combine RepeatMasker tabular files for all repeats - .out
# skip the headers by starting at the 4th line
export FULLOUT=${GENOME/fasta/full_mask.out}
cat 01_simple_out/${OUT/fasta/simple_mask} \
<(cat 02_tetrapoda_out/${OUT2/simple_mask.masked.fasta/vertebrata_mask} | tail -n +4) \
<(cat 03_known_out/${OUT3/vertebrata_mask.masked.fasta/known_mask} | tail -n +4) \
<(cat 04_unknown_out/${OUT4/known_mask.masked.fasta/uknown_mask} | tail -n +4) \
> 05_full_out/${FULLOUT}

# copy RepeatMasker tabular files for simple repeats - .out
cat 01_simple_out/${OUT/fasta/simple_mask} > 05_full_out/${OUT/fasta/simple_mask}

# combine RepeatMasker tabular files for complex, interspersed repeats - .out
export COMPLEXOUT=${GENOME/fasta/complex_mask.out}
cat 02_tetrapoda_out/${OUT2/simple_mask.masked.fasta/vertebrata_mask} \
<(cat 03_known_out/${OUT3/vertebrata_mask.masked.fasta/known_mask} | tail -n +4) \
<(cat 04_unknown_out/${OUT4/known_mask.masked.fasta/uknown_mask} | tail -n +4) \
> 05_full_out/${COMPLEXOUT}

# combine RepeatMasker repeat alignments for all repeats - .align
export FULLALIGN=${GENOME}/fasta/full_mask.align}
cat 01_simple_out/${ALIGN}/fasta/simple_mask} \
02_tetrapoda_out/${ALIGN2}/simple_mask.masked.fasta/vertebrata_mask} \
03_known_out/${ALIGN3}/vertebrata_mask.masked.fasta/known_mask} \
04_unknown_out/${ALIGN4}/known_mask.masked.fasta/uknown_mask} \
> 05_full_out/${FULLALIGN}

# recreate the tabular file using ProcessRepeats, which comes with RepeatMasker
ProcessRepeats -a -species vertebrata 05_full_out/${FULLCAT} 2>&1 | tee logs/05_fullmask.log

# create GFF files using Daren Card's custom script
./rmOutToGFF3custom -o 05_full_out/${FULLOUT} > 05_full_out/${FULLOUT}/out/gff3}
./rmOutToGFF3custom -o 05_full_out/${OUT}/fasta/simple_mask} > 05_full_out/${${OUT}/fasta/simple_mask}/out/gff3}
./rmOutToGFF3custom -o 05_full_out/${COMPLEXOUT} > 05_full_out/${COMPLEXOUT}/out/gff3}

###### Create masked genome files
# create masked genome FASTA files
# create simple repeat soft-masked genome
bedtools maskfasta -soft -fi ${GENOME} -bed 05_full_out/${${OUT}/fasta/simple_mask}/out/gff3} \
-fo 05_full_out/${GENOME}.simple_mask.soft.fasta

# create complex repeat hard-masked genome
bedtools maskfasta -fi 05_full_out/${GENOME}.simple_mask.soft.fasta \
-bed 05_full_out/${COMPLEXOUT}/out/gff3} \
-fo 05_full_out/${GENOME}/fasta/}simple_mask.soft.complex_mask.hard.fasta

