# Annotate and Extract Orthologs

Packages required:
* [Tiberius](https://github.com/Gaius-Augustus/Tiberius)
* [OrthoFinder](https://orthofinder.github.io/OrthoFinder/)
* [biopython](https://biopython.org/)
* [macse](https://www.agap-ge2pop.org/macse/)

`tiberius.py`: Runs Tiberius on all specified genomes 

With annotations from all species, next we ran OrthoFinder, first pulling the primary transcript per .faa

```
primary_transcript proteins/
orthofinder -s species_convergence.nwk -f proteins/primary_transcripts/
```

With the .cds sequences for each species, we extracted the nucleotide sequence for each single copy ortholog:
`cds_file.py`: Write the cds sequence with the proper gene id for the ortholog search 
`cds_ortholog.py`: Use the cds files to match the ID names to the orthologs and parse the nucleotide sequence
`edit_headers.py`: Remove anything but the species name from the header

Align nucleotide orthologs with codon-aware macse: 

```
parallel -j 40 --bar --joblog macse_jobs.log \
  '[ -s ../singlecopy_aligned/{/.}_NT.fa ] || java -jar /home/user/Analysis/macse_v2.07.jar -prog alignSequences -seq {} \
     -out_NT ../singlecopy_aligned/{/.}_NT.fa \
     -out_AA ../singlecopy_aligned/{/.}_AA.fa > ../singlecopy_aligned/{/.}.log 2>&1' \
  ::: *.fa
```

