# Lichen microbiome pipeline

![download_pipeline](https://img.shields.io/github/downloads/atom/atom/total.svg)
![build](https://img.shields.io/appveyor/ci/:user/:repo.svg)
![chat](https://img.shields.io/discord/:serverId.svg)

## Table of Contents

[Binning with metaWRAP](https://github.com/theo-llewellyn/Lichen-metagenomics-pipeline/blob/master/Lichen_microbiome_pipeline.md#Binning-with-metaWRAP)

[Identifying bins with GTDB-tk](https://github.com/theo-llewellyn/Lichen-metagenomics-pipeline/blob/master/Lichen_microbiome_pipeline.md#Identifying-bins-with-GTDB-tk)

[Annotating with DRAM](https://github.com/theo-llewellyn/Lichen-metagenomics-pipeline/blob/master/Lichen_microbiome_pipeline.md#Annotating-with-DRAM)
 * [BLANK](https://github.com/theo-llewellyn/Lichen-metagenomics-pipeline/blob/master/Lichen_microbiome_pipeline.md#)

## Binning with metaWRAP
Once the entire metagenome has been assembled we can separate the bacterial reads into bins which represent populations of individuals of a single bacterial species/OTU. To do this we use the metWRAP pipeline which incorporates various binning tools. The key steps are binning using three different binners (CONCOCT, MaxBin2 and metabat), refining the bins by consolidating the three binners, reassembled the reads using non-metagenome assembler SPAdes, and classifying the bins using the NCBI nt database. We already have the CONCOCT bins so dont need to redownload. We also need to download a local copy of the ncbi nt databse and unzip. This can downloaded from the FTP and should be saved in an ephemeral space as it is a few hundred GB large. metWRAP can then be run as follows:

```
ACCESSION=LIQ74CAUR
mkdir /rds/general/project/theollewellynproject/live/metaWRAP/${ACCESSION}
cd /rds/general/project/theollewellynproject/ephemeral/metaWRAP/${ACCESSION}


#metaWRAP binning \
 -o Initial_bins \
 -t 32 \
 -m 124 \
 -a /rds/general/user/tbl19/home/genomes/${ACCESSION}_megahit_contigs.fasta \
 --metabat2 \
 --maxbin2 \
 /rds/general/project/theollewellynproject/live/data/Trimmed_reads/${ACCESSION}_Filtered_1.fastq \
 /rds/general/project/theollewellynproject/live/data/Trimmed_reads/${ACCESSION}_Filtered_2.fastq

#metaWRAP bin_refinement \
 -o Refined_bins \
 -t 32 \
 -m 124 \
 -c 50 \
 -x 10 \
 --quick \
 -A /rds/general/project/theollewellynproject/live/metaWRAP/${ACCESSION}/Initial_bins/maxbin2_bins/ \
 -B /rds/general/project/theollewellynproject/live/metaWRAP/${ACCESSION}/Initial_bins/metabat2_bins/ \
 -C /rds/general/user/tbl19/home/concoct_output/${ACCESSION}/fasta_bins/

metaWRAP reassemble_bins \
 -o Reassembled_bins \
 -1 /rds/general/project/theollewellynproject/live/data/Trimmed_reads/${ACCESSION}_Trimmed_1.fastq \
 -2 /rds/general/project/theollewellynproject/live/data/Trimmed_reads/${ACCESSION}_Trimmed_2.fastq \
 -t 30 \
 -m 360 \
 -c 50 \
 -x 10 \
 -b Refined_bins/metawrap_50_10_bins

metaWRAP classify_bins \
 -b Reassembled_bins/reassembled_bins \
 -o Classified_bins \
 -t 30
 
rm -r */work_files

```
metaWRAP saves a LOT of working files which need to be deleted after otherwise they take up too much memory. For that reason its also best to save output to an ephemeral space and then move over after work_files have been deleted.

MetaWRAP paper: https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-018-0541-1
MetaWRAP github: https://github.com/bxlab/metaWRAP

You have to provide cut-offs for bin completion and contamination scores. If we want to further refine these after to say extract only the very best bins we can use an if statement like this:

```
# a for loop to filter the reassembled bins to only those of high quality
for i in LIQ*
do
 #changes into reassembled bins directory
 cd ${i}/Reassembled_bins
 #make new directory for high quality bins
 mkdir metawrap_70_5
 # read each line of the bin stats folder
 while read bin completeness contamination GC lineage N50 size
  do
    #need to use the bc command as normal bash if statements cant use decimals
    if echo "$completeness >= 70 && $contamination<= 5" | bc -l | grep -q 1
    then
    #if high quality copy to the new directory
    cp reassembled_bins/${bin}.fa metawrap_70_5
    #also save the stats about that bin in a new stats file
    echo $bin $completeness $contamination $GC $lineage $N50 $size >> reassembled_bins_70_5.stats
    fi
   done < reassembled_bins.stats
 cd ../..
done
```

## Identifying bins with GTDB-tk
Although metaWRAP gives an initial classification of the bins, it relies on NCBI nt database which is not fully curated and the IDs can be a bit shakey at times. Therefore we use the dedicated classifying tools GTDB-tk which incorporates over 200,000 whole bacterial genomes that have been curated and identified. It then extracts marker genes from your bins, aligns them to the database and places them within a phylogenetic context to ID. You need to download the GTDB database, untar and save in a place that can be accessed by the programme. It can then be run as follows
```
GTDBTK_DATA_PATH=/rds/general/user/tbl19/home/software/gtdbtk/release202/

mkdir /rds/general/project/theollewellynproject/ephemeral/gtdbtk/${ACCESSION}
cd /rds/general/project/theollewellynproject/ephemeral/gtdbtk/${ACCESSION}


#identify marker genes
gtdbtk identify \
 --genome_dir /rds/general/project/theollewellynproject/ephemeral/metaWRAP/${ACCESSION}/Reassembled_bins/metawrap_70_5 \
 --out_dir identify \
 --extension fa \
 --cpus 20

#align genomes
gtdbtk align \
 --identify_dir identify \
 --out_dir align \
 --cpus 20

#classify genomes
gtdbtk classify \
 --genome_dir /rds/general/project/theollewellynproject/ephemeral/metaWRAP/${ACCESSION}/Reassembled_bins/metawrap_70_5 \
 --align_dir align \
 --out_dir classify \
 -x fa \
 --cpus 20
```
The final stage is quite RAM heavy so probably need to give at least 120GB of RAM.

## Annotating with DRAM
We will now annotate the genomes using DRAM. This takes as input the reassembled bins from metaWRAP, the CheckM data from metaWRAP and the classification from gtdbtk. We need to edit the headers of the CheckM file just to give the colnames capitals otherwise DRAM wont recognise them. This can be done with a sed loop like so:
```
for i in metaWRAP/*; do sed 's/co/Co/g' ${i}/Reassembled_bins/reassembled_bins.stats >  ${i}/Reassembled_bins/reassembled_bins.stats1; done
```

They can then be run through DRAM with:
```
cd /rds/general/project/theollewellynproject/live/metaWRAP/${ACCESSION}/Reassembled_bins/
cp /rds/general/project/theollewellynproject/live/gtdbtk/${ACCESSION}/classify/gtdbtk.bac120.summary.tsv .

DRAM.py annotate -i 'reassembled_bins/*.fa' \
 -o /rds/general/project/theollewellynproject/ephemeral/DRAM/${ACCESSION}_annotation_1 \
 --min_contig_size 1000 \
 --gtdb_taxonomy gtdbtk.bac120.summary.tsv \
 --checkm_quality reassembled_bins.stats1 \
 --verbose --threads 32

cd /rds/general/project/theollewellynproject/ephemeral/DRAM

DRAM.py distill -i ${ACCESSION}_annotation_1/annotations.tsv -o ${ACCESSION}_genome_summaries_1 --trna_path ${ACCESSION}_annotation_1/trnas.tsv --rrna_path ${ACCESSION}_annotation_1/rrnas.tsv
```
