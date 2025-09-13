# Whole Genome and Transcriptome Sequence Analysis of Endophytic Fungi Aspergillus welwitschiae Ocstreb1 (AwOcstreb1)

This document outlines the steps performed for genome assembly, annotation, and downstream analyses of the fungal strain Aspergillus welwitschiae Ocstreb1 genome and transcriptome using various command-line tools. Web-based tools used for some downstream analyses are not included in this documentation.

---

**Fastp:** 

Filtering raw data for poly G and adapter contamination

````python
  fastp -i Forward_reads.fastq.gz -I Reverse_reads.fastq.gz -o fastp_trimmed/r1_fastp_trimmed -O fastp_trimmed/r2_fastp_trimmed --detect_adapter_for_pe --trim_poly_g
  ````
**Abyss:** 

Assembling Reads with Abyss

````python
  abyss-pe name=abyss_assembled k=96 B=50G in='./fastp_trimmed/r1_fastp_trimmed.fq ./fastp_trimmed/r2_fastp_trimmed.fq'
  ````
**Spades:** 

Assembling Reads with Spades using Trusted Contigs from Abyss

````python
  ./SPAdes-3.15.4-Linux/bin/spades.py --isolate -k 55,77,99,121 --trusted-contigs ./abyss_assembled-contigs.fa -1 ./fastp_trimmed/r1_fastp_trimmed.fq -2 ./fastp_trimmed/r2_fastp_trimmed.fq -o ./spades_output
  ````
**Quast:** 

Evaluating assembly quality

````python
  python quast.py spades_output/scaffolds.fasta
  ````
Removing scaffolds<500bp using custom perl script

````python
  perl removesmallcontigs.pl 500 scaffolds.fasta > scaffold_minlen500.fasta
  ````

**BUSCO:** 

Checking Completeness of Assembly 

````python
  busco -m genome -i contigs.fasta -o BUSCO_OUTPUT -l fungi_odb10
  ````

**Fungap:** 

Annoting the Genome

````python
  python $FUNGAP_DIR/fungap.py --output_dir ./fungap_output --trans_read_1 fungi_rnaseq_1.fq.gz --trans_read_2 fungi_rnaseq_2.fq.gz --genome_assembly ./fungap_on_scaffolds_min500/scaffolds_minlen500.fasta --augustus_species aspergillus_fumigatus --busco_dataset ascomycota_odb10 --sister_proteome ./fungap_on_scaffolds_min500/prot_db.faa --num_cores 30
  ````
**BUSCO:** 

Checking completeness of annotation

````python
  busco -m proteins -i fungap_output/fungap_out/fungap_out_prot.faa --lineage_dataset ascomycota_odb10 --out busco_on_fungap_result
  ````
**dREP:** 

Constructing phylogenetic tree based on Average Nucleotide Identity

````python
  dRep dereplicate --S_algorithm ANImf -d --ignoreGenomeQuality output_folder_name -g *
  ````
**Orthofinder:** 

Extracting Primary Transcripts

````python
  for f in *faa ; do python ~/Documents/bioinfo_softwares/OrthoFinder_source/tools/primary_transcript.py $f ; done
  ````
Constructing Phylogenetic Tree Based on Single Copy Ortholog

````python
  ./orthofinder -f /home/intern1/Downloads/7/FUNGI_AN/fungi_all_prot_file_aspergillus_niger_related_strains/drep_result/data/prodigal/primary_transcripts/
  ````
**Diamond:** 

Downloading swissprot database 

````python
  wget https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/swissprot.gz 
  ````
Downloading kog classification file "kog" and kog database kyva 

````python
  wget https://ftp.ncbi.nlm.nih.gov/pub/COG/KOG/kog 
  wget https://ftp.ncbi.nlm.nih.gov/pub/COG/KOG/kyva
````
Creating reference databases

````python
  ./diamond makedb --in kyva -d kog_reference
  ./diamond makedb --in swissprot -d swiss_reference
````
Running Diamond for functional annotation

````python
  ./diamond blastp -d kog_reference -q ./othofinder_proteome/Aspergillus_welwitschiae_ocstreb1.faa --outfmt 6 qseqid sseqid --max-target-seqs 1 --very-sensitive | sort -u > kyva_best_hit.txt
  ./diamond blastp -d swiss_reference -q ./othofinder_proteome/Aspergillus_welwitschiae_ocstreb1.faa --outfmt 6 qseqid sseqid --max-target-seqs 1 --very-sensitive | sort -u > swissprot_best_hit.txt
````
**EggNog-Mapper:**

Extracting transcripts from genome file for eggnog annotation of predicted transcripts:
````python
 gffread fungap_out_Aspwel_AwOcstreb1.gff3 -g AwOcstreb1_scaffolds_minlen500.fasta -w Awocstreb1.genes.fa
````
Running eggnog-mapper for diamond:
````python
 emapper.py -m diamond --data_dir ~/Documents/eggnog_database_with_diamond/ --itype CDS -i high_expressed_genes.fasta -o eggnog_high --cpu 32
````
*Note: For orthologous proteins, eggnog-mapper was run in Protein mode*

**MISA:** 

For Discovering SSR with [MISA](https://webblast.ipk-gatersleben.de/misa/misa_sourcecode_25082020.zip), downloading misa tool

````python
 perl misa.pl /Aspergillus_welwitsciae_strain1.fasta
````
**Command Line BLAST Tool:**

Making custom blast database and running blast for biosynthetic Gene Cluster Identification 

````python
  makeblastdb -in transcriptomic_analysis/cbs.fasta -dbtype nucl
  blastn -query fum_genes/fum_cluster.fasta -db transcriptomic_analysis/ocstreb1.fasta -out fum_cluster_ocstreb1
````
**Hisat2:**

Indexing Genome

````python
 hisat2-build ocstreb1.fasta ocstreb1_genome
````
Aligning Reads to Genome and converting to bam file with **Samtools**

````python
 hisat2 -q -x ocstreb1_genome -1 LCS9631_AwOcstreb1_18_Clean_Data1.fq -2 LCS9631_AwOcstreb1_18_Clean_Data2.fq | samtools sort -o ocstreb1.bam
````
**FeatureCounts:**

Generating count data table 

````python
 featureCounts -p -t exon -g gene_id -Q 10 -a ocstreb1.gff3 -o ocstreb1.fcount.txt ocstreb1.bam
````
