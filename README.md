This analysis was performed using Linux Command Line Tools
Fastp: 
fastp -i DU_PBTlab_S6_L001_R1_001.fastq.gz -I DU_PBTlab_S6_L001_R2_001.fastq.gz -o fastp_trimmed/r1_fastp_trimmed -O fastp_trimmed/r2_fastp_trimmed --detect_adapter_for_pe --trim_poly_g

Abyss:
abyss-pe name=abyss_assembled k=96 B=50G in='./fastp_trimmed/r1_fastp_trimmed.fq ./fastp_trimmed/r2_fastp_trimmed.fq'

Spades: 

Quast:
On spades scaffolds: python /home/cblast/miniconda3/envs/wgs_analysis/bin/quast.py /home/cblast/Documents/7_f\ proj/spades_output/scaffolds.fasta

Removed scaffolds<500bp with 
perl removesmallcontigs.pl 500 scaffolds.fasta > scaffold_minlen500.fasta


BUSCO (completeness of assembly): 
busco -m genome -i contigs.fasta -o BUSCO_OUTPUT -l fungi_odb10


Fungap:
python $FUNGAP_DIR/fungap.py --output_dir ~/Documents/7_f_proj/fungap_finalrun_on_spades_plus_abyss_scaffolds_min500/fungap_with_rnaseq19 --trans_read_1 ~/Documents/7_f_proj/fungap_finalrun_on_spades_plus_abyss_scaffolds_min500/fungi_19_rnaseq_1.fq.gz --trans_read_2 ~/Documents/7_f_proj/fungap_finalrun_on_spades_plus_abyss_scaffolds_min500/fungi_19_rnaseq_2.fq.gz --genome_assembly ~/Documents/7_f_proj/fungap_finalrun_on_spades_plus_abyss_scaffolds_min500/scaffolds_minlen500.fasta --augustus_species aspergillus_fumigatus --busco_dataset ascomycota_odb10 --sister_proteome ~/Documents/7_f_proj/fungap_finalrun_on_spades_plus_abyss_scaffolds_min500/prot_db.faa --num_cores 30


BUSCO (completeness of annotation):
busco -m proteins -i fungap_with_rnaseq19/fungap_out/fungap_out_prot.faa --lineage_dataset ascomycota_odb10 --out busco_on_fungap_result
dREP:
dRep dereplicate --S_algorithm ANImf -d --ignoreGenomeQuality output_folder_name -g * (on all genome fasta file)
Orthofinder:
for f in *faa ; do python ~/Documents/bioinfo_softwares/OrthoFinder_source/tools/primary_transcript.py $f ; done
Diamond: 
Downloading swissprot database- 
wget https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/swissprot.gz 
Downloading kog database kyva and kog classification kog file- 
wget https://ftp.ncbi.nlm.nih.gov/pub/COG/KOG/kog 
wget https://ftp.ncbi.nlm.nih.gov/pub/COG/KOG/kyva
Creating reference databases
./diamond makedb --in kyva -d kog_reference
./diamond makedb --in swissprot -d swiss_reference
Running Diamond
./diamond blastp -d kog_reference -q ~/Documents/7_f_proj/proteomes_for_orthofinder_run/Aspergillus_welwitschiae_ocstreb1.faa --outfmt 6 qseqid sseqid --max-target-seqs 1 --very-sensitive | sort -u > kyva_best_hit.txt
MISA: 
Download misa tool- 
wget https://webblast.ipk-gatersleben.de/misa/misa_sourcecode_25082020.zip
perl misa.pl /Aspergillus_welwitsciae_strain1.fasta
