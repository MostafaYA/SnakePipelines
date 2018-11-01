#define directories variables
#important not to forget the slash at the end of the name

#directories needed for the denovoassembly.Snakefile
link_dir="results/0_fastq_gz/"
fastq_dir="results/1_fastq/"
fastqc_dir="results/2_fastqc/"
sickle_dir="results/3_sickle/"
SPAdes_dir="results/4_SPAdes/"
FilteredContigs_dir="results/4_FilteredContigs/"
kraken_fastq_dir="mini_kraken/" #at the moment, dont change the name. The scripts need some modifications
coverage_info_dir="results/5_coverageInfo/"
kraken_contig_dir="results/kraken_contig/"
quast_dir="results/6_quast/"
qualimap_dir="results/7_qualimap/"
multiqc_dir="results/8_multiqc/"

#directories needed for the pangenome.Snakefile
prokka_dir="results/5_Prokka/"
prokka_assembly_dir="results/5_Prokka_assemblies/"
gff_dir="results/5_prokka_gff/"
Roary_dir="results/6_Roary/"

#
GENOMES_dir="assemblies"
temporary_todelete="temporary_todelete/"


#FASTQ - done
#fastqc - done
#SICKLE - done
#mini_kraken - done
#mini_kraken_summary
#SPAdes - done
#Coverage_info - doen
#kraken_contig - done
#quast -done
#qualimap
#denovo_assembly_report







#kraken="4_kraken/"
#metaphlan="4_metaphlan/"


refgenome="ref_genome/"
