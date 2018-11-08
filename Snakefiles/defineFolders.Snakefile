"""----------------------------
Mostafa.Abdel-Glil@fli.de, date: October, 15, 2018, Friedrich-Loeffler-Institut (https://www.fli.de/)
-------------------------------"""
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

#directories needed for the virulenceAndAMR.Snakefile
ariba_amr_dir="results/9_ariba_amr_dir/"
ariba_vf_dir="results/9_ariba_vf_dir/"
abricate_amr_dir="results/9_abricate_amr_dir/"
abricate_vf_dir="results/9_abricate_vf_dir/"
#
GENOMES_dir="assemblies"
temporary_todelete="temporary_todelete/"

#snakemake log and benchmarks
envs_folder="envs/"
benchmarks_folder="results/benchmarks/"
log_folder="results/log/"

#Enviroments
#ariba_env: source activate ariba






#metaphlan="4_metaphlan/"


refgenome="ref_genome/"
