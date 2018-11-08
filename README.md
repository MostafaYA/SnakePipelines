### RUN
#####generate the config file using the script generateConfigSnakemake.sh
'''
./Scripts/generateConfigSnakemake.sh /home/mostafa.abdel/ncbi/public/sra assembly1.config.file Campylobacter
'''
#####edit the snakemake_folder in the config file 
directories:
  snakemake_folder: /home/mostafa.abdel/aProjects/Campylobacter/snakemakeProject/Final-Snake-Project
#####run the pipeline in a dry run 
snakemake -np --snakefile master.Snakefile
