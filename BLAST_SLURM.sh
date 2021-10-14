#!/bin/bash

#SBATCH --job-name=blaster                              # name of the job submitted
#SBATCH -p brief-low                                        # name of the queue you are submitting to
#SBATCH -N 1                                            # number of nodes in this job
#SBATCH -n 72                                           # number of cores/tasks in this job, you get all 20 cores with 2 threads per core with hyperthreading
#SBATCH -t 2:00:00                                     # time allocated for this job hours:mins:seconds
#SBATCH -o "stdout.%j.%N"                               # standard out %j adds job number to outputfile name and %N adds the node name
#SBATCH -e "stderr.%j.%N"                               # optional but it prints our standard error
#SBATCH --mem=100G   
#SBATCH --mail-type=ALL
#SBATCH --mail-user=julestrachsel@gmail.com

# ENTER COMMANDS HERE:

set -e
module load parallel
module load blast+


collect_blast_results(){

if [ -z ${1+x} ]; then echo "usage: collect_blast_results DB_PREFIX";
                                return 1 ;
else


find . -name "${1}*.blast" | xargs -n 100 cat >> "$1".results
find . -name "${1}*.blast" | xargs rm 
fi
}

export -f collect_blast_results


generate_blast_commands2(){

if [ -z ${1+x} ]; then echo "usage: generate_blast_commands DB_PATH FASTA_PATH";
                                return 1 ;
else
	DB_PATH=$1
        FASTA_PATH=$2
        DB_NAME=$(basename "$DB_PATH")
	fasta_name1=$(basename $FASTA_PATH)
	fasta_name="${fasta_name1%.f*a}"
        output_name=$(printf '%s_%s.blast' $DB_NAME $fasta_name)
	echo "blastn -db $DB_PATH -query $FASTA_PATH -out $output_name -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send slen evalue bitscore qcovs qcovhsp' "
fi


}

export -f generate_blast_commands2


for x in ./blastdbs/*fasta
do
parallel "generate_blast_commands2 ${x%.fasta} {}" ::: $(find ../PROKKA/ -name '*fna')
done > all_blast_commands.txt



# which one?
# cat all_blast_commands.txt | parallel
parallel < all_blast_commands.txt

for x in ./blastdbs/*fasta
	do
	DB=$(basename -s .fasta $x)
	collect_blast_results $DB
done

