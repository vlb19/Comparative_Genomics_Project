cd NCBI_BLAST

# Update software list and install Python and NCBI BLAST+
sudo apt-get update && sudo apt-get -y install python ncbi-blast+

# Grab data to play with (cow and human RefSeq proteins)
wget ftp://ftp.ncbi.nih.gov/refseq/B_taurus/mRNA_Prot/cow.1.protein.faa.gz
wget ftp://ftp.ncbi.nih.gov/refseq/H_sapiens/mRNA_Prot/human.1.protein.faa.gz

# Unzip database files 
gunzip *gz
ls

# have a quick look at the data format
head cow.1.protein.faa
head human.1.protein.faa

# These are protein sequences in FASTA format (.faa instead of .fasta is NCBI's way of denoting that this a fasta file with amino acids instead of nucleotides)

# Count number of sequences in each file
grep -c '^>' cow.1.protein.faa
grep -c '^>' human.1.protein.faa

# Take a smaller sub-set for practice
head -5 cow.1.protein.faa > cow.small.faa

#############
### BLAST ###
#############

# Now we can blast the two cow sequences against the set of human sequences. 
    # First we need to tell blast about our database
    # BLAST needs to do some pre-work on the database file prior to searching 
    # This helps make the software work a lot faster 
    # make the BLAST database: 
makeblastdb -in human.1.protein.faa -dbtype protein
ls

# this makes a lot of extra files, with the same name as the database and new extensions
# To make BLAST work, these files (aka index files) must be in the same directory as the fasta file 
# Now we can run the blast job using blastp for protein to protein comparisons

blastp -query cow.small.faa -db human.1.protein.faa

# we can save this output to a text file

blastp -query cow.small.faa -db human.1.protein.faa -out cow_vs_human_blast_results.txt
ls

# Take a look at the results using less. There can be more than one match between the query and the same subject. These are referred to as high-scoring segment pairs (HSPs)

less cow_vs_human_blast_results.txt

# Parameters of interest include the -evalue and the -outfmt
# Lets filter for more statistically significant matches with a different output format 
blastp \
-query cow.small.faa \ 
-db human.1.protein.faa \
-out cow_vs_human_blast_results.tab \
-evalue 1e-5 \
-outfmt 7

# Check out results 
less cow_vs_human_blast_results.tab

# Same again with a medium-sized data set 
head -199 cow.1.protein.faa > cow.medium.faa

# Check db size 
grep -c '^>' cow.medium.faa

# Run the blast again but only return the best hit for each query
blastp \
-query cow.medium.faa \
-db human.1.protein.faa \
-out cow_vs_human_blast_results.tab \
-evalue 1e-5 \
-outfmt 6 \
-max_target_seqs 1

# Check out results 
less cow_vs_human_blast_results.tab