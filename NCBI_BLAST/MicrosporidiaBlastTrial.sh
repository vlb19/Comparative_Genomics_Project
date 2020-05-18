# copy data into new directory
cp ../data/Amphiamblys_sp._WSBS2006/MicrosporidiaDB-46_AspWSBS2006_AnnotatedProteins.fasta Amphiamblys_WSBS2006.faa 
cp ../data/Saccharomyces_cerevisiae_S288C/GCF_000146045.2_R64_protein.faa Saccharomyces_cerevisiae_S288C.faa

#head the microsporidia data
head ../data/Amphiamblys_sp._WSBS2006/MicrosporidiaDB-46_AspWSBS2006_AnnotatedProteins.fasta
head ../data/Anncaliia_algerae_PRA109/MicrosporidiaDB-46_AalgeraePRA109_AnnotatedProteins.fasta
head ../data/Anncaliia_algerae_PRA339/MicrosporidiaDB-46_AalgeraePRA339_AnnotatedProteins.fasta
head ../data/Edhazardia_aedis_USNM_41457/MicrosporidiaDB-46_EaedisUSNM41457_AnnotatedProteins.fasta

# head the Saccharomyces cerevisiae data
    # Unzip database files 
    gunzip ../data/Saccharomyces_cerevisiae_S288C/*gz
head ../data/Saccharomyces_cerevisiae_S288C/GCF_000146045.2_R64_protein.faa

# Count number of sequences in each file
grep -c '^>' ../data/Amphiamblys_sp._WSBS2006/MicrosporidiaDB-46_AspWSBS2006_AnnotatedProteins.fasta

grep -c '^>' ../data/Edhazardia_aedis_USNM_41457/MicrosporidiaDB-46_EaedisUSNM41457_AnnotatedProteins.fasta

grep -c '^>' ../data/Saccharomyces_cerevisiae_S288C/GCF_000146045.2_R64_protein.faa

# Take a smaller sub-set for practice
head -199 ../data/Amphiamblys_sp._WSBS2006/MicrosporidiaDB-46_AspWSBS2006_AnnotatedProteins.fasta > Ampiamblys_WSBS2006.small.faa
head -199 ../data/Edhazardia_aedis_USNM_41457/MicrosporidiaDB-46_EaedisUSNM41457_AnnotatedProteins.fasta > EaedisUSNM41457.small.faa

# Make blast database from S. cerevisiae
makeblastdb -in Saccharomyces_cerevisiae_S288C.faa -dbtype prot

# Now we can run the blast job using blastp for protein to protein comparisons
blastp -query EaedisUSNM41457.small.faa -db Saccharomyces_cerevisiae_S288C.faa -out Eaedis_vs_Yeast_testblast_results.txt
ls

# Take a look at the results using less. There can be more than one match between the query and the same subject. These are referred to as high-scoring segment pairs (HSPs)

less Eaedis_vs_Yeast_testblast_results.txt

# Parameters of interest include the -evalue and the -outfmt
# Lets filter for more statistically significant matches with a different output format 
blastp \
-query ../data/Edhazardia_aedis_USNM_41457/MicrosporidiaDB-46_EaedisUSNM41457_AnnotatedProteins.fasta \ 
-db Saccharomyces_cerevisiae_S288C.faa \
-out Eaedis_vs_Yeast_testblast_results.tab \
-evalue 1e-5 \
-outfmt 7

# Check out results 
less Eaedis_vs_Yeast_testblast_results.tab

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