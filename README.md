# Early blight manuscript (change to proper title of paper) data analysis

The analysis from Anderson T and Sudermann M et al. happens in two phases:

1) Generate haplotype cluster information from vcf file
2) Perform comparisons and visualize the results of the cluster analysis

Note that these notes pertain to an ubuntu linux platform. Code is written in R and python3 and will work on other OS platforms with minor changes

Further details are as follows:

PART 1:
- Generate chromosome-level VCFs that contain no missing data. This can be accomplished by imputing or by filtering
- Navigate to the folder containing the python3 script "cluster_haplotype.py" and run the script using the parameters of your choosing (see example command below)
- Running the script takes a command of the following usage and a total of 8 REQUIRED arguments

"python3 cluster_haplotypes.py [vcf_file] [chromosome_basename] [window_size] [window_step_size] [min_snps_cutoff] [min d] [max d] [step d]"

Argument 1: Specify location of uncompressed chromosome level vcf file
Argument 2: Specify output file basename (typically the chromosome number, i.e. ch09)
Argument 3: Set window size to for iterating through the genome - MUST BE AN EVEN NUMBER
Argument 4: Set step size for window iterations
Argument 5: Set the minimum number of SNPs in each window in order to perform calculations, otherwise outputs NaN
Argument 6: Set the minimum distance threshold for merging clusters
Argument 7: Set the maximum distance threshold for merging clusters
Argument 8: Set the distance (d) step size
See the associated publication and supplementary methods for information on these parameters

The following python dependencies are required:
scikit-allel (we are using version 1.2.1)
numpy
pandas
scikit-learn
system commands (os)

#EXAMPLE COMMMAND - try running the following command in this current directory:
python3 ./cluster_haplotypes.py './example_files/SL4.0ch09_subset.vcf' ch09 250000 100000 10 2 80 10

This example code sets the following parameters:
- a chromosome 9 file [will need to change to your input vcf file]
- an output filestem "ch09" [change to chromosome name you want as filestem]
- a window size of 250 Kb [even integer]
- a window step size of 100 Kb [even integer]
- a minimum of 10 SNPs per window [integer]
- a minimum d of 2 [even integer]
- a maximum d of 80 [even integer]
- a step size for d of 10 [integer]

You can modify the script "parallelize_haplotype_clustering.sh" in order to parallelize by chromosome and by windowing parameters.
Each instance of cluster_haplotype.py uses a single thread, so scatter-gather parallelization will save you time in your analysis.

PART 2:
Perform interactive visualizations and analyses of the cluster data in RStudio using the visualize_haplotypes.Rmd script as a starting point.
Inside the markdown script are useful contrast and plotting functions that will get you started.
See the visualize_haplotypes.html document for completed examples to gain insight into the possibilities of this analysis

If you have any questions, shoot me an email at tayloranthonyanderson@gmail.com
Happy introgression hunting!


