# Haplotype analysis and visualization 

The analysis from Anderson T and Sudermann M et al. happens in two phases:

1) Generate haplotype cluster information from VCF file
2) Perform comparisons and visualize the results of the cluster analysis

Code is written in R and Python3. The code was run on a Linux-based system.

This analysis method for the identification of homologous haplotypes is both simple and computationally efficient. 

Sequences were windowed and stepped according to user specifications. Next, a hierarchical agglomerative clustering model is built using Ward distance. Sequences in a window are then clustered together into homologous haplotypes if their pairwise distance falls under a distance threshold. The program determines the distance threshold, *d*, by iterating through a user-specified range and identifying the value of *d* that maximizes the mean silhouette coefficient for all sequences, thus allowing the clustering algorithm to account for unequal information content across genomic windows. The code is built on Scikit-learn’s clustering module.

## PART 1
- Generate chromosome-level VCFs that contain no missing data. This can be accomplished by imputing or by filtering
- Navigate to the folder containing the python3 script "cluster_haplotype.py" and run the script using the parameters of your choosing (see example command below)
- Running the script takes a command of the following usage and a total of 8 REQUIRED arguments
- The following command assues that the repository has been cloned and user is working within HaplotypeAnalysis_Visualization directory
```bash
python3 ./Haplotype_analysis_scripts/cluster_haplotypes.py [vcf_file] [chromosome_basename] [window_size] [window_step_size] [min_snps_cutoff] [min d] [max d] [step d]
```
-  `Argument 1: [vcf_file]` Specify location of uncompressed chromosome level vcf file
-  `Argument 2: [chromosome_basename]` Specify output file basename (typically the chromosome number, i.e. ch09) 
-  `Argument 3: [window_size]`  Set window size to for iterating through the genome - MUST BE AN EVEN NUMBER
-  `Argument 4: [window_step_size]`  Set step size for window iterations
-  `Argument 5: [min_snps_cutoff]`  Set the minimum number of SNPs in each window, otherwise outputs NaN
-  `Argument 6: [min d]`  Set the minimum distance threshold for merging clusters
-  `Argument 7: [max d]`  Set the maximum distance threshold for merging clusters
-  `Argument 8: [step d]`  Set the distance (d) step size

See the associated publication and supplementary methods for information on these parameters

The following python dependencies are required:
> - [`scikit-allel`](https://github.com/cggh/scikit-allel) we used version 1.2.1
> - [`NumPy`](https://numpy.org/)) 
> - [`pandas`](https://pandas.pydata.org/)) 
> - [`scikit-learn`](https://scikit-learn.org/stable/)
> - `system commands (os)`

## EXAMPLE COMMMAND 
The following command can be run in the `HaplotypeAnalysis_Visualization` directory:
```bash
python3 ./Haplotype_analysis_scripts/cluster_haplotypes.py ./Example_files/SL4.0ch09_subset.vcf ch09 250000 100000 10 2 80 10
```
This example code sets the following parameters:  
> - `Argument 1: a chromosome 9 file [will need to change to your input vcf file]`
> - `Argument 2: an output filestem "ch09" [change to chromosome name you want as filestem]`
> - `Argument 3: a window size of 250 Kb [even integer]`  
> - `Argument 4: a window step size of 100 Kb [even integer]`
> - `Argument 5: a minimum number of SNPs of 10 [integar]`
> - `Argument 6: a minimum d of 2 [even integer]`  
> - `Argument 7: a maximum d of 80 [even integer]`  
> - `Argument 8: a step size for d of 10 [integer]`

You can modify the script "parallelize_haplotype_clustering.sh" in order to parallelize by chromosome and by windowing parameters.
Each instance of `cluster_haplotype.py` uses a single thread, so scatter-gather parallelization will save you time in your analysis.

## PART 2
-  Perform interactive visualizations and analyses of the cluster data in RStudio using the [`visualize_haplotypes.Rmd`](Visualization/blob/main/Haplotype_analysis_scripts/visualize_haplotypes.Rmd) script as a starting point.  
-  Inside the Rmarkdown script are useful contrast and plotting functions that will get you started.  
-  Download and view the [`visualize_haplotypes.html`](https://github.com/masudermann/HaplotypeAnalysis_Visualization/blob/main/Example_files/visualize_haplotypes.html) document for completed examples to gain insight into the possibilities of this analysis

Happy introgression hunting!
