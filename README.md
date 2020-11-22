Codes and supplementary materials for the coming manuscript on mutation accumulation in DNA repair deficient *C. elegans*.

B. Meier, N. Volkova, Y Hong, S Bertolini, B Wang, V Gonzalez-Huici, S Boulton, PJ Campbell, M Gerstung and A Gartner

"Protection of the *C. elegans* germ cell genome depends on diverse DNA repair pathways during normal proliferation"
 
The manuscript is currently under preparation. For any questions or suggestions, please contact me via volk1189@gmail.com.

Script structure:

`Location_plots.Rmd` - contains the plotting of mutation distributions across the genome in multiple samples and clustering analysis. Requires samples annotation table, lists of VCFs for base substitutions and indels, as well as a table of structural variants.

`Mutation_rate_estimation_final.R`	- codes for estimation of mutation rates and mutational signatures, analysis of differences to wild-type, and figure plotting. Requires a matrix with mutation counts and sample annotation table.

`SV_info.R`	- in-depth analysis of structural variants

`dog_1_analysis.R`	- separate script producing plot with distributions of G4 and non-G4 associated mutations in dog-1 mutants.

`plotting_complex_SVs.R` - separate script for plotting complex structural variants, requires a table of SVs as well as .bw and bamstat files for each sample (including the reference samples).

`plotting_functions.R`	- supplementary fuctions for signature plotting.

`useful_functions.R`	- supplementary functions for analysis.
