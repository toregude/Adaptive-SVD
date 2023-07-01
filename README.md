# Adaptive-SVD
Reasearch internship at TUM, implementation of algorithms for updating and downdating the SVD of a matrix as used in Hakel factorization in a LTV system. Following files are included, four Matlab files and one report

main.m

Run this file to execute the program. Plots for timing comparison and reconstruction accuracy is shown.

downdate_SVD.m

Function to downdate the SVD based on removal of the top row.

update_SVD.m

Function to update the SVD based on adding of a column to the right. 

combined_SVD.m

Optimized function to first downdate then update the SVD, for speeding purposes.

Adaptive_SVD_for_use_in_Hankel_factorisation_in_a_LTV_system.pdf

Report for this project, including results
