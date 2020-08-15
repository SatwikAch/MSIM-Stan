# A monotone single index model for missing-at-random (MAR) longitudinal proportion data
R and stan codes to implement a monotone single index model connected with beta regression on MAR longitudinal proportion data. More details are discussed in the paper.   

Data_gen.R: To generate a simulated longitudinal proportion data set which is missing at random.

fels_beta_reg.stan: Stan code to execute the BR-MSIM which is referred as Model 4 in section 5 of the paper.

Beta_reg_stan_pdf.Rmd: R codes that contain the implementation of stan code on simulated data and plots the single index curve. 

Beta_reg_stan_pdf.pdf: PDF output of the Rmd file. 
 
