# A-New-Total-Variational-Algorithm-for-Image-Processing
# Codes are provided by Alireza Hosseini
# hosseini.alireza@ut.ac.ir
# University of Tehran 
# Permission is granted for anyone to copy, use, modify, or distribute this program and accompanying programs and documents for any purpose, provided this copyright notice is retained and prominently displayed, along with a note saying that the original programs are available from the supplementrary material of our paper
#	The programs and documents are distributed without any warranty, express or implied.  As the programs were written for reseach purposes only, they have not been tested to the degree that would be advisable in any important application.  All use of these programs is entirely at the user's own risk.
# How do run the codes?	Two problems are considered; denoising and upscaling. The first; "denoising" which contains Matlab codes for denoisng three test problems with four different total variation models; you can run BARBARARUN.m for denoising test problem "barbara". lenaruuun.m for denoising "lena" test problem and rungold2.m for denoising "goldhill" test image. The new proposed matlab function is main_newalgorithm2.m. The second; "upscaling" which contains Matlab codes for upscaling two test problems with four different total variation models; bikeupsrun.m for upscaling "bike" and fruitsrun.m for upscaling "fruits". The new proposed matlab function is main_new2_ups.m.
# mathlab m-files for "denoising" problem are as follows: 
# "main_isotropic" solves denoising problem for given grayscale image "x" and regularization parameter lambda via isotropic TV.
# "main_upwind" solves denoising problem for given grayscale image "x" and regularization parameter lambda via upwind TV.
# "main_condat" solves denoising problem for given grayscale image "x" and regularization parameter lambda via condat's TV.
# "main_newalgorithm2" solves denoising problem for given grayscale image "x" and regularization parameter lambda via new proposed TV (this code is original).
# "main_isotropic" solves denoising problem for given grayscale image "x" and regularization parameter lambda via isotropic TV.
# "main_lambda" finds optimal regularization parameter "lambda" for any couple of image and TV.
 # mathlab m-files for "upscaling" problem are as follows: 
# "main_isotropic_ups" solves upscaling problem for a given image via isotropic TV.
# "main_upwind_ups" solves upscaling problem for a given image via upwind TV.
# "main_cond_ups" solves upscaling problem for a given image via Condat's TV.
# "main_new2_ups" solves upscaling problem for a given image via new proposed TV (this code is original).



