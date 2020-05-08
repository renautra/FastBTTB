# FastBTTB
Fast multiplications for matrices with block Toeplitz Toeplitz block structure
This MATLAB software provides functions that generate the kernels used in gravity and magnetic data forward modeling.

The full matrices are generated for comparison

The transform matrices that are required for fast BTTB implementation are also generated

A script that tests the software  for gravity and magnetic kernels is provided. 

Run Testing_Script.m and make sure all functions and scripts are in the same directory

Run Test_plot to obtain plots of efficiency

Suggestion initial run with lowerscale = 1, upperscale = 2, padding = 1. 

  This will verify that the code is correctly running.
  
  Then run with lowerscale = 1, upperscale = 2, padding = 2. 
  
  This the equivalent run with padding added. 
  
  Finally run (expensive) with lowerscale = 1, upperscale = 12 
  
  (Pick padding = 1 for no padding) 



This is all described in the paper that can be found on arxiv (https://arxiv.org/abs/1912.06976)
A Tutorial and Open Source Software for the Efficient Evaluation of Gravity and Magnetic Kernels (2019)
Jarom Hogue, Rosemary Renaut and Saeed Vatankhah

Scripts: (show how to use the functions)

  Testing_Script.m        : Tests the codes for chosen parameter selections. 

  Test_Efficiency.m       : Calculates timings for use of FFT and without FFT. 

  Test_plot.m             : Produce efficiency and error plots after running Testing_Script.m

  figure_properties.m     : Just makes plots prettier

Functions:

  forward_magnetic.m      : Calculate entries of magnetic kernel matrix
  
  forward_gravity.m       : Calculate entries of gravity kernel matrix
  
  forward_magnetic_bttb.m : Calculate magnetic transform 
  
  forward_gravity_bttb.m  : Calculate gravity transform 
  
  OneMagSliceResponse.m   : Calculate entries of gravity kernel for fixed layer in depth
  
  OneGravSliceResponse.m  : Calculate entries of gravity kernel for fixed layer in depth
  
  matrix_mult_bttb.m      : Matrix multiply for FFT
  
  matsplit.m.             : Helper for identifying parameters
