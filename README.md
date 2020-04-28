Addtional notes about the codes on electron ptychography
2020/02/22 Dr. Zhen Chen, Prof. David A. Muller Group, Cornell University
	 
----------------------
----------------------
To run reconstruction, run the main drive script: ptychography_clean_drive.m. Before doing this, please refer to the notes below:
 
1. read README.md for details about the algorithms and settings. 

2. Experimental data (supple. Fig. 4 in the manuscript) is assumed to be put in the same folder as the main folder of the codes with the path: 
    .\rawdata_21\rawdata_1x_crop.mat
	You need to change the drive script accordingly if you put the data somewhere else.
		
3. Look into the script ptychography_clean_drive.m and check / modify parameters.

4. Important parameters are parameters for data related and some of the reconstruction ones.

5. For the example dataset and using the default parameters (300 iterations with 128 x 128 pixels for each diffraction), 
   it will take less than 10 minutes on a decent GPU card (memory > 1 GB) and about two hours on CPU.
   It is the top-left corner and only about contains 1/4 diffractoins. 
   
   if no computational GPU card is available, please changes to use CPU:
   param.use_gpu = false; % use GPU
   
6. The outputs are in a folder defined by variable 'dir_base' within data folder and the phase image is under name:
     MLs_backgroundremove_final_crop_phase.png
	 All input parameters are stored in sample_pty_inputs.mat 
	 Clean reconstructions are in MLs_backgroundremove_final_crop_data.mat
	 raw reconstructions are in sample_pty_refine_outputs.mat
	 
----------------------
----------------------

License agreements
----------------------

The main ptychography toolkit developed at Paul Scherrer Institut, Switzerland is available at
 [www.psi.ch/sls/csaxs/software](http://www.psi.ch/sls/csaxs/software)
Copyright and license issues should follow the agreements in their codes and/or refer to their website. 

-----------------------
The interface and data handling parts of the codes were writen by Zhen Chen. 
Details on the data set and collection conditions can be found in
Zhen Chen, Michal Odstrcil, Yi Jiang, Yimo Han, Ming-Hui Chiu, Lain-Jong Li, David A. Muller
Mixed-state electron ptychography enables sub-angstrom resolution imaging with picometer precision at low dose 
under review on Nature Communications.

This paper should be cited whenever this code, dataset or their derivatives are used.
 
----------------------
----------------------

Contact information
----------------------
Zhen Chen (zhen.chen@conell.edu) or David A. Muller (david.a.muller@cornell.edu)
