% main script for running reconstruction 
% adopted for electron ptychography from x-ray ptychograpy
% By Zhen Chen, Cornell University, 2/2020

%%  
clear;

addpath('utils')
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%% SETTINGS of data related parameters  %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% file pointer
fout='sample_pty_'; % strings for saving
dirbase='sample_test1x_crop_run\'; % folder for saving
dir_rawdata='rawdata_21'; % folder for raw data
fn_raw='rawdata_1x_crop.mat'; % file name of raw data

dirbase=fullfile(dir_rawdata,dirbase);
mkdir(dirbase);
fout=fullfile(dirbase,fout);
%% data parameters
exins.fname = fullfile(dir_rawdata,fn_raw);

exins.cutoff = 64;
exins.ADU = 151; % counts to electron, 80 kV
exins.rot_ang = -30 ; % relative rotation between scan and diffraction, exp: -30

exins.df=-500; % focuse
exins.Np_p = [128,128]; % can also pad to 256
exins.scanStepSize_x=0.85; % scan step size in Angstrom
exins.scanStepSize_y=exins.scanStepSize_x; % scan step size 
exins.alpha0 = 21.4;  % convergence angle in mrad
exins.voltage = 80; % electron energy in keV
exins.dataclean = 1;  % 

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%% SETTINGS of reconstruction parameters  %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% parameters for ptychography 
param.Nprobes = 2; % number of probes 
param.probe_reconstruct = 1;  % iteration when the probe reconstruction is started
param.object_reconstruct = 1;% iteration when the object reconstruction is started
param.method = 'MLs';   %% ePIE,ePIE_OPR, DM, MLs (ePIE like partially parallel method), MLc (ML like partially parallel code)
param.Niter = 300;
param.probe_pos_search = 10;            % iteration number from which position correction is started 
param.plot_results_every = 100;  % plot interval

param.beta_object = 1;
param.beta_probe = 1;  
param.likelihood = 'amplitude';  % choose likelihood: amplitude, poisson
param.grouping = 60;  % size of the parallel groups, if this is larger than total, error occurs 

%% GPU or CPU
param.use_gpu = true; % use GPU
param.keep_on_gpu = true; 
param.gpu_id = []; % used GPU id, [] means chosen by matlab

%% Common OPTIONS, usually not change       
param.apply_subpix_shift = true;       % apply FFT-based subpixel shift, important for good position refinement but it is slow

param.variable_probe = true;           % Use SVD to account for variable illumination during a single scan (Orthogonal probe relaxation)
param.variable_probe_modes = 1;         % Number of SVD modes used for the orthogonal probe relaxation, ML method allows only one variable mode
param.variable_intensity = false;      % Account for variable illumination intensity during a single scan

param.beta_LSQ = true;   % use LSQ preditive step for ML method 

param.apply_multimodal_update = false; % apply higher Thibault's modes to the object update, may result in lower object quality for real datasets

param.verbose_level = 0;  % verbosity of the solver , 0 = quiet, 4 = all 
param.delta_p = 0.1;     % LSQ damping constant, improves convergence for the MLc method

% apply farfield support, can be useful constraint in for zoneplate focused beam 
param.probe_support_radius = [];  % radius of circular probe suppport, [] = none, useful for DM code 
param.probe_support_distance = 0; % distance of the probe support from the object plane, 0 = none, inf = farfield

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%% running and saving  %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% prepare data
[inputs,exins] = Data('real', exins);
%% run
param.fout = strcat(fout,'temp_');
tic
[outputs, fourier_error] = core.ptycho_solver(inputs, param);
toc
%% save outputs
inputs.diffraction = [];
save(strcat(fout,'inputs.mat'),'inputs','exins','param','-v7.3');
save(strcat(fout,'refine_outputs.mat'),'fourier_error','outputs','-v7.3');
[O_phase,O_amp]=pty_ramprm_f(outputs,exins,dirbase);
% run pty_ramprm.m;
%%
% Academic License Agreement
%
% Source Code
%
% Introduction 
% •	This license agreement sets forth the terms and conditions under which the PAUL SCHERRER INSTITUT (PSI), CH-5232 Villigen-PSI, Switzerland (hereafter "LICENSOR") 
%   will grant you (hereafter "LICENSEE") a royalty-free, non-exclusive license for academic, non-commercial purposes only (hereafter "LICENSE") to use the cSAXS 
%   ptychography MATLAB package computer software program and associated documentation furnished hereunder (hereafter "PROGRAM").
%
% Terms and Conditions of the LICENSE
% 1.	LICENSOR grants to LICENSEE a royalty-free, non-exclusive license to use the PROGRAM for academic, non-commercial purposes, upon the terms and conditions 
%       hereinafter set out and until termination of this license as set forth below.
% 2.	LICENSEE acknowledges that the PROGRAM is a research tool still in the development stage. The PROGRAM is provided without any related services, improvements 
%       or warranties from LICENSOR and that the LICENSE is entered into in order to enable others to utilize the PROGRAM in their academic activities. It is the 
%       LICENSEE's responsibility to ensure its proper use and the correctness of the results.�?
% 3.	THE PROGRAM IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR 
%       A PARTICULAR PURPOSE AND NONINFRINGEMENT OF ANY PATENTS, COPYRIGHTS, TRADEMARKS OR OTHER RIGHTS. IN NO EVENT SHALL THE LICENSOR, THE AUTHORS OR THE COPYRIGHT 
%       HOLDERS BE LIABLE FOR ANY CLAIM, DIRECT, INDIRECT OR CONSEQUENTIAL DAMAGES OR OTHER LIABILITY ARISING FROM, OUT OF OR IN CONNECTION WITH THE PROGRAM OR THE USE 
%       OF THE PROGRAM OR OTHER DEALINGS IN THE PROGRAM.
% 4.	LICENSEE agrees that it will use the PROGRAM and any modifications, improvements, or derivatives of PROGRAM that LICENSEE may create (collectively, 
%       "IMPROVEMENTS") solely for academic, non-commercial purposes and that any copy of PROGRAM or derivatives thereof shall be distributed only under the same 
%       license as PROGRAM. The terms "academic, non-commercial", as used in this Agreement, mean academic or other scholarly research which (a) is not undertaken for 
%       profit, or (b) is not intended to produce works, services, or data for commercial use, or (c) is neither conducted, nor funded, by a person or an entity engaged 
%       in the commercial use, application or exploitation of works similar to the PROGRAM.
% 5.	LICENSEE agrees that it shall make the following acknowledgement in any publication resulting from the use of the PROGRAM or any translation of the code into 
%       another computing language:
%       "Data processing was carried out using the cSAXS ptychography MATLAB package developed by the Science IT and the coherent X-ray scattering (CXS) groups, Paul 
%       Scherrer Institut, Switzerland."
%
% Additionally, any publication using the package, or any translation of the code into another computing language should cite for difference map:
% P. Thibault, M. Dierolf, A. Menzel, O. Bunk, C. David, F. Pfeiffer, High-resolution scanning X-ray diffraction microscopy, Science 321, 379–382 (2008). 
%   (doi: 10.1126/science.1158573),
% for mixed coherent modes:
% P. Thibault and A. Menzel, Reconstructing state mixtures from diffraction measurements, Nature 494, 68–71 (2013). (doi: 10.1038/nature11806),
% for LSQ-ML method 
% M. Odstrcil, A. Menzel, M. Guizar-Sicairos, Iterative least-squares solver for generalized maximum-likelihood ptychography, Opt. Express, in press (2018). (doi: ).
% for OPRP method 
%  M. Odstrcil, P. Baksh, S. A. Boden, R. Card, J. E. Chad, J. G. Frey, W. S. Brocklesby,  "Ptychographic coherent diffractive imaging with orthogonal probe relaxation." Opt. Express 24, 8360 (2016). (doi: 10.1364/OE.24.008360).
% 6.	Except for the above-mentioned acknowledgment, LICENSEE shall not use the PROGRAM title or the names or logos of LICENSOR, nor any adaptation thereof, nor the 
%       names of any of its employees or laboratories, in any advertising, promotional or sales material without prior written consent obtained from LICENSOR in each case.
% 7.	Ownership of all rights, including copyright in the PROGRAM and in any material associated therewith, shall at all times remain with LICENSOR, and LICENSEE 
%       agrees to preserve same. LICENSEE agrees not to use any portion of the PROGRAM or of any IMPROVEMENTS in any machine-readable form outside the PROGRAM, nor to 
%       make any copies except for its internal use, without prior written consent of LICENSOR. LICENSEE agrees to place the following copyright notice on any such copies: 
%       © All rights reserved. PAUL SCHERRER INSTITUT, Switzerland, Laboratory for Macromolecules and Bioimaging, 2017. 
% 8.	The LICENSE shall not be construed to confer any rights upon LICENSEE by implication or otherwise except as specifically set forth herein.
% 9.	DISCLAIMER: LICENSEE shall be aware that Phase Focus Limited of Sheffield, UK has an international portfolio of patents and pending applications which relate 
%       to ptychography and that the PROGRAM may be capable of being used in circumstances which may fall within the claims of one or more of the Phase Focus patents, 
%       in particular of patent with international application number PCT/GB2005/001464. The LICENSOR explicitly declares not to indemnify the users of the software 
%       in case Phase Focus or any other third party will open a legal action against the LICENSEE due to the use of the program.
% 10.	This Agreement shall be governed by the material laws of Switzerland and any dispute arising out of this Agreement or use of the PROGRAM shall be brought before 
%       the courts of Zürich, Switzerland. 

