%Copyright (c) 2018 Norwegian Institute for Water Research (NIVA)
%Created by Dr. Saer Samanipour (saer.samanipour@niva.no)
%Please cite: "Combining a Deconvolution and a Universal Library Search
%Algorithm for the Nontarget Analysis of Data-Independent Acquisition Mode 
%Liquid Chromatography?High-Resolution Mass Spectrometry Results", Saer Samanipour,
%Malcolm J. Reid, Kine Bæk, and Kevin V. Thomas, 2018, ES&T, 10.1021/acs.est.8b00259
%
%Before starting your analysis, please make sure that your data is
%configured properly, see the user manual. 
%Please assure to read the licensing agreement called License.txt. 
close ,clear ,clc

cd('C:\GitHub\Library_search_fun')                                          %The main folder where all the fun are stored
Project_name='NMR_test_5';                                                  %This is the project name (i.e. the folder for storing the results)
File_name='DENMARK';                                                        %This is the the name of the inner folder
Path_Chrom='C:\GitHub\Library_search_fun\Chrom.mat';                        %Path to the chromatogram data saved as a "mat" file (please see the user manual for the structure of the chromatogram)
Path_Peak_List='C:\GitHub\Library_search_fun\Targets_1.xlsx';               %Path to the list of features to be identified. An example of that file is provided.
Path_Library='C:\GitHub\Library_search_fun\MassBank_matlab.mat';            %Path to the reference library saved as a "mat" file (please see the user manual for the structure of that data)
path_adducts='C:\GitHub\Library_search_fun\Pos_adducts.xlsx';               %Path to the adducts and isotops file. This file is stored in the main folder
source='ESI';                                                               %This parameter sets the ion source.
mode='POSITIVE';                                                            %This parameter sets the polarity.
W=15;                                                                       %This is the retention window for your analysis. This parameter should be set to around two times the peak width for the widest peak.
min_int=300;                                                                %This is the minimum intensity for an MS2 peak to be considered.                                                                %This 
R_min=0.80;                                                                 %This is to define the minimum correlation coefficient for an XIC2 to be considered.
P_max=0.05;                                                                 %This is the maximum p value for the corrlation coefficent
r_t=3;                                                                      %This is the retention tollerance of the XIC2
ms_w=15;                                                                    %Maximum number of points in an MS peak
ms_p_w=0.05;                                                                %Maximum accepatble MS peak width
S2N=3;                                                                      %Minimum signal to noise ratio for XIC2 (noise is calculated using the median of the retention window).
Parent=1;                                                                   %Defines if the analyzed feature is the parent ion or unknow (1 for parent and 0 for unknown).                                                                   
ms_tol_w=0.01;                                                              %The mass window to filter the candidates in the library.
Mass_tollerance=0.01;                                                       %The mass tollerance for accepting the fragments.
mz_precision=3;                                                             %The number of digits in the MS data when doing the processing
Weigth_fun=[0.8,0.8,1,1,1,1,1];                                             % The weight values [ion_match1,ion_match2,MZ_e1,MZ_e2,MZ_e3,D_match,R_match]
time_for_pause=1;                                                           %Pause time for the figure to be displayed. Set to zero for batch analysis

DeconPlusULSA(Project_name,File_name,Path_Chrom,...
    Path_Peak_List,Path_Library,path_adducts,mode,source,W,...
    min_int,Parent,ms_tol_w,Mass_tollerance,...
    R_min,P_max,r_t,S2N,ms_p_w,ms_w,mz_precision,Weigth_fun,time_for_pause)

