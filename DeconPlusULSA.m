function DeconPlusULSA(Project_name,File_name,Path_Chrom,...
    Path_Peak_List,Path_Library,path_adducts,mode,source,W,...
    min_int,Parent,ms_tol_w,Mass_tollerance,...
    R_min,P_max,r_t,S2N,ms_p_w,ms_w,mz_precision,Weigth_fun,time_for_pause)
%Copyright (c) 2018 Norwegian Institute for Water Research (NIVA)
%Created by Dr. Saer Samanipour (saer.samanipour@niva.no)
%Please cite: "Combining a Deconvolution and a Universal Library Search
%Algorithm for the Nontarget Analysis of Data-Independent Acquisition Mode 
%Liquid Chromatography?High-Resolution Mass Spectrometry Results", Saer Samanipour,
%Malcolm J. Reid, Kine Bæk, and Kevin V. Thomas, 2018, ES&T, 10.1021/acs.est.8b00259

%This is the function for the complete workflow. The steps in the workflow
%are separated in order to increase readability.  
%%
%Create the analysis folder 
p=pwd();
Foldername=strcat(p,'\',Project_name);

if (exist(Foldername,'dir')==0)
    mkdir(Foldername)
    cd(Foldername)
    mkdir(strcat(File_name,'_Spec'))
    mkdir(strcat(File_name,'_IDs'))
else 
    cd(Foldername)
    
end

%%
%Import the Chromatogram 

load(Path_Chrom)
%%
%Import the Peak List
List=readtable(Path_Peak_List);

%%
%Import the Reference Library

load(Path_Library,'data');

[ MassBank ] = comp_extractor( data,source,mode);

%%
%Deconvolution and ULSA

for i=1:size(List.ID,1)
    if (isnan(List.ID(i))==0)
        Scan_num=round(List.Retention_time(i));
        Ions=[round(List.MZ(i),3),100];
        
        [ spec ] = Frag_extractor_v1(data_c,Scan_num,W,min_int,Mass_tollerance,R_min,P_max,r_t,ms_w,Ions,ms_p_w,S2N,time_for_pause);
        if (isempty(spec.ms_v2)==0)
            
            Tab=table(spec.ms_v2',spec.ms_in2','VariableNames',{'Mass','Intensity'});
            writetable(Tab,strcat(p,'\',Project_name,'\',strcat(File_name,'_spec'),'\Target_',num2str(List.ID(i)),'.txt'))
        end
        
    end
    
    %%
    %Library search and save the results
    if (isempty(spec.ms_v2)==0)
        
        [ Final_table ] = Lib_search(Parent,ms_tol_w,MassBank,spec,mz_precision,Mass_tollerance,Weigth_fun,mode,path_adducts);
        disp(Final_table)
        if (isempty(Final_table)==0)
            writetable(Final_table,strcat(p,'\',Project_name,'\',strcat(File_name,'_IDs'),'\Target_',num2str(List.ID(i)),'.txt'))
            
            
        end
        
    end
    pause(time_for_pause)
end


end

