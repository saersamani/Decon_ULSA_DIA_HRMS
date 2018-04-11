function [ data_c ] = comp_extractor( data,source,mode)
%Copyright (c) 2018 Norwegian Institute for Water Research (NIVA)
%Created by Dr. Saer Samanipour (saer.samanipour@niva.no)
%Please cite: "Combining a Deconvolution and a Universal Library Search
%Algorithm for the Nontarget Analysis of Data-Independent Acquisition Mode 
%Liquid Chromatography-High-Resolution Mass Spectrometry Results", Saer Samanipour,
%Malcolm J. Reid, Kine Bæk, and Kevin V. Thomas, 2018, ES&T, 10.1021/acs.est.8b00259
%
%This fuction performs the filtering of the refernce library for the parent ion.
%
m=size(data.Access,2);

ind=zeros(size(data.Access));




c=1;
for i=1:m
    if (isempty(data.Ion_mode{i})==0&&isempty(data.Ion_source{i})==0)
        if (strcmp(data.Ion_source{i},source)==1&&strcmp(data.Ion_mode{i},mode)==1)
            data_c.mz_values(:,c)=data.mz_values(:,i);
            data_c.mz_Int(:,c)=data.mz_Int(:,i);
            data_c.mz_rel_Int(:,c)=data.mz_rel_Int(:,i);
            data_c.Exact_mass(:,c)=data.Exact_mass(:,i);
            data_c.Name{1,c}=data.Name{i};
            data_c.MS_Type{1,c}=data.MS_Type{i};
            data_c.Access{1,c}=data.Access{i};
            c=c+1;
            
        end
    end
    
end






end

