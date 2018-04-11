function [ Final_table ] = Lib_search(Parent,ms_tol_w,MassBank,spec,mz_precision,Mass_tollerance,Weigth_fun,mode,path)
%Copyright (c) 2018 Norwegian Institute for Water Research (NIVA)
%Created by Dr. Saer Samanipour (saer.samanipour@niva.no)
%Please cite: "Combining a Deconvolution and a Universal Library Search
%Algorithm for the Nontarget Analysis of Data-Independent Acquisition Mode 
%Liquid Chromatography-High-Resolution Mass Spectrometry Results", Saer Samanipour,
%Malcolm J. Reid, Kine Bæk, and Kevin V. Thomas, 2018, ES&T, 10.1021/acs.est.8b00259
%
%This function performs the library sereach and produces the final list of candidates. 
%
Meas_Mass=spec.ms_v1;
if (Parent==1)
    
    ind_mz=find(abs(Meas_Mass-MassBank.Exact_mass(1,:))<=ms_tol_w);
    if (isempty(ind_mz)==0)
        
        for i=1:length(ind_mz)
            ref_spec=MassBank.mz_values(:,ind_mz(i));
            ref_spec_int=MassBank.mz_Int(:,ind_mz(i));
            user_spec=spec.ms_v2;
            user_spec_int=spec.ms_in2;
            [ S1,ions,Mass_error ] = ion_count( ref_spec,user_spec,MassBank,Meas_Mass,ind_mz(i),Mass_tollerance );
            [ MatchFactor ] = Match_factor_calc( ref_spec,ref_spec_int,user_spec,user_spec_int,mz_precision,Mass_tollerance );
            A={'nr','Name','Ref_ion_matched','User_ion_matched','Parent_ion_error','Average_frag_error',...
                'Std_frag_error','Direct_match','Reverse_match','Score'};
            final_results(i,:)=table(i,MassBank.Name{ind_mz(i)},ions{1},ions{2},Mass_error(1),Mass_error(2),Mass_error(3),MatchFactor(1),MatchFactor(2),...
                sum(Weigth_fun.*[S1,MatchFactor]),'VariableNames',A);
            
        end
        
        Final_results=sortrows(final_results,10,'descend');
        Final_results.nr=linspace(1,length(Final_results.nr),length(Final_results.nr))';
        m=find(Final_results.Score<=0,1);
        Final_table=Final_results(1:m-1,:);
    else 
        A={'nr','Name','Ref_ion_matched','User_ion_matched','Parent_ion_error','Average_frag_error',...
                'Std_frag_error','Direct_match','Reverse_match','Score'};
            final_results=cell2table({1,'-','-','-',-99,-99,-99,-99,-99,...
                -99});
            final_results.Properties.VariableNames=A;
        Final_table=final_results;    
    end
    
elseif (Parent==0)
    
    
    if (strcmp(mode,'POSITIVE')==1)
        adducts=readtable(path);
        Adducts=unique(adducts.Mass);
        res=cell(size(Adducts,1),1);
        for i=1:size(Adducts,1)
            tv1=str2num(cell2mat(Adducts(i)));
            
            ind_mz=find(abs((Meas_Mass-tv1+1.007276)-MassBank.Exact_mass(1,:))<=ms_tol_w);
            A={'nr','Name','Ref_ion_matched','User_ion_matched','Parent_ion_error','Average_frag_error',...
                        'Std_frag_error','Direct_match','Reverse_match','Score'};
            if (isempty(ind_mz)==0)
                for j=1:length(ind_mz)
                    %disp(i)
                    ref_spec=MassBank.mz_values(:,ind_mz(j));
                    ref_spec_int=MassBank.mz_Int(:,ind_mz(j));
                    user_spec=spec.ms_v2;
                    user_spec_int=spec.ms_in2;
                    [ S1,ions,Mass_error ] = ion_count( ref_spec,user_spec,MassBank,Meas_Mass,ind_mz(j),Mass_tollerance );
                  
                    if (isempty(S1(S1>0))==0)
                        [ MatchFactor ] = Match_factor_calc( ref_spec,ref_spec_int,user_spec,user_spec_int,mz_precision,Mass_tollerance );
                        
                        final_results(j,:)=table(j,MassBank.Name{ind_mz(j)},ions{1},ions{2},Mass_error(1),Mass_error(2),Mass_error(3),MatchFactor(1),MatchFactor(2),...
                            sum(Weigth_fun.*[S1,MatchFactor]),'VariableNames',A);
                    else 
                        final_results(j,:)=table(j,MassBank.Name{ind_mz(j)},ions{1},ions{2},-99,-99,-99,-99,-99,...
                            0,'VariableNames',A);
                        
                    end
                    
                end
                
                Final_results=sortrows(final_results,10,'descend');
                Final_results.nr=linspace(1,length(Final_results.nr),length(Final_results.nr))';
                res{i}=Final_results;
                clear Final_results final_results
            else
                res{i}=[];
            end
            
            
        end
        
        cellsz=zeros(1,size(res,1));
        IND=cellsz;
        for i=1:size(res,1)
            tv1=res{i};
            if (isempty(tv1)==0)
                tv2=tv1(tv1.Score>0,:);
                if (isempty(tv2)==0)
                    cellsz(i)=size(tv2,1);
                    IND(i)=i;
                end
            end
            
        end
        
        Cellsz=sum(cellsz);
        num_nz=length(cellsz(cellsz>0));
        IND_1=IND(IND>0);
        Resu=table(1,{'a'},{'a'},{'a'},1,1,1,1,1,...
            -99,'VariableNames',A);
        
        for i=1:length(IND_1)
            tv1=res{IND_1(i)};
            Resu=vertcat(Resu,tv1);
            
            
            
        end
        
        Final_table=sortrows(Resu,10,'descend');
        Final_table=Final_table(Final_table.Score>0,:);
        Final_table.nr=linspace(1,length(Final_table.nr),length(Final_table.nr))';
        
        
        
        
    elseif (strcmp(mode,'NEGATIVE')==1)
        adducts=readtable(path);
        Adducts=unique(adducts.Mass);
        res=cell(size(Adducts,1),1);
        for i=1:size(Adducts,1)
            tv1=str2num(cell2mat(Adducts(i)));
            
            ind_mz=find(abs((Meas_Mass+tv1-1.007276)-MassBank.Exact_mass(1,:))<=ms_tol_w);
            A={'nr','Name','Ref_ion_matched','User_ion_matched','Parent_ion_error','Average_frag_error',...
                        'Std_frag_error','Direct_match','Reverse_match','Score'};
            if (isempty(ind_mz)==0)
                for j=1:length(ind_mz)
                    %disp(i)
                    ref_spec=MassBank.mz_values(:,ind_mz(j));
                    ref_spec_int=MassBank.mz_Int(:,ind_mz(j));
                    user_spec=spec.ms_v2;
                    user_spec_int=spec.ms_in2;
                    [ S1,ions,Mass_error ] = ion_count( ref_spec,user_spec,MassBank,Meas_Mass,ind_mz(j),Mass_tollerance );
                  
                    if (isempty(S1(S1>0))==0)
                        [ MatchFactor ] = Match_factor_calc( ref_spec,ref_spec_int,user_spec,user_spec_int,mz_precision,Mass_tollerance );
                        
                        final_results(j,:)=table(j,MassBank.Name{ind_mz(j)},ions{1},ions{2},Mass_error(1),Mass_error(2),Mass_error(3),MatchFactor(1),MatchFactor(2),...
                            sum(Weigth_fun.*[S1,MatchFactor]),'VariableNames',A);
                    else 
                        final_results(j,:)=table(j,MassBank.Name{ind_mz(j)},ions{1},ions{2},-99,-99,-99,-99,-99,...
                            0,'VariableNames',A);
                        
                    end
                    
                end
                
                Final_results=sortrows(final_results,10,'descend');
                Final_results.nr=linspace(1,length(Final_results.nr),length(Final_results.nr))';
                res{i}=Final_results;
                clear Final_results final_results
            else
                res{i}=[];
            end
            
            
        end
        
        cellsz=zeros(1,size(res,1));
        IND=cellsz;
        for i=1:size(res,1)
            tv1=res{i};
            if (isempty(tv1)==0)
                tv2=tv1(tv1.Score>0,:);
                if (isempty(tv2)==0)
                    cellsz(i)=size(tv2,1);
                    IND(i)=i;
                end
            end
            
        end
        
        Cellsz=sum(cellsz);
        num_nz=length(cellsz(cellsz>0));
        IND_1=IND(IND>0);
        Resu=table(1,{'a'},{'a'},{'a'},1,1,1,1,1,...
            -99,'VariableNames',A);
        
        for i=1:length(IND_1)
            tv1=res{IND_1(i)};
            Resu=vertcat(Resu,tv1);
            
            
            
        end
        
        Final_table=sortrows(Resu,10,'descend');
        Final_table=Final_table(Final_table.Score>0,:);
        Final_table.nr=linspace(1,length(Final_table.nr),length(Final_table.nr))';
        
    end
    
    
    
end


end

