function [ MS_IN ] = centroiding_MS_data( ms_in,min_int,ms_w)
%Copyright (c) 2018 Norwegian Institute for Water Research (NIVA)
%Created by Dr. Saer Samanipour (saer.samanipour@niva.no)
%Please cite: "Combining a Deconvolution and a Universal Library Search
%Algorithm for the Nontarget Analysis of Data-Independent Acquisition Mode 
%Liquid Chromatography-High-Resolution Mass Spectrometry Results", Saer Samanipour,
%Malcolm J. Reid, Kine Bæk, and Kevin V. Thomas, 2018, ES&T, 10.1021/acs.est.8b00259

%
%This function converts the profile MS data to centeroided data.
%

MS_IN=zeros(size(ms_in));
for i=1:size(ms_in,1)
    %tv1=smooth(ms_in(i,:),n_p,'sgolay',order)';
    tv1=ms_in(i,:);
    tv2=find(ms_in(i,:)>=min_int);
    if (isempty(tv2)==0)
        for j=1:size(tv2,2)
            [M1_I]=find(tv1==tv1(tv2(j)),1);
            if (M1_I-ms_w>1&&M1_I+ms_w+1<size(tv1,2))
                tv3=find(tv1(M1_I-ms_w:M1_I)>min_int,1);
                tv4=find(tv1(M1_I:M1_I+ms_w+1)<=min_int,1);
            elseif (M1_I-ms_w<=1&&M1_I+ms_w+1<size(tv1,2))
                tv3=find(tv1(1:M1_I)>min_int,1);
                tv4=find(tv1(M1_I:M1_I+ms_w+1)<=min_int,1);
            elseif (M1_I-ms_w>1&&M1_I+ms_w+1>=size(tv1,2))
                %disp(M1_I)
                tv3=find(tv1(1:M1_I)>min_int,1);
                tv4=find(tv1(M1_I:size(tv1,1))<=min_int,1);
                
            end
            if (isempty(tv3)==0&&isempty(tv4)==0&&tv3<=ms_w+1&&tv4>=1&&M1_I>1)
                if (tv1(M1_I-1)<=tv1(M1_I)&&tv1(M1_I+1)<=tv1(M1_I))
                    MS_IN(i,M1_I)=tv1(M1_I);
                    %disp(j)
                end
                
            elseif (isempty(tv3)==0&&isempty(tv4)==0&&tv3<=ms_w+1&&tv4>=1&&M1_I<1)
                if (tv1(M1_I+1)<=tv1(M1_I))
                    MS_IN(i,M1_I)=tv1(M1_I);
                    
                end
                
                
                
            end
            
            
            
        end
    end
    
       
    
    
end

% plot(ms_in(i,:))
% hold on
% plot(MS_IN(i,:),'--o')
% plot(1000*ones(size(ms_in(i,:))))
% %plot(tv1,'--')
% legend('raw data','center')
% hold off
% pause




end

