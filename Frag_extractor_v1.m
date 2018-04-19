function [ spec ] = Frag_extractor_v1(data_c,Scan_num,W,min_int,Mass_tol,R_min,P_max,r_t,ms_w,Ions,ms_p_w,S2N,time_for_pause)
%Copyright (c) 2018 Norwegian Institute for Water Research (NIVA)
%Created by Dr. Saer Samanipour (saer.samanipour@niva.no)
%Please cite: "Combining a Deconvolution and a Universal Library Search
%Algorithm for the Nontarget Analysis of Data-Independent Acquisition Mode 
%Liquid Chromatography-High-Resolution Mass Spectrometry Results", Saer Samanipour,
%Malcolm J. Reid, Kine Bæk, and Kevin V. Thomas, 2018, ES&T, 10.1021/acs.est.8b00259

%
%This function performs the deconvolution and produces the clean pseudo MS2 spectra of the features.

if (Scan_num-W>1&&Scan_num+W<size(data_c.Mz_values,1))
    ms_v1=data_c.Mz_values(Scan_num-W:Scan_num+W,:);
    ms_in1=data_c.Mz_intensity(Scan_num-W:Scan_num+W,:);
    ms_v2=data_c.Mz_values_2(Scan_num-W:Scan_num+W,:);
    ms_in2=data_c.Mz_intensity_2(Scan_num-W:Scan_num+W,:);
    
elseif (Scan_num-W<=1&&Scan_num+W<size(data_c.Mz_values,1))
    ms_v1=data_c.Mz_values(1:Scan_num+W,:);
    ms_in1=data_c.Mz_intensity(1:Scan_num+W,:);
    ms_v2=data_c.Mz_values_2(1:Scan_num+W,:);
    ms_in2=data_c.Mz_intensity_2(1:Scan_num+W,:);
    
elseif (Scan_num-W>1&&Scan_num+W>=size(data_c.Mz_values,1))
    ms_v1=data_c.Mz_values(Scan_num-W:size(data_c.Mz_values,1),:);
    ms_in1=data_c.Mz_intensity(Scan_num-W:size(data_c.Mz_values,1),:);
    ms_v2=data_c.Mz_values_2(Scan_num-W:size(data_c.Mz_values_2,1),:);
    ms_in2=data_c.Mz_intensity_2(Scan_num-W:size(data_c.Mz_values_2,1),:);
    
    
end



%%
%main ion extraction
% plot(ms_v1(ms_w+1,:),ms_in1(ms_w+1,:))
% pause


ttv1=abs(Ions(1)-ms_v1(ms_w+1,:));
[~,min_l]=min(ttv1);
if (min_l+ms_w>=length(ttv1))
    [~,ttv2]=max(ms_in1(ms_w+1,min_l-ms_w:end));
    Ub=length(ttv1);
    Lb=min_l-ms_w;
elseif (min_l-ms_w<=0)
    [~,ttv2]=max(ms_in1(ms_w+1,1:min_l+ms_w));
    Ub=min_l+ms_w;
    Lb=1;
else
    [~,ttv2]=max(ms_in1(ms_w+1,min_l-ms_w:min_l+ms_w));
    Ub=min_l+ms_w;
    Lb=min_l-ms_w;
end

if (isempty(ttv2)==0)
    ttv3=ms_v1(ms_w+1,Lb:Ub);
    Ions(1)=ttv3(ttv2);
    
    %plot(ms_v1(ms_w+1,Lb:Ub),ms_in1(ms_w+1,Lb:Ub))
    
    
    tv1=abs(Ions(1)-ms_v1);
    
%     figure
%     hold on
    ref_spec=zeros(size(ms_v1,1),1);
    for i=1:size(ms_v1,1)
        [V,J]=min(tv1(i,:));
        if (V<=Mass_tol)
            ref_spec(i)=ms_in1(i,J);
%             plot(ms_v1(i,:),ms_in1(i,:))
%             xlim([103 105])
        end
        clear V J
        
    end
    
    figure
    plot(ref_spec)
    title('XIC of MS1')
    xlabel('Scan number')
    ylabel('Intensity')
    %pause
    pause(time_for_pause)
    close
    %Low energy channel
    Low_E=reshape(Ions,[],2)';
    if (size(Low_E,2)<2)
        MS_v1=Low_E(1,:)';
        MS_In1=Low_E(2,:)';
    else
        MS_v1=Low_E(:,1)';
        MS_In1=Low_E(:,2)';
    end
    
    
    %  plot(ms_v1(W+1,:),ms_in1(W+1,:))
    %  pause
    % hold on
    
    % MS_IN1=zeros(1,size(ms_in1,2));
    % for i=1:size(P_frag,1)
    %     tv1=abs(P_frag(i)-ms_v1);
    %     P_frag_xic=zeros(size(ms_v1,1),1);
    %     for j=1:size(ms_v1,1)
    %         [V,J]=min(tv1(j,:));
    %     if (V<=Mass_tol)
    %         P_frag_xic(j)=ms_in1(j,J);
    %     end
    %     clear V J
    %
    %     end
    %     [M_v,M_i]=max(P_frag_xic);
    %     if (abs(M_i-W+1)<=r_t&&M_v>=min_int)
    %         [R,P] = corrcoef(P_frag_xic,ref_spec);
    %         if (R(2,1)>=R_min&&P(2,1)<=P_max)
    %
    %            MS_IN1(ms_v1(W+1,:)==ms_v1(W+1,abs(ms_v1(W+1,:)-P_frag(i))<=Mass_tol))=ms_in1(W+1,abs(ms_v1(W+1,:)-P_frag(i))<=Mass_tol);
    %
    %
    % %            plot(ms_v1(W+1,:),MS_IN1,'--')
    % %            pause(2)
    %         end
    %
    %
    %     end
    %
    %
    %    clear P_frag_xic
    %
    % end
    
    
    
    
    %%
    %High energy channel
    
    
    Tv2=(find(ms_in2(W+1,:)>=min_int));
    MS_IN2=zeros(1,size(ms_in2,2));
    
    for i=1:size(Tv2,2)
        Tv3=ms_in2(W+1,Tv2);
        Tv4=find(ms_in2(W+1,:)==Tv3(i),1);
        Tv5=abs(ms_v2(W+1,Tv4)-ms_v2);
        Tv6=zeros(size(ms_v2,1),1);
        for j=1:size(ms_v2,1)
            [V,J]=min(Tv5(j,:));
            if (V<=Mass_tol)
                Tv6(j)=ms_in2(j,J);
            end
            clear V J
        end
        
        [M_v,M_i]=max(Tv6);
        if (abs(M_i-W+1)<=r_t&&M_v>=min_int)
            E_p=W+find(ref_spec(W+1:end)<=(1/S2N)*ref_spec(W),1);
            I_p=(W+1)-find(flipud(ref_spec(1:W))<=(1/S2N)*ref_spec(W),1);
            
            
            [R,P] = corrcoef(Tv6(I_p:E_p),ref_spec(I_p:E_p));
            %disp(R)
            if (R(2,1)>=R_min&&P(2,1)<=P_max&&Tv6(W)>=S2N*median(Tv6))
                
                MS_IN2(ms_v2(W+1,:)==ms_v2(W+1,Tv4))=ms_in2(W+1,ms_v2(W+1,:)==ms_v2(W+1,Tv4));
                
                %             plot( ref_spec./max(ref_spec))
                %             hold on
                %             plot(Tv6./max(Tv6),'--')
                %             hold off
                %             %         legend(num2str(ms_v2(W+1,Tv4)),num2str(i))
                %             pause
                %             %
                %             disp(R(1,2))
                %disp(S2N*median(Tv6(I_p:E_p)))
                
                
                
                
            end
            
            
        end
        
        
        
    end
    
%     stem(ms_v2(W+1,:),MS_IN2)
%     hold on
%     stem(ms_v2(W+1,:),ms_in2(W+1,:),'--r')
%     hold off
%     pause
    %%
    %Centeroiding the spectra
    
    %[ ms_int1 ] = centroiding_MS_data( MS_IN1,min_int,ms_w );
    
    [ ms_int2 ] = centroiding_MS_data( MS_IN2,min_int,ms_w);
%     n_p=5;
%     order=3;
    
    
    %[ ms_int2 ] = centroiding_MS_data_WithSmoothing(MS_IN2,min_int,ms_w,n_p,order);
    
    %%
    %mass alignment
    
    for i=1:size(MS_v1,2)
        
        [tv3,tv4]=min(abs(MS_v1(i)-ms_v2(W+1,:)));
        if (ms_int2(tv4)<=0&&tv4>W&&tv4+W<=size(ms_int2,2))
            tv5=find(ms_int2(tv4-W:tv4+W)>=min_int);
            for k=1:size(tv5,2)
                tv6=tv4-W+tv5(k)-1;
                if (isempty(tv6)==0&&abs(tv6-tv4)<=W&&abs(ms_v2(W+1,tv4)-ms_v2(W+1,tv6))<=ms_p_w)
                    ms_int2(tv4)=ms_int2(tv6);
                    ms_int2(tv6)=0;
                    
                end
                
            end
            
        elseif (ms_int2(tv4)<=0&&tv4>W&&tv4+W>size(ms_int2,2))
            tv5=find(ms_int2(tv4-W:end)>=min_int);
            for k=1:size(tv5,2)
                tv6=tv4-W+tv5(k)-1;
                if (isempty(tv6)==0&&abs(tv6-tv4)<=W&&abs(ms_v2(W+1,tv4)-ms_v2(W+1,tv6))<=ms_p_w)
                    ms_int2(tv4)=ms_int2(tv6);
                    ms_int2(tv6)=0;
                    
                end
                
            end
            
            
        end
        
        clear tv2 tv3 tv4 tv5 tv6
        
    end
    
    
    
    
    %%
    
    
    
    spec.ms_v1=MS_v1;
    spec.ms_v2=ms_v2(W+1,ms_int2>0);
    spec.ms_in1=MS_In1;
    spec.ms_in2=ms_int2(ms_int2>0);
    
    
    %disp(spec)
    
else 
    spec.ms_v1=[];
    spec.ms_v2=[];
    spec.ms_in1=[];
    spec.ms_in2=[];
    disp('This target mass does not exist in this sample')
    
    
end


end

