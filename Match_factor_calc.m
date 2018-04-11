function [ MatchFactor ] = Match_factor_calc( ref_spec,ref_spec_int,user_spec,user_spec_int,mz_precision,Mass_tollerance )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

%%Mass alignment
total_mz_values=unique(round([ref_spec(ref_spec>0)',user_spec(user_spec>0)],mz_precision));
Ref_int=zeros(1,length(total_mz_values));
User_int=Ref_int;

for i=1:length(total_mz_values)
    
    tv1=total_mz_values(i);
    [M_v,M_i]=min(abs(tv1-ref_spec));
    if (M_v<=Mass_tollerance)
        Ref_int(i)=ref_spec_int(M_i);
    end
    
    
end


for i=1:length(total_mz_values)
    
    tv1=total_mz_values(i);
    [M_v,M_i]=min(abs(tv1-user_spec));
    if (M_v<=Mass_tollerance)
        User_int(i)=user_spec_int(M_i);
    end
    
    
end

%%Normalization

Norm_ref_spec=((total_mz_values.*sqrt(Ref_int))./...
    sum(total_mz_values.*sqrt(Ref_int)));


Norm_User_spec=((total_mz_values.*sqrt(User_int))./...
    sum(total_mz_values.*sqrt(User_int)));

%%Dot-product calculations

Dot_ref=dot(Norm_ref_spec,Norm_ref_spec);

Dot_User=dot(Norm_User_spec,Norm_User_spec);


Dot_direct=dot(Norm_User_spec,Norm_ref_spec);

Dot_revers=dot(Norm_ref_spec,Norm_User_spec);

%%MatchFactor calc
if (Dot_direct<0.01*Dot_User||Dot_revers<0.01*Dot_ref)
    MatchFactor(1)=0;
    MatchFactor(2)=0;
elseif (Dot_direct>0&&Dot_revers>0)
    
    
    Match_Direct=1-abs(Dot_ref-Dot_direct);
    Match_Inverse=1-abs(Dot_User-Dot_revers);
    
    MatchFactor(1)=Match_Direct;
    MatchFactor(2)=Match_Inverse;
else
    MatchFactor(1)=0;
    MatchFactor(2)=0;
end


end

