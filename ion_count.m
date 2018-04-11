function [ scores,ions,Mass_error ] = ion_count( ref_spec,user_spec,MassBank,Meas_Mass,ind,Mass_tollerance )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

m_ref=size(ref_spec(ref_spec>0),1);
m_user=size(user_spec(user_spec>0),2);

if (isfield(MassBank,'Exact_mass')==1)
    
    P_mass_e=MassBank.Exact_mass(1,ind)-Meas_Mass;
    
else
    P_mass_e=MassBank-Meas_Mass;
    
end


errors=zeros(1,m_ref);

for i=1:m_ref
    
    errors(i)=min(abs(ref_spec(i)-user_spec));
    
end



Explian_ions=length(errors(errors<=Mass_tollerance));
if (Explian_ions<1)
   scores(1)=-99;
    scores(2)=-99;
    scores(3)=-99;
    scores(4)=-99;
    scores(5)=-99;
    Mass_error=[-99,-99,-99];
elseif (Explian_ions==1)
    ave_error=mean(errors(errors<=Mass_tollerance));
    std_error=std(errors(errors<=Mass_tollerance));
    Mass_error=[P_mass_e,ave_error,-99];
    scores(1)=round(Explian_ions/m_ref,2);
    scores(2)=round(Explian_ions/m_user,2);
    scores(3)=round(Explian_ions/m_ref,2)*round((Mass_tollerance-ave_error)/Mass_tollerance,2);
    scores(4)=0;
    scores(5)=round((Mass_tollerance-abs(P_mass_e))/Mass_tollerance,2);
    
else
    
    ave_error=mean(errors(errors<=Mass_tollerance));
    std_error=std(errors(errors<=Mass_tollerance));
    Mass_error=[P_mass_e,ave_error,std_error];
    scores(1)=round(Explian_ions/m_ref,2);
    scores(2)=round(Explian_ions/m_user,2);
    scores(3)=round(Explian_ions/m_ref,2)*round((Mass_tollerance-ave_error)/Mass_tollerance,2);
    scores(4)=round(Explian_ions/m_ref,2)*round((2*Mass_tollerance-abs(std_error))/(2*Mass_tollerance),2);
    scores(5)=round((Mass_tollerance-abs(P_mass_e))/Mass_tollerance,2);
end



ions{1}=cellstr(strcat(num2str(Explian_ions),'--',num2str(m_ref)));
ions{2}=cellstr(strcat(num2str(Explian_ions),'--', num2str(m_user)));

end

