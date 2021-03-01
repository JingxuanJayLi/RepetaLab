%Function Accumulation

function [tpeak,t, EnerIN,t_percentage]= Accummulation (data)
data1=data(:,1);
data2=data(:,2);
[m,n]=find(data1<0); 
data1(m,:)=[]; 
data2(m,:)=[];
X=data2;
t=data1;
j=1; 
for j=2:length(X)    
    EnerIN(1)=X(1).^2;
    EnerIN(j)=X(j).^2+EnerIN(j-1);
end
EnerIN=EnerIN/(max(EnerIN)); 
[INM05,INM05I]=min(abs(EnerIN-0.05));
[INM95,INM95I]=min(abs(EnerIN-0.95));
delta_tIN=t(INM95I)-t(INM05I);   
halfdelta_t=delta_tIN/2;
tpeak=t(INM05I)+halfdelta_t;
t_total=t(length(t))-t(1);
t_percentage=(delta_tIN/t_total)*100;
end
