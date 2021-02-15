%Function for CMA-C, which relies on MultiMSAlignGa to align the spectra, and rawEIC
clear; close all; clc
  load('_orbtimeStdMix10nM.mat');
load('_ms1spectraStdMix10nM.mat');


%here we specify the time window
tfake=find(orbtime<2195 | orbtime>2395);
orbtime(tfake)=[];
ms1spectra(tfake)=[];


tic
for i=1:length(orbtime)

for j= 1:length(ms1spectra{i})
%here we specify the mass window
   if ms1spectra{i}(j,1)<500 || ms1spectra{i}(j,1)>1000;
       ms1spectra{i}(j,:)=0;
   end
end

index=find(ms1spectra{i}(:,1)==0);
ms1spectra{i}(index,:)=[];
end
toc


mat=[0 0 0 0 0];

tic
for i=1:length(orbtime)
scan=cell2mat(ms1spectra(i));
scan=double(scan);

for j= 1:length(ms1spectra{i})
    mm=ms1spectra{i}(j,1);mi=ms1spectra{i}(j,2);
    
    if mi>0
    mmin=mm+1.0034-0.01;mmax=mm+1.0034+0.01;
    imin=0.2*mi*(mm*0.5/12*1.0816)/100; imax=1.8*mi*(mm*0.5/12*1.0816)/100;
%here we specify the tolerance in delta m,
%and in the ratio between M and M+1 peak
    
    index=find(scan(:,1)>=mmin & scan(:,1)<=mmax & scan(:,2)>=imin & scan(:,2)<=imax &scan(:,2)>0);
    ind=length(index);
%here we look for the isotopologue of M and M+1 peak
    if ind>0
      smat=zeros(ind,5);
      smat(:,1)=mm;  smat(:,2)=mi; 
      smat(:,3)=ms1spectra{i}(index,1); smat(:,4)=ms1spectra{i}(index,2); 
      smat(:,5)=orbtime(i); 
      mat=[mat;smat];
%here we 'virtually label' the M peak as blue and the M+1 peak as red
    end
   end
 end
end
toc

mat(1,:)=[];
lowmass=mat(:,1);



tic
%here we specify the mass window again
mass=[500:0.01:1000];

m=0;
for i=1:length(mass) 
 binmin=mass(i)-0.005;binmax=mass(i)+0.005;
 binbar=find(lowmass>=binmin & lowmass<=binmax);
 bb=length(binbar);
 
%here we examine each EIC to see how many 'blue' isotopologues are there,
%and we specify how many datapoints need to be there for an EIC to pass the filter
if bb>40
    m=[m,mass(i)];
 end
    
end
m(1)=[];
toc


for i=1:length(m)
figure;
fechcts=rawEIC(m(i),orbtime,ms1spectra);
plot(orbtime,fechcts,':b','LineWidth',4);hold on

fechcts=rawEIC(m(i)-1.995,orbtime,ms1spectra)*15;
plot(orbtime,fechcts,':r','LineWidth',4);hold on

title(m(i)); legend('56Fe','54Fe'); legend boxoff

end