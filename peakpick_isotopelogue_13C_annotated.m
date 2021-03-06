%Function for CMA-C, which shares the data import and alignment as CMA-T
%here we start with the organized data, which was generated by line 1-43 in 'peakpick_isotopelogue_v1'
clear; close all; clc
  load('_orbtimeStdMix10nM.mat');
load('_ms1spectraStdMix10nM.mat');

%here we specify the time window 
tfake=find(orbtime<2195 | orbtime>2395);
orbtime(tfake)=[];
ms1spectra(tfake)=[];

tic
%here we cut the dataset and discard the data beyond mass window of interest
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
%here we go over each MS scan, and in each scan we look for pairs of peaks representing 12C and 13C isotopologue
for i=1:length(orbtime)

%extract the data for each scan
scan=cell2mat(ms1spectra(i));
scan=double(scan);

%for each ion in a MS scan, we extract its m/z and intensity
%then we look for its 13C isotopologue within the same MS scan
%two ions are regarded as 12C-13C isotopologue only if their masses and itensities follow specific pattern
%in summary, here we are trying to find all pairs of 12C-13C isotopologues for metabolites in each scan

for j= 1:length(ms1spectra{i})
    mm=ms1spectra{i}(j,1);mi=ms1spectra{i}(j,2);
    
%here we specify the pattern in delta m, and the ratio between M and M+1 peak
    if mi>0
    mmin=mm+1.0034-0.01;mmax=mm+1.0034+0.01;
%we assume a carbon content of 50%±40% in the molecule, and estimate the number of carbon atoms in the molecule
%therefore, isotopologue representing most of the organic compounds would be found by CMA-C
%the compounds with a carbon content < 10% or > 90% would be discarded
%mi represents the intensity or the 12C ion
%mm represents the m/z of the 12C ion
%0.5 represents a carbon content of 50%
%0.2*0.5 represents a carbon content of 10%, 1.8*0.5 represents a carbon content of 90%
%12 represents the atomic mass of 12C, 1.0816 represents the relative abundance of 13C
%for example, for an ion which has m/z of 600, and intensity of 1000, the intensity of its M+1 peak should be 54 to 487 
%the tolerance is large here as we are trying to be conservative
%we try to avoid rejection of true positives that have a carbon content significantly different from 50%

imin=0.2*mi*(mm*0.5/12*1.0816)/100; imax=1.8*mi*(mm*0.5/12*1.0816)/100;

%here we look for the isotopologue of M and M+1 peak
    index=find(scan(:,1)>=mmin & scan(:,1)<=mmax & scan(:,2)>=imin & scan(:,2)<=imax &scan(:,2)>0);
    ind=length(index);

%here we 'virtually label' the M peak as 'blue' and the M+1 peak as 'red' dataopoints
    if ind>0
      smat=zeros(ind,5);
      smat(:,1)=mm;  smat(:,2)=mi; 
      smat(:,3)=ms1spectra{i}(index,1); smat(:,4)=ms1spectra{i}(index,2); 
      smat(:,5)=orbtime(i); 

      mat=[mat;smat];
    end
   end
 end
end
toc

mat(1,:)=[];
lowmass=mat(:,1);

%here we specify the mass window again
tic
mass=[500:0.01:1000];

m=0;
for i=1:length(mass) 
 binmin=mass(i)-0.005;binmax=mass(i)+0.005;
 binbar=find(lowmass>=binmin & lowmass<=binmax);
 bb=length(binbar);
 
%here we examine each EIC to see how many 'blue' datapoints are there
%and we specify how many datapoints need to be there for an EIC to pass the filter
if bb>40
    m=[m,mass(i)];
 end
    
end
m(1)=[];
toc

%SECOND SECTION OF ISOTOPELOGUE OVERLAY STARTS HERE
%we assume each 'm' output by CMA-C represents an 56FeL compound, and look for its hypothetically existing 54FeL version
%this is an examination step to see which of the target ions contain Fe
for i=1:length(m)
figure;
fechcts=rawEIC(m(i),orbtime,ms1spectra);
plot(orbtime,fechcts,':b','LineWidth',4);hold on

fechcts=rawEIC(m(i)-1.995,orbtime,ms1spectra)*15;
plot(orbtime,fechcts,':r','LineWidth',4);hold on

title(m(i)); legend('56Fe','54Fe'); legend boxoff

end