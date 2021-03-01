%Function for CMA-T, a feature detection algorithm 
%it classifies each EIC as either a target ion or a noise
%we use a mixture of 10 nM Fe-desferrichrome and 10 nM Fe-desferrioxamine E for demo
%the sample also contains 10 nM Ga-desferrioxamine E as internal standard
%the LC-ICPMS trace can be found in the same folder as data for demo
%here we show how to find the mass of Fe-ferrichrome from LC-ESIMS data
%this algorithm calls for three functions: 'multiMSalignGa', 'rawEIC', and 'accumulation'
clear; close all; clc

%you need to re-specify dir
cd /Users/jingxuanli/Dropbox/matlab

%the esims data is a .raw file, we need to convert that to mzXML using MSConvert
%VERY IMPORTANTLY, when you convert the raw data, please follow the figure
%https://www.mathworks.com/matlabcentral/answers/277803-error-using-mzxml-read
%TWO NOTES: select mslevel 1-2, not 1-1; select 30000 most intense counts if needed
Orbifile='StdMix_10nM_20200517191203.mzXML';

%this function import mzXML data in a form managable by Matlab
mzXMLstruct = mzxmlread(Orbifile);

%now we need to align the esi dataset towards icp dataset
%import the icp dataset first
Icpfile='03_200819_GP15_Sta12RE_StdMix_10nM';
[data,txt,raw]=xlsread(Icpfile);

%the alignment is done using the function 'multiMSalignGa'
%it finds the retention time of Ga-DFE on both MS and align them
%'aa' means 'after alignment'
%the only difference between 'esiaa' and 'mzXMLstruct' is the RT column
esiaa=multiMSalignGa(mzXMLstruct, raw);

%extract the data for each MS scan to an accessible form
%WHEN YOU HAVE NO internal standard, don't do the multiMSalignGa, instead do
%[ms1spectra, orbtime] = mzxml2peaks(mzXMLstruct);
%so the data can be extracted, though not aligned
[ms1spectra, orbtime] = mzxml2peaks(esiaa);

%save two important items: orbtime and ms1spectra (the MS scan at each orbtime)
matname = ['/Users/jingxuanli/Dropbox/matlab/_orbtimeStdMix10nM' ];
save(matname,'orbtime');
matname = ['/Users/jingxuanli/Dropbox/matlab/_ms1spectraStdMix10nM' ];
save(matname,'ms1spectra');





%now that the mzXML has been imported and aligned, feature detection begins
clear; close all; clc

  load('_orbtimeStdMix10nM.mat');
load('_ms1spectraStdMix10nM.mat');

%lets look for the ferrichrome (741.24, but pretending we dont know the mass)
%we call it the first Fe peak as we pretend that its identity is unknown
%specify where the ICP peak center is, should be available via ICP 56Fe
T=2295;

%cut the esi dataset around the peak, so that the run is faster
%we look for the mass of the first Fe peak within this 200 s time window in esi dataset 
tfake=find(orbtime<2195 | orbtime>2395);
orbtime(tfake)=[];
ms1spectra(tfake)=[];

tic

%go over from m/z 500 to m/z 1000, to classify each EIC as a target ion or noise
%the mass window could be customised, just like the time window of interest
%using 0.01 for binning, can be refined to 0.001, depending on resolution
%finer search would significantly increase the run time
mass=[500:0.01:1000];

for i=1:length(mass) 
 orbt=orbtime;
 
 %rawEIC is the function that bins the data within +/- 0.005 m/z, 
%which is why the step is 0.01 in for loop; the binning can be changed in the rawEIC function
 orbcts=rawEIC(mass(i),orbtime,ms1spectra);
 
 frame=[orbt orbcts];
 %for each of the mass, generate the EIC 
 
 %if there is too few datapoints on the EIC, we regard the EIC as 'NO PEAK'
 %this will make the overall run faster, as less EIC would need to be examined
 %the threshold of 30 SHALL VARY BETWEEN DATASET, 
 %consult peak for internal standard to see how many datapoints we expect to be in the EIC
 pos=length(find(orbcts>0));
 if pos>30
 
 %Accumulcation is the function culculating delta T from T5 (5% integral) to T95 (95% integral)
 %and normalize that to the time window (200 s here) for the t_percentage
    [tpeak( i,1),t, EnerIN,t_percentage(i,1)]= Accumulation(frame);
 
 %this is a 'noiseremove' step
 %tpeak is the average between T5 and T95, defined by esi data
 %T is the peak center defined by icp data
 %if the peak is perfectly symmatric, tpeak should be identical to the center of ICP peak (T)
 %in reality we always have tailing, tpeak is few seconds off EIC summit
 %we could also have shift in retention time between lc-icpms and lc-esims
 %that is why we have to allow for an offset, such as 60 s in this example
  %if the tpeak is too far from ICP summit, T, this EIC is no good
 %the threshold needs to be increased if no alignment is performed        
 if abs(tpeak(i,1)-T)>60
             t_percentage(i,1)=nan;tpeak(i,1)=nan;
             end

 else t_percentage(i,1)=nan;tpeak(i,1)=nan;
 end

end
toc

%if delta T spanning the 5-95% intergral < 20% of total T (200 s), we register the EIC as a targe ion
%OUTPUTTING 'm' as the list of target ions
index=find(t_percentage<20);
m=mass(index);





%SECOND SECTION OF ISOTOPELOGUE OVERLAY STARTS HERE
%we assume each 'm' output by CMA-T represents an 56FeL compound, and look for its hypothetically existing 54FeL version
%this is an examination step to see which of the target ions contain Fe

for i=1:length(m)
figure;
fechcts=rawEIC(m(i),orbtime,ms1spectra);
plot(orbtime,fechcts,':b','LineWidth',4);hold on

%go over each qualified mass, assume that is 56FeL
%plot the EIC of it, look for the EIC of 54FeL version, overlay
fechcts=rawEIC(m(i)-1.995,orbtime,ms1spectra)*15;
plot(orbtime,fechcts,':r','LineWidth',4);hold on

title(m(i)); legend('56Fe','54Fe'); legend boxoff
set (gca, 'fontsize',16);

end     
