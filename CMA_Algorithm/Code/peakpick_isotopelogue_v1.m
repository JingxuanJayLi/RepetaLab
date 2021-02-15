%Function for CMA-T, which calls the previous three function
%this script includes 2 sections: peak picking and isotopelogue overlay
%therefore, it should be used for Fe, Cu, etc
%in contrast to the 'bottom-up' searching algorithm by Boiteau et al 2015, 
%this is a 'top-down' searching, made possible by the peak picking,
%which diminishes the dataset from 10000s of frames (EICs) to 100s of frames
%I will use the StdMix as demo

clear; close all; clc
cd /Users/jingxuanli/Dropbox/matlab%you need to re-specify dir
Orbifile='StdMix_10nM_20200517191203.mzXML';
mzXMLstruct = mzxmlread(Orbifile);
%VERY IMPORTANTLY, when you convert the raw data, please follow the figure
%https://www.mathworks.com/matlabcentral/answers/277803-error-using-mzxml-read
%TWO NOTES: select mslevel 1-2, not 1-1; select 30000 most intense counts if needed

Icpfile='03_200819_GP15_Sta12RE_StdMix_10nM';
[data,txt,raw]=xlsread(Icpfile);
esiaa=multiMSalignGa(mzXMLstruct, raw);
%this is to align the esi dataset towards icp dataset, 'aa' means 'after alignment'
%using the function 'multiMSalign' I wrote, it will be included in the folder
%the only difference between 'esiaa' and 'mzXMLstruct' is the RT column

[ms1spectra, orbtime] = mzxml2peaks(esiaa);
%extract the data for each scan (retention time) to an accessible form
%WHEN YOU HAVE NO B12, don't do the multiMSalign, instead do
%[ms1spectra, orbtime] = mzxml2peaks(mzXMLstruct);


matname = ['/Users/jingxuanli/Dropbox/matlab/_orbtimeStdMix10nM' ];
save(matname,'orbtime');
matname = ['/Users/jingxuanli/Dropbox/matlab/_ms1spectraStdMix10nM' ];
save(matname,'ms1spectra');

%save two important items: orbtime and ms1spectra (the scan at each orbtime)

clear; close all; clc

  load('_orbtimeStdMix10nM.mat');
load('_ms1spectraStdMix10nM.mat');

%for the test data, the mzXML data is too large, 
%so we previde the orbtimeStdMix10nM.mat and ms1spectraStdMix10nM.mat,
%together with the ICPMS data. 
%we are happy to provide the raw mzXML data upon contact

T=2295;
%lets look for the ferrichrome (741.24, but pretending we dont know the mass)
%specify where the ICP peak center is, should be available via ICP 56Fe
%this is the only place icpdata is used
%the icpdata is not available but I know where the peak is

tfake=find(orbtime<2195 | orbtime>2395);
orbtime(tfake)=[];
ms1spectra(tfake)=[];
%cut the esi dataset around the peak, so that the run is faster

tic
mass=[500:0.01:1000];
%go over from m/z 500 to m/z 1000, can be expanded from 200 to 1800
%using 0.01 for binning


for i=1:length(mass) 
 orbt=orbtime;
 orbcts=rawEIC(mass(i),orbtime,ms1spectra);
 
%rawEIC is a function I wrote, to bin the data within +/- 0.005 m/z, 
%which is why the step is 0.01 in for loop; the binning can be changed in the rawEIC function

 frame=[orbt orbcts];
 %for each of the mass, generate the EIC 
 
 pos=length(find(orbcts>0));
 if pos>30
 %if there is too few datapoints on the EIC (<60), we regard the EIC as 'NO PEAK',
 %this will make the overall run faster
 %the threshold of 60 SHALL VARY BETWEEN DATASET, 
 %consult B12 peak (or Ga-DFE peak) to see how many datapoints we expect to be in the EIC
    [tpeak( i,1),t, EnerIN,t_percentage(i,1)]= Accumulation(frame);
 %Accumulcation is the function culculating delta T from T5 to T95
 %and compare that to the time window
             if abs(tpeak(i,1)-T)>60
 %this is a 'noiseremove' step
 %tpeak is the average between T5 and T95
 %if the peak is perfectly symmatric, tpeak should be identical to the center of ICP peak (T)
 %in reality we always have tailing, tpeak is few seconds off EIC summit
 %it does not matter, as we allow for a 20s offset
             t_percentage(i,1)=nan;tpeak(i,1)=nan;
             end
 %if the tpeak is too far from ICP summit, T, this EIC is no good
 %the threshold needs to be increased if no B12 alignment is performed
 else t_percentage(i,1)=nan;tpeak(i,1)=nan;
 end

end
toc


index=find(t_percentage<20);
m=mass(index);
%if delta T spanning the 5-95% energy < 30% of total T, we register that mass
%we believe these are real peaks and ignore other masses
%FIRST SECTION OF PEAK PICKING FINISHES HERE, OUTPUTTING 'm' 


%SECOND SECTION OF ISOTOPELOGUE OVERLAY STARTS HERE

for i=1:length(m)
figure;
fechcts=rawEIC(m(i),orbtime,ms1spectra);
plot(orbtime,fechcts,':b','LineWidth',4);hold on

fechcts=rawEIC(m(i)-1.995,orbtime,ms1spectra)*15;
plot(orbtime,fechcts,':r','LineWidth',4);hold on
%go over each qualified mass, assume that is 56FeL, 
%plot the EIC of it, look for the EIC of 54FeL version, overlay
%IT IS IMPORTANT THAT WHEN DOING COPPER, we are assuming it to be 63CuL
%so we need to look up for 65Cu and scale the 2nd EIC up

title(m(i)); legend('56Fe','54Fe'); legend boxoff
set (gca, 'fontsize',16);


end     
