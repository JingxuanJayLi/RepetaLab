%Function multiMSalignGa

%the input of this function is the original esi and icp data
%the output of this function is an updated esi data
%the alignment internal standard here is Ga-desferrioxamine E
%but can be adjusted for any internal standards, such as B12

function esiAA=multiMSalign(ESIDATA, ICPDATA)

hd=ICPDATA(1,:);
Coindex=find(contains(hd,'69Ga'));
%it's named Coindex as this was originally designed to use B12
%the naming does not matter, as long as we select Ga columns from ICPMS
%and correct m/z for Ga-DFE from Orbitrap
ICPDATA(1,:)=[]; 
ICPDATA=cell2mat(ICPDATA);
%ICPDATA is a matrix without headers 
icpcotime=ICPDATA(:,Coindex(1));
icpcocps=ICPDATA(:,Coindex(2));
IcpCopeak=icpcotime(find(icpcocps == max(icpcocps)));
%extract ICPMS Ga data and find the retention time of the Ga-DFE peak in seconds

%when using msconvert, check figure at https://www.mathworks.com/matlabcentral/answers/277803-error-using-mzxml-read
[ms1specba, orbtimeba] = mzxml2peaks(ESIDATA);
[ms2specba, orbtime2ba] = mzxml2peaks(ESIDATA,'levels',2);
%ms1spectra/ms2spectra lists all peaks at any specific scan
%a peak list is a 2 column matrix with the mass/charge (mz) and ion intensity for every peak
%orbtimeba/orbtime2ba is a vector with the time of every scan
%note these four variables are before alignment

%Then need to find Ga-DFE (we originally used B12) peak center of its EIC within Orbitrap data
b12counts=rawEIC(667.2586,orbtimeba,ms1specba);
%use this FUNCTION to generate EIC for Ga-DFE peak, and then find peak center
OrbiCopeak=orbtimeba(find(b12counts == max(b12counts)));
offset=IcpCopeak-OrbiCopeak;
%use original esi and icp data to calculate the offset at peak center 

subplot (2,2,1);
%plot icp Ga before alignment
plot(icpcotime,icpcocps);hold on; plot([IcpCopeak IcpCopeak],[0 max(icpcocps)],'r');
title('ICPCo Before Alignment');set(gca,'xlim',[0 5000]);
word = num2str(IcpCopeak);
text(IcpCopeak+100,  max(icpcocps)/2 , word, 'color','r', 'fontsize',16);

%plot orbi Ga-DFE before alignment
subplot (2,2,3);
plot(orbtimeba,b12counts);hold on; plot([OrbiCopeak OrbiCopeak],[0 max(b12counts)],'r');
title('OrbB12 Before Alignment'); set(gca,'xlim',[0 5000]);
word = num2str(OrbiCopeak);
text(OrbiCopeak+100,  max(b12counts)/2 , word, 'color','r', 'fontsize',16);

%Ga alignment for MS1 data
%Need to do it for MS2 data, when working on MS2, not here
orbtimeba=orbtimeba+offset; OrbiCopeak=OrbiCopeak+offset;

subplot (2,2,2);
%plot icp Ga after alignment
plot(icpcotime,icpcocps);hold on; plot([IcpCopeak IcpCopeak],[0 max(icpcocps)],'g');
title('ICPCo After Alignment');set(gca,'xlim',[0 5000]);
word = num2str(IcpCopeak);
text(IcpCopeak+100,  max(icpcocps)/2 , word, 'color','g', 'fontsize',16);


%plot orbi Ga after alignment
subplot (2,2,4);
plot(orbtimeba,b12counts);hold on; plot([OrbiCopeak OrbiCopeak],[0 max(b12counts)],'g');
title('OrbB12 After Alighnment'); set(gca,'xlim',[0 5000]);
word = num2str(OrbiCopeak);
text(OrbiCopeak+100,  max(b12counts)/2 , word, 'color','g', 'fontsize',16);



%total scan = ms1 scan + ms2 scan
totalesiscan=length(orbtimeba)+length(orbtime2ba);
%go over each scan including ms1 and ms2, and change the time within esi data
for oo=1:totalesiscan
      chr=ESIDATA.scan(oo).retentionTime;
  chr(1)=[]; chr(1)=[]; chr(end)=[];
    timeaa=num2str(str2num(chr)+offset);   
     ESIDATA.scan(oo).retentionTime=strcat('P', 'T',timeaa,'S' ); 
end
   esiAA=ESIDATA;
%return an updated ESIMS dataset, the only changed variable is time 


end