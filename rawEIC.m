%Function rawEIC

%This function generates EIC for a given m/z
%with a tolerance of +/- 0.005, it can be speficied using other tolerance
function EICcounts = rawEIC(EICmz,EICtime,EICspectra)
%inputs: EICmz: m/z; EICtime: time vector for all scans; EICspectra: intrensity-m/z pair for all scans
EICcounts=zeros(length(EICtime),1);
mzmin=EICmz-0.005;mzmax=EICmz+0.005;
    for i=1:length(EICtime)
   EIscan = cell2mat(EICspectra(i)); 
    EIN=find(mzmin<=EIscan(:,1) & EIscan(:,1)<=mzmax);
EICcounts(i)=sum(EIscan(EIN,2));
    end
end
