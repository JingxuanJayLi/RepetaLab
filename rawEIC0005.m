%this function is equivalent to R function 'rawEIC', basically generating EIC for a given m/z
%with a tolerance of +/- 0.005, of course it can be speficied using other tolerance
function EICcounts = rawEIC(EICmz,EICtime,EICspectra)
%inputs: m/z; the time vector for all scans; and the peak vector
EICcounts=zeros(length(EICtime),1);
mzmin=EICmz-0.0005;mzmax=EICmz+0.0005;
    for i=1:length(EICtime)
   EIscan = cell2mat(EICspectra(i)); 
    EIN=find(mzmin<=EIscan(:,1) & EIscan(:,1)<=mzmax);
EICcounts(i)=sum(EIscan(EIN,2));
    end
end