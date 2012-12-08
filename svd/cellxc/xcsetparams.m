% function params=xcsetparams(cellid,batchid);
%
% load params structure for a given cell/batch analysis.
% this consists of:
% 1. call cellfiledata to find relevant files in sCellFile, figure
%    out a times structure for exp/fit/val, and get the raw
%    information from batchdata.
% 2. construct params.stimfiles, .respfiles out of cellfiledata
% 3. set a few other random parameters and expand strings in
%    batchdata into cell arrays for things like
%    params.stimfilterparms{}
%
% CREATED SVD 5/29/03 ripped off of cellxcqueue.m
% 
function params=xcsetparams(cellid,batchid);

% figure out what files to use for what stage of the analysis
[cellfiledata,times,batchdata]=cellfiletimes(cellid,batchid);

% fill up params structure for passing to cellxcnodb
params=batchdata;
params.times=times;
params.cellid=cellid;
params.batch=batchid;

params.stimfiles={};
params.respfiles={};
resplens=zeros(length(cellfiledata),1);
stimpix=zeros(length(cellfiledata),1);
params.stimcrfs=zeros(length(cellfiledata),1);
totlen=0;
for ii=1:length(cellfiledata),
   params.stimfiles{ii}=[cellfiledata(ii).stimpath,...
                    cellfiledata(ii).stimfile];
   params.respfiles{ii}=[cellfiledata(ii).path,...
                    cellfiledata(ii).respfile];
   resplens(ii)=cellfiledata(ii).resplen;
   tstimpix=strsep(cellfiledata(ii).stimiconside,',');
   if length(tstimpix)>0,
      stimpix(ii)=tstimpix{1};
   end
   if cellfiledata(ii).stimfilecrf>0,
      params.stimcrfs(ii)=cellfiledata(ii).stimfilecrf;
   else
      sql=['SELECT * FROM gCellMaster WHERE cellid="',cellid,'"'];
      celldata=mysql(sql);
      params.stimcrfs(ii)=stimpix(ii)./celldata.rfsize;
   end
end

% predbatch--batches containing other stim classes
if isempty(batchdata.predbatch),
   params.predbatch={batchdata.id};
else
   params.predbatch=strsep(batchdata.predbatch,',');
end

% these entries in batchdata need to be parsed
params.resploadparms=strsep(batchdata.resploadparms,',');
params.respfilterparms=strsep(batchdata.respfilterparms,',');   
params.stimloadparms=strsep(batchdata.stimloadparms,',');
params.stimfilterparms=strsep(batchdata.stimfilterparms,',');
if ~isnumeric(params.sffiltsigma),
   params.sffiltsigma=strsep(batchdata.sffiltsigma,',');
end

params.docellfit2=0;
params.shrinkage=1;   % rather than just thresholding
if ismember(params.batch,[66 70 72 23 29 32 37]),
   params.repexclude=1;
else
   params.repexclude=0;
end

params.zipoutfile=1;
params.maxlag=[params.minlag params.maxlag];

if length(batchdata.parmstring)>0,
   eval(char(batchdata.parmstring));
end

