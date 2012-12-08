% function kerncompmaster(cellid,batch)
%
% single queue entry wrapper for kerncomp4.m
% attentional modulation testing
%
% basic flow: 
% 1. read specified "kerncomp" entry in tQueue
% 2. load up parameters and call kerncomp.m
%
% created SVD 07/10/03  branched off of kerncompqueue.m
%
function kerncompmaster(cellid,batch)

disp('kerncompmaster.m: STARTING');

% set rand seed for the hell of it
rand('state',sum(100*clock));

global BATQUEUEID
dbopen;

sql=['SELECT * from sRunData WHERE cellid="',cellid,'"',...
     ' AND batch=',num2str(batch)];
rundata=mysql(sql);

fprintf('CELLID=%s BATCH=%d\n',cellid,rundata.batch);

[cellfiledata,times,params]=cellfiletimes(cellid,rundata.batch);

params.times=times;
params.cellid=cellid;

% predbatch--batches containing other stim classes
if isempty(params.predbatch),
   params.predbatch={params.id};
else
   params.predbatch=strsep(params.predbatch,',');
end

params.batch=rundata.batch;

% these entries in params need to be parsed
params.resploadparms=strsep(params.resploadparms,',');
params.respfilterparms=strsep(params.respfilterparms,',');   
params.stimloadparms=strsep(params.stimloadparms,',');
params.stimfilterparms=strsep(params.stimfilterparms,',');
if ~isnumeric(params.sffiltsigma),
   params.sffiltsigma=strsep(params.sffiltsigma,',');
end

params.docellfit2=0;
params.shrinkage=1;   % rather than just thresholding
params.repexclude=0;

params.zipoutfile=1;
params.outfile=[rundata.respath,rundata.resfile];
%params.outfile='/auto/k5/david/test.mat';
%params.meansub=0;

if length(params.parmstring)>0,
   eval(char(params.parmstring));
end

postanal='kerncomp';

% this is for the old paired attention comparison. may want
% to go back
%resampcount=1;
%noisecount=300;
%kerncomp3;  % kerncomp or kerncomp2 or kerncomp3 or kerncomp4?

params.noisecount=500;
%params.meansub=0;

if length(params.parmstring)>0,
   eval(char(params.parmstring));
end

% kerncomp or kerncomp2 or kerncomp3 or kerncomp4 or kerncomp5?
if isfield(params,'altkva'),
   eval(params.altkva);
else
   kernvsatt2;
end


% save the results
daterun=date;
skernfile=rundata.kernfile;
sfullfile=[rundata.respath,skernfile];
disp(sprintf('Saving kerncomp results to %s...',skernfile));

% then save everything else (inclusive to avoid future screw-ups
% involving accidentally not saving new important stuff)
tpath=getenv('TMP');
save([tpath,'/',skernfile]);

% compress and copy over network if necessary
disp(['alt save! tpath=',tpath]);
disp(['gzip -f ',tpath,'/',skernfile]);
unix(['gzip -f ',tpath,'/',skernfile]);
disp(['mv -f ',tpath,'/',skernfile,'.gz ',sfullfile,'.gz']);
unix(['mv -f ',tpath,'/',skernfile,'.gz ',sfullfile,'.gz']);
%disp('deleting temp file');
%delete([tpath,'/',skernfile]);

% construct a string to dump to the results table for later summarizing
sres=sprintf('preddata(%d).cellid=''%s'';',rundata.id,cellid);
sres=sprintf('%s preddata(%d).predxc=%s;',...
             sres,rundata.id,mat2string(predxc));
sres=sprintf('%s preddata(%d).pxc=%s;',...
             sres,rundata.id,mat2string(pxc));
sres=sprintf('%s preddata(%d).ncount=%s;',...
             sres,rundata.id,mat2string(ncount));
sres=sprintf('%s preddata(%d).outnlmode=%d;',...
             sres,rundata.id,OUTNLMODE);
if exist('predinf','var'),
   sres=sprintf('%s preddata(%d).predinf=%s;',...
                sres,rundata.id,mat2string(predinf));
end

% from kerncompX:
if exist('nondampedcount','var'),
   sres=sprintf('%s preddata(%d).nondampedcount=%s;',...
                sres,rundata.id,mat2string(nondampedcount));
end
% from kernvsatt:
if exist('randxc','var'),
   sres=sprintf('%s preddata(%d).randxc=%s;',...
                sres,rundata.id,mat2string(randxc));
end
if exist('predxcpair','var'),
   sres=sprintf('%s preddata(%d).predxcpair=%s;',...
                sres,rundata.id,mat2string(predxcpair));
   sres=sprintf('%s preddata(%d).pxcpair=%s;',...
                sres,rundata.id,mat2string(pxcpair));
end

% determine whether there's already a results record and delete
% it so that it can be replaced.  maybe this should be
% superceded by an archiving scheme someday?
mysql(['DELETE FROM sResults WHERE runid=',num2str(rundata.id)]);

sqlinsert('sResults',...
          'runid',rundata.id,...
          'batch',rundata.batch,...
          'matstr',sres);




