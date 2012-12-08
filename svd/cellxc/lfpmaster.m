% function lfpmaster(siteid,batch)
%
% wrapper for lfp analysis program
%
% basic flow: 
% 1. figure out current queue id from enviroment (set by dbqueuemaster)
% 2. load up parameters from sBatch
% 3. run the analysis
%
% created SVD 10/16/07 -- ripped off of cellxcmaster
%
function lfpmaster(siteid,batch)

disp('cellxcmaster.m:');

if ~exist('batch','var'),
   disp('syntax error: lfpmaster(siteid,batch) parameters required');
   return
end

% set rand seed for the hell of it
rand('state',sum(100*clock));

global BATQUEUEID

dbopen;

sql=['SELECT * from sRunData WHERE cellid="',siteid,'"',...
     ' AND batch=',num2str(batch)];
rundata=mysql(sql);

if length(rundata)==0,
   fprintf('(siteid,batch)=(%s,%d) not found!\n',cellid,batch);
   return
end

fprintf('SITEID=%s BATCH=%d\n',siteid,batch);




% figure out what files to use for what stage of the analysis
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
params.meansub=1;

params.zipoutfile=0;
params.outfile=[rundata.respath,rundata.resfile];

if onseil,
   params.finaloutfile=params.outfile;
   params.outfile=strrep(params.outfile,'/auto/data/','/homes/svd/data/');
   [bb,pp]=basename(params.outfile);
   unix(['mkdir -p ',pp]);
end

if length(params.parmstring)>0,
   eval(char(params.parmstring));
end

% actually do the strf estimation:
cellxcnodb(params);

disp('.');

z=zload([params.outfile,'.gz']);

disp('.');

% construct a string to dump to the results table for later summarizing
sres=sprintf('preddata(%d).cellid=''%s'';',rundata.id,cellid);
sres=sprintf('%s preddata(%d).predxc=%s;',...
             sres,rundata.id,mat2string(z.predxc));
sres=sprintf('%s preddata(%d).predp=%s;',...
             sres,rundata.id,mat2string(z.predp));
sres=sprintf('%s preddata(%d).prederr=%s;',...
             sres,rundata.id,mat2string(z.prederr));
sres=sprintf('%s preddata(%d).predinf=%s;',...
             sres,rundata.id,mat2string(z.predinf));
sres=sprintf('%s preddata(%d).predfix=%s;',...
             sres,rundata.id,mat2string(z.predfix));
sres=sprintf('%s preddata(%d).predfixerr=%s;',...
             sres,rundata.id,mat2string(z.predfixerr));
sres=sprintf('%s preddata(%d).sigfit=%s;',...
             sres,rundata.id,mat2string(z.sigfit));
sres=sprintf('%s preddata(%d).sfsfit=%s;',...
             sres,rundata.id,mat2string(z.sfsfit));
sres=sprintf('%s preddata(%d).expxc=%s;',...
             sres,rundata.id,mat2string(z.expxc));

% determine whether there's already a results record and delete
% it so that it can be replaced.  maybe this should be
% superceded by an archiving scheme someday?
sql=['SELECT * FROM sResults WHERE runid=',num2str(rundata.id)];
resdata=mysql(sql);
if length(resdata)>0,
   mysql(['DELETE FROM sResults WHERE id=',num2str(resdata(1).id)]);
end

disp('.');

sqlinsert('sResults',...
          'runid',rundata.id,...
          'batch',rundata.batch,...
          'matstr',sres);

if rundata.batch>=51 & rundata.batch<=62,
   xcfittune(cellid,rundata.batch);
end

if onseil,
   ['scp ',params.outfile,' svd@bhangra.isr.umd.edu:',params.finaloutfile]
   unix(['scp ',params.outfile,' svd@bhangra.isr.umd.edu:',params.finaloutfile]);
end
   
% record that we're done with this queue entry
dbsetqueue(BATQUEUEID,1000,1);

disp('.');

