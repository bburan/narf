% function cellxcmaster(cellid,batch)
%
% wrapper for cellxcnodb.m reverse correlation routine
%
% basic flow: 
% 1. figure out current queue id from enviroment (set by dbqueuemaster)
% 2. load up parameters from sBatch 
% 3. cellxcno db.m
%
% created SVD 6/7/03 -- ripped off of cellxcqueue
%
function cellxcmaster(cellid,batch)

disp('cellxcmaster.m:');
baphy_path_analysis;

if ~exist('batch','var'),
   disp('syntax error: cellxcmaster(cellid,batch) parameters required');
   return
end

% set rand seed for the hell of it
rand('state',sum(100*clock));

global BATQUEUEID

dbopen;

sql=['SELECT * from sRunData WHERE cellid="',cellid,'"',...
     ' AND batch=',num2str(batch)];
rundata=mysql(sql);

if length(rundata)==0,
   fprintf('(cellid,batch)=(%s,%d) not found!\n',cellid,batch);
   keyboard
   return
end

fprintf('CELLID=%s BATCH=%d\n',cellid,batch);

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

params.showres=1;
params.zipoutfile=0;
params.outfile=[rundata.respath,rundata.resfile];

ss_stat=onseil;
if ss_stat,
   params.finaloutfile=params.outfile;
   if ss_stat==1,
       params.outfile=strrep(params.outfile,'/auto/data/','/homes/svd/data/');
   elseif ss_stat==2.
       % keep the same;
       disp('keeping output file same location at OHSU');
   end
end

% make the final path if it doesn't exist yet
[bb,pp]=basename(params.outfile);
unix(['mkdir -p ',pp]);

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
   w=unix(['scp ',params.outfile,' svd@bhangra.isr.umd.edu:',params.finaloutfile]);
   if ~w && onseil==1,
      disp('copied from seil sucessfully, deleting local copy to save space');
      unix(['rm -f ',params.outfile]);
   end
   
end


% record that we're done with this queue entry
dbsetqueue(BATQUEUEID,1000,1);

disp('.');

