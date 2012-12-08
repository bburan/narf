% function reconmaster(cellid,batch)
%
% wrapper for cellxcnodb.m reverse correlation routine
%
% basic flow: 
% 1. figure out current queue id from enviroment (set by dbqueuemaster)
% 2. load up parameters from sBatch 
% 3. cellxcno db.m
%
% created SVD 6/7/11 -- ripped off of cellxcmaster
%
function reconmaster(cellid,batch)

disp('reconmaster.m:');
baphy_set_path

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
   return
end

fprintf('CELLID=%s BATCH=%d\n',cellid,batch);


% figure out what files to use for what stage of the analysis
[cellfiledata,times,params]=cellfiletimes(cellid,rundata.batch);

params.cellid=cellid;
params.times=times;

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

params.smoothwin=getparm(params,'smoothwin',5);
params.mlindep=getparm(params,'mlindep',0);

if strcmp(params.stimloadcmd,'loadsiteraster'),
   % figure out unique units that exist across all files
   for ss=1:length(params.stimfiles),
      bb=basename(params.stimfiles{ss});
      sql=['SELECT channum,unit FROM sCellFile WHERE respfile="',bb,'";'];
      fdata=mysql(sql);
      if ss==1,
         unitset=[cat(1,fdata.channum).*10+cat(1,fdata.unit)];
      else
         unitset=intersect(unitset,[cat(1,fdata.channum).*10+cat(1,fdata.unit)]);
      end
   end
   if isfield(params.stimloadparms,'lfp') &&...
         params.stimloadparms.lfp>0,
      cc=unique(floor(unitset/10));
      params.stimloadparms{1}.channel=cc;
      params.stimloadparms{1}.unit=ones(size(cc));
   else
      params.stimloadparms{1}.channel=floor(unitset/10);
      params.stimloadparms{1}.unit=mod(unitset,10);
   end
end

params.keepneg=1;
params.docellfit2=0;
params.shrinkage=1;   % rather than just thresholding
params.repexclude=0;
if strcmp(params.altcore,'mlrcore'),
   params.meansub=1;
else
   params.meansub=1;
end
params.showres=0;

params.zipoutfile=0;
params.outfile=[rundata.respath,rundata.resfile];

params.allpaircounts=getparm(params,'allpaircounts',0);
params.adddelay=0;

if params.allpaircounts==1,
   % test subsets of cell pairs from 2 to total
   global SUBXC
   SUBXC=[];
   
   cset=params.stimloadparms{1}.channel;
   uset=params.stimloadparms{1}.unit;
   cellcount=length(cset);
   params.finaloutfile=params.outfile;      
   params.outfile=[tempdir basename(params.finaloutfile)];
   subxc=zeros(cellcount,3);
   for cc=2:(cellcount-1),
      fprintf('testing set size %d/%d:\n',cc,cellcount);
      ff=round(linspace(1,length(cset),cc));
      [cset(ff) uset(ff)]
      params.stimloadparms{1}.channel=cset(ff);
      params.stimloadparms{1}.unit=uset(ff);
      
      cellxcnodb(params);
      
      tr=load(params.outfile);
      subxc(cc,:)=tr.predxc;
   end
   
   SUBXC=subxc;
   params.outfile=params.finaloutfile;
   params.stimloadparms{1}.channel=cset;
   params.stimloadparms{1}.unit=uset;
elseif ismember(params.allpaircounts,[2 3]),
   %analyze each pair or triplet of cells in the site
   
   cset=params.stimloadparms{1}.channel;
   uset=params.stimloadparms{1}.unit;
   pairset=nchoosek(1:length(cset),params.allpaircounts);
   paircount=size(pairset,1);
   params.finaloutfile=params.outfile;      
   params.outfile=[tempdir basename(params.finaloutfile)];
   subxc=zeros(paircount,params.maxlag-params.minlag+1,length(params.predbatch));
   predres=[];
   for cc=1:paircount,
      fprintf('testing pair %d/%d\n',cc,paircount);
      ff=pairset(cc,:);
      [cset(ff) uset(ff)]
      params.stimloadparms{1}.channel=cset(ff);
      params.stimloadparms{1}.unit=uset(ff);
      
      cellxcnodb(params);
      
      % same important results from this loop
      tr=load(params.outfile);
      subxc(cc,:,:)=permute(tr.predxc,[2 3 1]);
      strf(cc)=tr.strf;
      predres=cat(2,predres,tr.predres);
      
   end
   
   predxc=subxc;
   params.outfile=params.finaloutfile;
   params.stimloadparms{1}.channel=cset;
   params.stimloadparms{1}.unit=uset;
   
elseif params.allpaircounts==4,
   %analyze each cell individually with 2 time lags
   
   cset=params.stimloadparms{1}.channel;
   uset=params.stimloadparms{1}.unit;
   setcount=length(cset);
   params.adddelay=1;
   params.finaloutfile=params.outfile;      
   params.outfile=[tempdir basename(params.finaloutfile)];
   subxc=zeros(setcount,params.maxlag-params.minlag+1,length(params.predbatch));
   predres=[];
   for cc=1:setcount,
      fprintf('testing cell %d/%d\n',cc,setcount);
      ff=cc;
      [cset(ff) uset(ff)]
      params.stimloadparms{1}.channel=cset(cc);
      params.stimloadparms{1}.unit=uset(cc);
      
      cellxcnodb(params);
      
      % same important results from this loop
      tr=load(params.outfile);
      subxc(cc,:,:)=permute(tr.predxc,[2 3 1]);
      strf(cc)=tr.strf;
      predres=cat(2,predres,tr.predres);
   end
   
   predxc=subxc;
   params.outfile=params.finaloutfile;
   params.stimloadparms{1}.channel=cset;
   params.stimloadparms{1}.unit=uset;
   
end

if onseil,
   params.finaloutfile=params.outfile;
   if onseil==1,
       params.outfile=strrep(params.outfile,'/auto/data/','/homes/svd/data/');
   end
end
   
% make the final path if it doesn't exist yet
[bb,pp]=basename(params.outfile);
unix(['mkdir -p ',pp]);

if length(params.parmstring)>0,
   eval(char(params.parmstring));
end
params.smoothwin=getparm(params,'smoothwin',5);

% actually do the strf estimation:
if params.allpaircounts==2,
   figure;
   clf
   imagesc(subxc(:,:));
   colorbar
   title(sprintf('%s - %d\n',cellid,batch));
   save(params.outfile,'strf','predxc','pairset','paircount','params',...
        'times','cellfiledata','predres');
elseif params.allpaircounts==3,
    disp('triplet analysis under construction');
    keyboard
elseif params.allpaircounts==4,
   figure;
   clf
   imagesc(subxc(:,:));
   colorbar
   title(sprintf('%s - %d\n',cellid,batch));
   
   save(params.outfile,'strf','predxc','setcount','params',...
        'times','cellfiledata','predres');
else
   cellxcnodb(params);
end

disp('Saving quick results to celldb...');

z=zload([params.outfile,'.gz']);

% construct a string to dump to the results table for later summarizing
sres=sprintf('preddata(%d).cellid=''%s'';',rundata.id,cellid);
sres=sprintf('%s preddata(%d).predxc=%s;',...
             sres,rundata.id,mat2string(z.predxc));
if params.allpaircounts<2,
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
end
if isfield(z,'subxc'),
   sres=sprintf('%s preddata(%d).subxc=%s;',...
                sres,rundata.id,mat2string(z.subxc));
end

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

