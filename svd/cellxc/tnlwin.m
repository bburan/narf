% function tnlwin(cellid,batch)
%
% ripped off of tnl.m SVD 11/26
%
function tnlwin(cellid,batch)

fprintf('starting tnlwin.m:\n');

dbopen;

global BATQUEUEID

sql=['SELECT * from sRunData WHERE cellid="',cellid,'"',...
     ' AND batch=',num2str(batch)];
rundata=mysql(sql);
outfile=[rundata.respath,rundata.resfile];

% figure out what files to use for what stage of the analysis
[cellfiledata,times,params]=cellfiletimes(cellid,batch);

% fill up params structure
params.times=times;
params.cellid=cellid;
params.batch=batch;

params.repcount=cat(1,cellfiledata.repcount);
params.stimcrfs=zeros(length(cellfiledata),1);
resplens=cat(1,cellfiledata.resplen);
params.stimfiles={};
params.respfiles={};
stimpix=zeros(length(cellfiledata),1);
totlen=0;
for ii=1:length(cellfiledata),
   params.stimfiles{ii}=[cellfiledata(ii).stimpath,...
                    cellfiledata(ii).stimfile];
   params.respfiles{ii}=[cellfiledata(ii).path,...
                    cellfiledata(ii).respfile];
   
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
if isempty(params.predbatch),
   params.predbatch={params.id};
else
   params.predbatch=strsep(params.predbatch,',');
end

% these entries in batchdata need to be parsed
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
params.outfile=outfile;

if isfield(params,'maxlag') & length(params.maxlag)>=2,
   % do nothing, this is a good format for running cellxc
else
   params.maxlag=[getparm(params,'minlag',-6) getparm(params,'maxlag',13)];
end
params.tbinms=getparm(params,'tbinms',14);
params.meansub=getparm(params,'meansub',1);
params.smoothtime=getparm(params,'smoothtime',0);
params.sharpspacenorm=getparm(params,'sharpspacenorm',0);
params.nrandxcov=getparm(params,'nrandxcov',200);
params.maxlag=[0 0];

if length(params.parmstring)>0,
   eval(char(params.parmstring));
end

% 
% LOAD DATA
% 

% call xcfilefracs to get times
if isfield(params,'times'),
   times=params.times;
   
elseif params.fitboot,
   tparams=params;
   tparams.fitfrac=0;
   %tparams.predfrac=0;
   times=xcfilefracs(tparams);
   %times(1).stop(end)=times(1).stop(end)+1;
else
   times=xcfilefracs(params);
end

starttimes=times(1).start;
stoptimes=times(1).stop;
fitfile=times(2).fileidx;
fitstartframe=times(2).start;
fitstopframe=times(2).stop;
predfile=times(3).fileidx;
predstartframe=times(3).start;
predstopframe=times(3).stop;

attcount=size(times(1).start,2);  % ignore attcount>1 for the time being
attidx=1;

%
% LOAD FIT DATA
%

% skip stim filter... gets applied after variable window
stimfiltercmdbak=params.stimfiltercmd;
params.stimfiltercmd='';

xcloadfiles;

% trim fixations to shorter length if there are few longer ones
fixlen=size(resp,2);
fcount=sum(~isnan(resp),1);
lastgoodlag=max(find(fcount>2));
if lastgoodlag<fixlen,
   fixlen=lastgoodlag;
   resp=resp(:,1:fixlen);
end

% exclude fixations that are too short
if fixlen>10,
   okfixidx=find(sum(isnan(resp(:,1:5)),2)==0);
else
   okfixidx=find(sum(isnan(resp(:,1:2)),2)==0);
end

% use mikie's padding program to pad out acceptable fixations
resp(find(isnan(resp)))=-1;
resp(okfixidx,:)=padpsth(resp(okfixidx,:)')';
resp(find(resp==-1))=nan;

% save full response matrix for fitting output nls
respbak=resp;

% convert pfth to transient and sustained components
mresp=nanmean(resp);

g=[.1 .8 .1];
sm=conv2(mresp,g,'valid');

sustendidx=fixlen;
if fixlen>10,
   transstartidx=min(find(mresp>min(sm)+std(mresp).*1.5));
   transendidx=max([transstartidx+1 find(mresp>min(sm)+std(mresp).*1.5)]);
   if transendidx>fixlen-5,
      transendidx=fixlen-5;
   end
   if 0&fixlen>20,
      sustendidx=20;
      if transendidx>=sustendidx-1,
         transendidx=sustendidx-2;
      end
   end
   
else
   transstartidx=1;
   transendidx=3;
end

fprintf('\n transframes: %d-%d  sustframes: %d-%d\n\n',...
        transstartidx,transendidx,transendidx+1,sustendidx);

transresp=nansum(resp(:,transstartidx:transendidx)');
sresp=nansum(resp(:,transendidx+1:sustendidx)');
aresp=nansum(resp(:,transstartidx:sustendidx)');

nanidx=find(isnan(transresp) | isnan(sresp));
aresp(nanidx)=nan;
resp=[transresp' sresp' aresp'];

rsize=size(resp);
resplen=rsize(1);
respcount=rsize(2);

%
% LOAD VAL DATA
%
tpredstartframe=times(3).start;
tpredstopframe=times(3).stop;
tpredfile=times(3).fileidx;

[cdata.stim,cdata.resp]=xcloadstimresp(tpredfile,tpredstartframe,...
                                       tpredstopframe,params);

cdata.resp=cdata.resp(:,1:fixlen);
cdata.respbak=cdata.resp;

% use mikie's padding program
cdata.resp(find(isnan(cdata.resp)))=-1;
cdata.resp=padpsth(cdata.resp')';
cdata.resp(find(cdata.resp==-1))=nan;

vtransresp=nansum(cdata.resp(:,transstartidx:transendidx)');
vsresp=nansum(cdata.resp(:,transendidx+1:end)');
varesp=nansum(cdata.resp(:,transstartidx:end)');
nanidx=find(isnan(vtransresp) | isnan(vsresp));
varesp(nanidx)=nan;
cdata.resp=[vtransresp' vsresp' varesp'];

vallen=size(cdata.resp,1);

%
% DO RC
%

VCELLXC=3;
clear ESTIMATIONPHASE VALIDATIONPHASE
global ESTIMATIONPHASE VALIDATIONPHASE
ESTIMATIONPHASE=1;
VALIDATIONPHASE=0;

angles=[0 0 pi/4 pi/2 3*pi/4];
rats=[1 0.5 0.5 0.5 0.5];
stimbak=permute(reshape(stim,[size(stim,1) iconside]),[2 3 1]);
cstimbak=permute(reshape(cdata.stim,[size(cdata.stim,1) iconside]),[2 3 1]);

clear savestrf predres
xcsave=[];

for filtidx=1:length(angles),
   mask=gmask(iconside(1),angles(filtidx),rats(filtidx));
   stim=stimbak.*repmat(mask,[1 1 size(stimbak,3)]);
   stim=movpower(stim,0,0,0,1,0)';
   
   resp=aresp';
   rsize=size(resp);
   resplen=rsize(1);
   respcount=rsize(2);
   params.maxlag=[-2 1];
   resp0idx=-params.maxlag(1)+1;
   framecount=diff(params.maxlag)+1;
   
   cdata.resp=varesp';

   nnidx=find(~isnan(aresp'));
   vnnidx=find(~isnan(varesp'));
   
   % set parms to get xccore to run correctly
   spacecount=size(stim,2);
   firstseg=1;
   
   params.smoothtime=0;
   params.decorrtime=0;
   params.maxlag=[0 0];
   
   xccore;
   
   params.nloutparm=1;  % only do linear fit here
   xcfit;
   
   strf.mask=mask;
   savestrf(filtidx)=strf;
   xcsave=cat(3,xcsave,xc);
   
   mask=gmask(iconside(1),angles(filtidx),rats(filtidx));
   cdata.stim=cstimbak.*repmat(mask,[1 1 size(cstimbak,3)]);
   cdata.stim=movpower(cdata.stim,0,0,0,1,0)';
   cdata.resp=varesp';
   
   predres(filtidx)=xcval(strf(:),params,cdata);
end
predxc=cat(1,predres.predxc)

xcshowkern(savestrf);

disp('just did rc for each window size. need to deal with preds');

keyboard


xcfixcore;

% save results, zipped if necessary
daterun=date;
fprintf('cellxcnodb.m: SAVING to %s\n',params.outfile);
if params.zipoutfile,
   skernfile=basename(params.outfile);
   
   % then save everything (inclusive to avoid future screw-ups
   % involving accidentally not saving new important stuff)
   save([getenv('TMP'),'/',skernfile]);
   
   % compress and copy over network if necessary
   if checkbic,
      unix(['gzip ',getenv('TMP'),'/',skernfile]);
      if checkbic==1,
         unix(['jsend ',getenv('TMP'),'/',skernfile,'.gz ',...
               params.outfile,'.gz']);
      elseif checkbic==2,
         unix(['bsend ',getenv('TMP'),'/',skernfile,'.gz ',...
               params.outfile,'.gz']);
      end
      delete([getenv('TMP'),'/',skernfile,'.gz']);
   else
      unix(['gzip -cf ',getenv('TMP'),'/',skernfile,' > ',...
            params.outfile,'.gz']);
      delete([getenv('TMP'),'/',skernfile]);
   end
else
   save(params.outfile);
end

sres=sprintf('preddata(%d).cellid=''%s'';',rundata.id,cellid);
sres=sprintf('%s preddata(%d).predxc=%s;',...
             sres,rundata.id,mat2string(predxc));
sres=sprintf('%s preddata(%d).predp=%s;',...
             sres,rundata.id,mat2string(predp));
sres=sprintf('%s preddata(%d).prederr=%s;',...
             sres,rundata.id,mat2string(prederr));
sres=sprintf('%s preddata(%d).predmse=%s;',...
             sres,rundata.id,mat2string([]));
sres=sprintf('%s preddata(%d).predfix=%s;',...
             sres,rundata.id,mat2string(predfix));
sres=sprintf('%s preddata(%d).predfixerr=%s;',...
             sres,rundata.id,mat2string(predfixerr));
sres=sprintf('%s preddata(%d).sigfit=%s;',...
             sres,rundata.id,mat2string(modeloptridx));
sres=sprintf('%s preddata(%d).sfsfit=%s;',...
             sres,rundata.id,mat2string(modeloptsfs));
sres=sprintf('%s preddata(%d).expxc=%s;',...
             sres,rundata.id,mat2string([]));

% determine whether there's already a results record and delete
% it so that it can be replaced.  maybe this should be
% superceded by an archiving scheme someday?
sql=['SELECT * FROM sResults WHERE runid=',num2str(rundata.id)];
resdata=mysql(sql);
if length(resdata)>0,
   mysql(['DELETE FROM sResults WHERE id=',num2str(resdata(1).id)]);
end

sqlinsert('sResults',...
          'runid',rundata.id,...
          'batch',rundata.batch,...
          'matstr',sres);

% record that we're done with this queue entry
dbsetqueue(BATQUEUEID,1000,1);


return


