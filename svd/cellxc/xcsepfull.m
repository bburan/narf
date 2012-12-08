% function xcsepfull(cellid,batch)
%
% space-time separable kernel esimation, interchangeable with
% cellxcmaster, which does full kernel estimation.
%
% plus features for doing jackknifed kernels and measuring the
% estimation noise ceiling on predictions (bfracs, bootcount)
% NOTE: bfracs should be matched to bootcount if doing
% jackknifing. (check details of sorting through bfrac in code)
%
% created SVD 1/04  ripped off xcstsep
%
function xcsepfull(cellid,batch)

disp('xcsepfull.m: STARTING');

global BATQUEUEID
dbopen;

sql=['SELECT * from sRunData WHERE cellid="',cellid,'"',...
     ' AND batch=',num2str(batch)];
rundata=mysql(sql);

if length(rundata)==0,
   disp('no db entry!');
   return
end

%%
%% THIS SECTION HACKED FROM CELLXCMASTER
%%

% figure out what files to use for what stage of the analysis
[cellfiledata,times,params]=cellfiletimes(cellid,rundata.batch);

% fill up params structure for passing to cellxcnodb
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

if isfield(params,'maxlag') & length(params.maxlag)>=2,
   % do nothing, this is a good format for running cellxc
else
   params.maxlag=[getparm(params,'minlag',-6) getparm(params,'maxlag',13)];
end

params.docellfit2=0;
params.shrinkage=1;   % rather than just thresholding
params.repexclude=0;
params.fitboot=0;
params.smoothtime=0;
params.tbinms=getparm(params,'tbinms',16);
params.nrandxcov=getparm(params,'nrandxcov',200);
params.sharpspacenorm=getparm(params,'sharpspacenorm',0);

params.zipoutfile=1;
params.outfile=[rundata.respath,rundata.resfile];

%% Maybe want to change these at some point?
params.nosupp=0;
params.nloutparm=2;
params.nlidxsave=2;
params.meansub=getparm(params,'meansub',1);

if params.predfrac>0,
   params.bfracs=[0.05 0.1 0.25 0.60 0.00];
else
   params.bfracs=[0.05 0.1 0.25 0.55 0.05];
end

if length(params.parmstring)>0,
   eval(char(params.parmstring));
end

% call xcfilefracs to get times
if params.fitfrac>0 & params.fitboot>1,
   tparams=params;
   tparams.fitfrac=0;
   tparams.predfrac=0;
   times=xcfilefracs(tparams);
end

%%
%% NEXT SECTION HACKED FROM CELLXCNODB
%%

VCELLXC=3;

clear ESTIMATIONPHASE VALIDATIONPHASE
global ESTIMATIONPHASE VALIDATIONPHASE
ESTIMATIONPHASE=1;
VALIDATIONPHASE=0;

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

xcloadfiles;

if params.fitboot,
   
   resp0=resp;
   stim0=stim;
   
   filelens=stimstop-stimstart+1;
   filecount=max(find(filelens>1));
   filelens=filelens(1:filecount);
   
   totlen=sum(filelens);
   clen=[0; cumsum(filelens)];
   
   filegoodidx={};
   filegoodlens=zeros(filecount,1);
   for fidx=1:filecount,
      filegoodidx{fidx}=find(~isnan(resp(clen(fidx)+(1:filelens(fidx)))));
      filegoodlens(fidx)=length(filegoodidx{fidx});
   end
   bootcount=params.fitboot;
   bootstep=filegoodlens./bootcount;
   bcount=length(params.bfracs);
   predmtx=nan.*ones(length(resp),bcount-1);
   
   for bootidx=1:bootcount,
      
      % segregate out the exploratory, fit and validation data
      % in order to keep things balanced, take the same fraction
      % from each file included in the stim and resp matrices
      %if params.fitfrac>0,
      %   fdata.stim=[];
      %   fdata.resp=[];
      %end
      
      % if predfrac>0 use the same validation data for each
      % bootstrap. this is done in the standard section far
      % below. otherwise create anew validation set during each boot
      % iteration   
      cdata.stim=[];
      cdata.resp=[];
      
      resp=resp0;
      
      bresp=nan.*ones(size(resp,1),bcount);
      bstart=zeros(filecount,bcount);
      bstop=zeros(filecount,bcount);
      
      for fidx=1:filecount,
         
         bbase=(bootidx-1)*bootstep(fidx);
         for segidx=1:bcount,
            bstart(fidx,segidx)=round(bbase+1);
            bstop(fidx,segidx)=round(bbase+filegoodlens(fidx)*...
                                     params.bfracs(segidx));
            if bstart(fidx,segidx)>filegoodlens(fidx),
               bstart(fidx,segidx)=bstart(fidx,segidx)-filegoodlens(fidx);
               bstop(fidx,segidx)=bstop(fidx,segidx)-filegoodlens(fidx);
            end
            
            if bstop(fidx,segidx)>filegoodlens(fidx),
               bstop(fidx,segidx)=bstop(fidx,segidx)-filegoodlens(fidx);
               bbase=bstop(fidx,segidx);
               
               bstart(fidx,segidx)=clen(fidx)+...
                   filegoodidx{fidx}(bstart(fidx,segidx));
               bstop(fidx,segidx)=clen(fidx)+...
                   filegoodidx{fidx}(bstop(fidx,segidx));
               
               bresp(bstart(fidx,segidx):clen(fidx+1),segidx)=...
                   resp0(bstart(fidx,segidx):clen(fidx+1),1);
               bresp(clen(fidx)+1:bstop(fidx,segidx),segidx)=...
                   resp0(clen(fidx)+1:bstop(fidx,segidx),1);
            else
               bbase=bstop(fidx,segidx);
               if bstop(fidx,segidx)==0,
                  bstop(fidx,segidx)=1;
               end
               
               bstart(fidx,segidx)=clen(fidx)+...
                   filegoodidx{fidx}(bstart(fidx,segidx));
               bstop(fidx,segidx)=clen(fidx)+...
                   filegoodidx{fidx}(bstop(fidx,segidx));
               
               bresp(bstart(fidx,segidx):bstop(fidx,segidx),segidx)=...
                   resp0(bstart(fidx,segidx):bstop(fidx,segidx),1);
            end
         end
      end
      
      if params.predfrac==0,
         valmask=...
             union(find(~isnan(shift(bresp(:,bcount),-params.maxlag(1)))),...
                   find(~isnan(shift(bresp(:,bcount),-params.maxlag(2)))));
         cdata.stim=stim0(valmask,:);
         cdata.resp=bresp(valmask,bcount);
      end
      
      for segidx=1:bcount-1,  % last bcount is validation data!
         
         fprintf('\nbootidx=%d segidx=%d (%.0f%%):\n',...
                 bootidx,segidx,params.bfracs(segidx)*100);
         
         resp=bresp(:,segidx);
         
         
         fitmask=union(find(~isnan(shift(resp,-params.maxlag(1)))),...
                       find(~isnan(shift(resp,-params.maxlag(2)))));
         resp=resp(fitmask);
         
         % make sure there aren't any early stray nans that snuck
         % in ... seems to be a problem when there are a bunch of
         % dropped frames.
         resp(1:(params.maxlag(2)+1),:)=nan;
         resp((end+params.maxlag(1)-1):end,:)=nan;
         
         stim=stim0(fitmask,:);
         
         firstseg=1;
         if diff(params.maxlag)>0 & size(stim,2)>1,
            xcsepcore;
         else
            xccore;
            xcfit;
         end
         
         % save important stuff from this bootidx and test preds
         if bootidx==1 & segidx==1,
            vexpxc=expxc(:);
            vstrf=strf;
            if params.predfrac==0,
               predres=xcval(strf,params,cdata);
            end
         else
            if params.predfrac==0,
               predres(bootidx,segidx)=xcval(strf,params,cdata);
            end
            if bootidx>1 & bootidx<bootcount & isfield(strf,'mH'),
               for ii=1:length(strf(:))-1,
                  strf(ii).mH=[];
                  strf(ii).eH=[];
                  strf(ii).sSA2=[];
               end
            end
            %strf
            vexpxc(:,bootidx,segidx)=expxc(:);
            vstrf(:,bootidx,segidx)=strf;
         end
         
         if params.predfrac==0,
            predmtx(valmask,segidx)=...
                predres(bootidx,segidx).mod_psth{1}(:,end);
         end
      end
      
   end
   
   stim=stim0;
   resp=resp0;
   clear stim0 resp0 bresp
   
   strf=vstrf;
   clear vstrf
   
else
   %
   % NON-BOOTSTRAPPED KERNEL ESTIMATION SECTION
   %
   xcsepcore;
end

if diff(params.maxlag)>0 & size(stim,2)>1,
   params.nlidxsave=3;
end

%%
%% CREATE DUMMY DIRECT MODEL CELL STRF
%%

sepfiles=strsep(params.respfiles{end},'+');
if length(sepfiles)>1,
   rr=load(sepfiles{1});
else
   rr=load(params.respfiles{end});
end

if strcmp(cellid(1:4),'mode') & strcmp(params.stimfiltercmd,'') & ...
      isstruct(rr) & isfield(rr,'STRF'),
   if size(strf,3)==bcount-1,
      strf(:,:,bcount)=strf(:,:,bcount-1);
      h=rr.STRF.Phase0;
      h=reshape(h(:,:,1:14),16*16,14)./256;
      strf(end,1,bcount).h=h;
      strf(end,1,bcount).mS(:)=128;
      strf(end,1,bcount).nlparms=0;
      hspace=h(:,4);
      [u,s,v]=svd(h);
      if u(:,1)'*hspace > 0,
         strf(end,1,bcount).hspace=u(:,1);
         strf(end,1,bcount).tempresp=v(:,1);
      else
         strf(end,1,bcount).hspace=u(:,1);
         strf(end,1,bcount).tempresp=v(:,1);
      end         
      
      bcount=bcount+1;
   end
end

%%
%% VALIDATION SECTION - PREDICT NOVEL RESPONSES
%%

ESTIMATIONPHASE=0;
VALIDATIONPHASE=1;

predbatchcount=length(params.predbatch);
if predbatchcount==0,
   params.predbatch{1}=params.batch;
   predbatchcount=1;
end

if predbatchcount>0 & params.predfrac>0,
   
   % multiple validation run classes
   for predidx=1:predbatchcount,
      
      fprintf('Predicting responses to batch %d.\n',...
              params.predbatch{predidx});
      
      % figure out pred files for the current batch
      [pcellfiledata,ptimes,pbatchdata]=...
          cellfiletimes(params.cellid,params.predbatch{predidx});
      
      % does this cell have data for batchid=predidx?
      if length(pcellfiledata)>0,
         predparams=params;
         predparams.stimfiles=pbatchdata.stimfiles;
         predparams.respfiles=pbatchdata.respfiles;
         predparams.stimcrfs=pbatchdata.stimcrfs;
         
         predparams.stimloadcmd=pbatchdata.stimloadcmd;
         predparams.stimloadparms=strsep(pbatchdata.stimloadparms,',');
         predparams.stimfiltercmd=pbatchdata.stimfiltercmd;
         predparams.stimfilterparms=strsep(pbatchdata.stimfilterparms,',');
         predparams.resploadcmd=pbatchdata.resploadcmd;
         predparams.resploadparms=strsep(pbatchdata.resploadparms,',');
         predparams.respfiltercmd=pbatchdata.respfiltercmd;
         predparams.respfilterparms=strsep(pbatchdata.respfilterparms,',');
         predparams.times=ptimes;
         
         tpredstartframe=ptimes(3).start;
         tpredstopframe=ptimes(3).stop;
         tpredfile=ptimes(3).fileidx;
         
         [cdata.stim,cdata.resp]=xcloadstimresp(tpredfile,tpredstartframe,...
                                                tpredstopframe,predparams);
      else
         cdata.stim=[];
         cdata.resp=[];
         predparams=params;
         if predidx>1,
            predres(predidx).predxc=[];
         end
      end
      nlcount=size(strf,1);
      attcount=size(strf,2);
      respcount=size(strf,3);
      
      if params.fitboot,
         predres(predidx)=xcval(strf(:),predparams,cdata);
      else
         predres(predidx)=xcval(reshape(permute(strf,[1 3 2]),...
                                        nlcount*respcount,attcount),...
                                predparams,cdata);
      end
   end
   
elseif params.predfrac>0 & params.fitboot<=1,
   
   % just assume current batch
   if ~isempty(params.batch),
      params.predbatch={params.batch};
   else
      params.predbatch={0};  % "not a batch"
   end
   
   nlcount=size(strf,1);
   attcount=size(strf,2);
   respcount=size(strf,3);
   
   if params.fitboot,
      predres=xcval(strf(:),params,times(3));
   else
      predres=xcval(reshape(permute(strf,[1 3 2]),nlcount*respcount,...
                            attcount),params,times(3));
   end
end

% consoldiate predres into compact matrices for saving to sResults
% in the database.

predxc=cat(4,predres.predxc);
predinf=cat(4,predres.predinf);
predone=cat(4,predres.predone);
predp=cat(4,predres.predp);
prederr=cat(4,predres.prederr);
predfix=cat(4,predres.predfix);
predfixerr=cat(4,predres.predfixerr);
predmse=cat(4,predres.predmse);

if params.predfrac>0 & ~params.fitboot,
   
   predxc=permute(reshape(predxc,strfcount,respcount,attcount,...
                          predbatchcount),[4 1 3 2]);
   predinf=permute(reshape(predinf,strfcount,respcount,attcount,...
                          predbatchcount),[4 1 3 2]);
   predone=permute(reshape(predone,strfcount,respcount,attcount,...
                          predbatchcount),[4 1 3 2]);
   predp=permute(reshape(predp,strfcount,respcount,attcount,...
                         predbatchcount),[4 1 3 2]);
   prederr=permute(reshape(prederr,strfcount,respcount,...
                           attcount,predbatchcount),[4 1 3 2]);
   predmse=permute(reshape(predmse,strfcount,respcount,...
                           attcount,predbatchcount),[4 1 3 2]);
   predfix=permute(reshape(predfix,strfcount,respcount,...
                           attcount,predbatchcount),[4 1 3 2]);
   predfixerr=permute(reshape(predfixerr,strfcount,respcount,...
                              attcount,predbatchcount),[4 1 3 2]);
   expxc=bexpxc;
else
   
   % predres ( batch x 1) . predxc ( nlcount*bootcount*respcount x 1)
   
   % vpredxc should be batch X nl X boot X resp
   % predxc should be batch X nl X 1 X resp   (averaged over
   % boots!)
   
   strfcount=size(strf,1);
   if bcount>1,
      respcount=1;
   end
   if bootcount>1,
      attcount=1;
   end
   
   vpredxc=permute(reshape(predxc,strfcount,attcount,respcount,...
                           bootcount,bcount-1),[5 1 4 2 3]);
   vpredinf=permute(reshape(predinf,strfcount,attcount,respcount,...
                           bootcount,bcount-1),[5 1 4 2 3]);
   vpredone=permute(reshape(predone,strfcount,attcount,respcount,...
                           bootcount,bcount-1),[5 1 4 2 3]);
   vpredp=permute(reshape(predp,strfcount,attcount,respcount,...
                           bootcount,bcount-1),[5 1 4 2 3]);
   vprederr=permute(reshape(prederr,strfcount,attcount,respcount,...
                           bootcount,bcount-1),[5 1 4 2 3]);
   vpredmse=permute(reshape(predmse,strfcount,attcount,respcount,...
                           bootcount,bcount-1),[5 1 4 2 3]);
   vpredfix=permute(reshape(predfix,strfcount,attcount,respcount,...
                           bootcount,bcount-1),[5 1 4 2 3]);
   vpredfixerr=permute(reshape(predfixerr,strfcount,attcount,respcount,...
                           bootcount,bcount-1),[5 1 4 2 3]);
   
   predxc=mean(vpredxc,3);
   predinf=mean(vpredinf,3);
   predone=mean(vpredone,3);
   predp=mean(vpredp,3);
   prederr=std(vpredxc,1,3).*sqrt(bootcount-1);
   predmse=mean(vpredmse,3);
   predfix=mean(vpredfix,3);
   predfixerr=std(vpredfixerr,1,3).*sqrt(bootcount-1);
   
end

disp('vexpxc is nl X boot X bfrac')
disp('(v)predxc is bfrac X nl ( X boot)')
for estidx=1:size(predxc,2);
   fprintf('batch %d est: valxc=',estidx);
   for bidx=1:size(predxc,1),
      fprintf('%7.3f %7.3f',predxc(bidx,estidx),predfix(bidx,estidx));
   end
   fprintf('\n');
end


% save results, zipped if necessary
daterun=date;

clear z cdata tfitres 
clear tstim tr th tH seplinpred rgoodidx respsave linpred tstrf
clear stim sSA2 u tsSA2 H savestrf fdata

if params.zipoutfile,
   skernfile=basename(params.outfile);
   fprintf('cellxcnodb.m: SAVING to %s\n',params.outfile);
   
   % then save everything (inclusive to avoid future screw-ups
   % involving accidentally not saving new important stuff)
   save([tempdir,skernfile]);
   
   % compress and copy over network if necessary
   unix(['gzip -cf ',tempdir,skernfile,' > ',...
         params.outfile,'.gz']);
   delete([tempdir,skernfile]);
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
sres=sprintf('%s preddata(%d).predinf=%s;',...
             sres,rundata.id,mat2string(predinf));
sres=sprintf('%s preddata(%d).predfix=%s;',...
             sres,rundata.id,mat2string(predfix));
sres=sprintf('%s preddata(%d).predfixerr=%s;',...
             sres,rundata.id,mat2string(predfixerr));
sres=sprintf('%s preddata(%d).sigfit=%s;',...
             sres,rundata.id,mat2string(sigfit));
sres=sprintf('%s preddata(%d).sfsfit=%s;',...
             sres,rundata.id,mat2string(sfsfit));
sres=sprintf('%s preddata(%d).expxc=%s;',...
             sres,rundata.id,mat2string(expxc));

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

disp('xcstsep.m: DONE');


