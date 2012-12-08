% function xcstsep(cellid,batch)
%
% load results from resfile and do separable time fit
%
% r=0 if no entries found in db, =1 otherwise
%
function xcstsep(cellid,batch)

disp('xcstsep.m: STARTING');

global BATQUEUEID
dbopen;

sql=['SELECT * from sRunData WHERE cellid="',cellid,'"',...
     ' AND batch=',num2str(batch)];
superrundata=mysql(sql);

sql=['SELECT * from sBatch WHERE id=',num2str(batch)];
superbatchdata=mysql(sql);

sparams=superbatchdata;
eval(char(superbatchdata.parmstring));
sparams.estbatch=getparm(sparams,'estbatch',[23 29 32]);
sparams.nlmatch=getparm(sparams,'nlmatch',[1 1 1 1]);
sparams.fitbatchidx=getparm(sparams,'fitbatchidx',1);
sparams.nosupp=getparm(sparams,'nosupp',0);
sparams.predbatch=strsep(sparams.predbatch,',');
if isfield(sparams,'maxlag') & length(sparams.maxlag)>=2,
   % do nothing, this is a good format for running cellxc
else
   sparams.maxlag=[getparm(sparams,'minlag',-6) getparm(sparams,'maxlag',13)];
end

VCELLXC=3;

clear ESTIMATIONPHASE VALIDATIONPHASE
global ESTIMATIONPHASE VALIDATIONPHASE
ESTIMATIONPHASE=1;
VALIDATIONPHASE=0;

goodbatch=zeros(1,length(sparams.estbatch));
batchcount=length(sparams.estbatch);
for ii=1:batchcount,
   sql=['SELECT * from sRunData WHERE cellid="',cellid,'"',...
           ' AND batch=',num2str(sparams.estbatch(ii))];
   trd=mysql(sql);
   if ~isempty(trd),
      rundata(ii)=trd;
      fprintf('loading %s\n',trd.resfile);
      z{ii}=zload([rundata(ii).respath,rundata(ii).resfile,'.gz']);
      
      %keyboard
      
      if isfield(z{ii},'goodbatchrange'),
         nlidx=find(sparams.nlmatch(ii)==z{ii}.goodbatchrange);
      elseif size(z{ii}.strf,1)>=sparams.nlmatch(ii);
         nlidx=sparams.nlmatch(ii);
      else
         nlidx=[];
      end
      
      if ~isempty(nlidx),
         goodbatch(ii)=1;
      end
   end
end

if ~goodbatch(sparams.fitbatchidx),
   disp('no fitbatchidx entry found in db!');
   return
end

% set parameters for rc and file organization
params=z{sparams.fitbatchidx}.params;
params.cellid=cellid;
times=z{sparams.fitbatchidx}.times;
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

params.maxlag=sparams.maxlag;
params.sfscount=8;
%params.sffiltsigma=5;
params.sffiltsigma=[1.0 1.2 1.5 1.8];
params.nloutparm=2;
params.nlidxsave=2;
params.meansub=1;
params.fitboot=0;
params.sharpspacenorm=0;
params.smoothtime=0;
%params.stimfilterparms={0,0,1,1,20};

% load the stimulus & response data (from the first--ie,review--batch)
xcloadfiles;
stimbak=stim;
movlen=size(stim,1);

goodbatchrange=find(goodbatch);

if ~isfield(sparams,'fitboot'),
   bootcount=1;
else
   bootcount=sparams.fitboot;
end

% do rc in time for spatial kernels from each estimation stim class
for bootidx=1:bootcount,
   for bidx=1:length(goodbatchrange),
      fprintf('batch %d, nlmatch=%d, boot %d/%d:\n',...
              sparams.estbatch(goodbatchrange(bidx)),...
              sparams.nlmatch(goodbatchrange(bidx)),bootidx,bootcount);
      
      % pull out spatial resposne function by summing over kernel...
      
      % use specified spatial kernel
      if isfield(z{goodbatchrange(bidx)},'goodbatchrange'),
         nlidx=find(sparams.nlmatch(goodbatchrange(bidx))==...
                    z{goodbatchrange(bidx)}.goodbatchrange);
      else
         nlidx=sparams.nlmatch(goodbatchrange(bidx));
      end
      mstimbak=zeros(size(stimbak,2),1);
      if size(z{goodbatchrange(bidx)}.strf,2)>=bootidx,
         hfull=z{goodbatchrange(bidx)}.strf(nlidx,bootidx).h(:,1:end);
         hpub=z{goodbatchrange(bidx)}.strf(nlidx,bootidx).powunbiased;
         if params.meansub,
            mstimbak=z{goodbatchrange(bidx)}.strf(nlidx,bootidx).mS;
         end
      else
         hfull=z{goodbatchrange(bidx)}.strf(nlidx).h(:,1:end);
         hpub=z{goodbatchrange(bidx)}.strf(nlidx).powunbiased;
         if params.meansub,
            mstimbak=z{goodbatchrange(bidx)}.strf(nlidx).mS;
         end
      end
      if isfield(z{goodbatchrange(bidx)},'mSA2'),
         mSA2save=z{goodbatchrange(bidx)}.mSA2;
      else
         mSA2save=z{goodbatchrange(bidx)}.strf(nlidx).sSA2;
      end
      
      hsum=sum(hfull,2);
      htime=sum(hfull,1);
      hsump=sum(hfull.*(hfull>0),2);
      
      % spatial kernel is first pc of the full space-time kernel
      [u,s,v]=svd(hfull);
      
      % make sure the sign is correct... for display, really. 
      % (well, also for some of the adaptation index measurements too)      
      if (htime*v(:,1) > 0 & hsum'*u(:,1) > 0) | hsump'*u(:,1) > 0,
         hspace=u(:,1);
         tempresp0=v(:,1)'.*s(1);
      else
         hspace=-u(:,1);
         tempresp0=-v(:,1)'.*s(1);
      end
      
      eigrat=s(1).^2./sum(diag(s).^2);
      fprintf('eigrat=%.2f\n',eigrat);
      
      % remove negative coefficient values from spatial kernel
      if sparams.nosupp,
         
         %if isfield(z{goodbatchrange(bidx)}.strf(nlidx,bootidx),...
         %           'hspacebiased');
         %   hspace=z{goodbatchrange(bidx)}.strf(nlidx,bootidx).hspacebiased;
         %end
         
         if sparams.nosupp==1,
            
            disp('no suppression allowed!');
            hspace=hspace.*(hspace>0);
         elseif sparams.nosupp==2,
            if bidx>1,
               hspace2=savestrf(1).hspace;
               hspace=hspace.*(hspace>0) + hspace2.*(hspace2<0);
            end
         elseif sparams.nosupp==-1,
            disp('no un-suppression allowed!');
            hspace=hspace.*(hspace<0);
         end
      end
      
      if var(hspace)>0,
         
         % project stim onto spatial response function
         stim=(stimbak-repmat(mstimbak',size(stimbak,1),1))*hspace;
         
         % do rc
         firstseg=1;
         spacecount=1;
         params.maxlag=sparams.maxlag;
         xccore;
      else
         mH=zeros(1,diff(params.maxlag)+1,params.sfscount);
         eH=zeros(1,diff(params.maxlag)+1,params.sfscount);
         mSall=0;
      end
      
      % choose regularization and shrinkage for temporal kernel
      clear strf
      
      if isfield(sparams,'altfit'),
         eval(sparams.altfit);
      else
         xcfit;
      end
      
      nlidx=params.nlidxsave;
      %strf(nlidx).h((end-1):end)=0;  %zero out end
      strf(nlidx).tempresp=strf(nlidx).h;
      strf(nlidx).tempresp0=tempresp0;
      strf(nlidx).hspace=hspace;
      strf(nlidx).h=hspace*strf(nlidx).h;
      strf(nlidx).mS=mstimbak;
      strf(nlidx).powunbiased=hpub;
      strf(nlidx).mH=mH;
      strf(nlidx).eH=eH;
      strf(nlidx).sSA2=mSA2save;
      
      %keyboard
      
      if 1,
         % re-measure threshold -- this may help if the last couple
         % bins of the temporal response were non-zero -- a la the
         % zeroing out that happened after xcfit
         
         % scale linear kernel to match stim space -- this isn't
         % really necessary (the threshold is) but why not be tidy?
         linpred=kernpredict(strf(nlidx).h,...
                    (stimbak'-repmat(mstimbak,1,size(stimbak,1))),1,0);
         tgoodidx=find(~isnan(resp));
         r0=resp(tgoodidx)-mean(resp(tgoodidx));
         r1=linpred(tgoodidx)-mean(linpred(tgoodidx));
         d1=sum(r1.^2);
         if d1>0,
            scf=sum(r0.*r1)./d1;
         else
            scf=1;
         end
         
         % adjust kernel by scaling factor
         strf(nlidx).h=strf(nlidx).h .* scf;
         strf(nlidx).hspace=strf(nlidx).hspace .* scf;
         
         % find optimal threshold, given the new scaling
         linpred=kernpredict(strf(nlidx).h,...
                   (stimbak'-repmat(mstimbak,1,size(stimbak,1))),1,0);
         strf(nlidx).nltype='post-sum thresh';
         strf(nlidx).nlparms=findthresh(linpred(tgoodidx),resp(tgoodidx),0);
      end
      
      savestrf(bidx,bootidx)=strf(nlidx);
      
      bexpxc(bidx,bootidx)=expxc(nlidx);
   end
end


%
% now use optimal temporal response to get optimal spatial kernel!
%
fprintf('Finding optimal spatial kernel, given optimal timecourse:\n');
respbak=resp;
bootstep=size(resp,1)./bootcount;
for bootidx=1:bootcount,
   refitidx=find(goodbatchrange==sparams.fitbatchidx);
   htime=savestrf(refitidx,bootidx).tempresp;
   spacecount=size(stimbak,2);
   th=ones(spacecount,1)*htime;
   
   if params.meansub,
      mstimbak=mean(stimbak,1)';
   else
      mstimbak=zeros(size(stimbak,2),1);
   end
   
   stim=kernpredict(th,stimbak'-repmat(mstimbak,[1 size(stimbak,1)]),...
                    spacecount,0,1);
   
   if bootcount>1,
      fitstart=round((bootidx-1)*bootstep+1);
      fitstop=round(bootidx*bootstep);
      fprintf('BOOT %d/%d: fit-excl: %d-%d  ',...
              bootidx,bootcount,fitstart,fitstop);
      
      % reset resp and nan out the ranges reserved for fit and pred
      resp=respbak;
      resp(fitstart:fitstop,:)=nan;
   end
   
   % do ic
   firstseg=1;
   params.maxlag=[0 0];
   params.sfscount=sparams.sfscount;
   params.sfsstep=sparams.sfsstep;
   params.sffiltsigma=sparams.sffiltsigma;
   params.resampcount=sparams.resampcount;
   if length(params.sffiltsigma)==1,
      fprintf(' sfs: %d/%.1f sffilt: %d\n',...
              params.sfscount,params.sfsstep,params.sffiltsigma);
   else
      fprintf(' sfs: %d/%.1f sffilt: %.1f-%.1f (%d)\n',...
              params.sfscount,params.sfsstep,params.sffiltsigma(1),...
              params.sffiltsigma(end),length(params.sffiltsigma));
   end
   
   xccore;
   
   %keyboard
   
   % choose regularization and shrinkage for temporal kernel
   clear tstrf strf fdata
   if isfield(sparams,'altfit'),
      eval(sparams.altfit);
   else
      xcfit;
   end

   REDOTIME=0;
   if REDOTIME,
      % now go back and optimally fit time to new spatial kernel!
      nlidx=params.nlidxsave;
      hspace=strf(nlidx).h;
      
      % project stim onto spatial response function
      stim=stimbak*hspace;
      
      % do rc
      firstseg=1;
      spacecount=1;
      params.maxlag=sparams.maxlag;
      params.sfscount=5;
      xccore;
      
      % choose regularization and shrinkage for temporal kernel
      clear strf
      if isfield(sparams,'altfit'),
         eval(sparams.altfit);
      else
         xcfit;
      end
      
      nlidx=params.nlidxsave;
      strf(nlidx).tempresp=strf(nlidx).h;
      strf(nlidx).tempresp0=htime;
      strf(nlidx).hspace=hspace;
      strf(nlidx).h=hspace*strf(nlidx).h;
      strf(nlidx).mH=mH;
      strf(nlidx).eH=eH;
      strf(nlidx).sSA2=mSA2;
   else
      % skip re-fitting temporal kernel ... seems to mess things up
      nlidx=params.nlidxsave;
      strf(nlidx).tempresp=htime;
      strf(nlidx).tempresp0=htime;
      strf(nlidx).hspace=strf(nlidx).h;
      strf(nlidx).h=strf(nlidx).h*htime;
      strf(nlidx).mS=mstimbak;
      sampidx=round(linspace(1,params.sfscount,10));
      strf(nlidx).mH=mH(:,:,sampidx);
      strf(nlidx).eH=eH(:,:,sampidx);
      strf(nlidx).sSA2=mSA2;
   end
   
   if 1,
      % scale linear kernel to match stim space
      linpred=kernpredict(strf(nlidx).h,...
                          (stimbak'-repmat(mstimbak,1,size(stimbak,1))),1,0);
      tgoodidx=find(~isnan(resp));
      r0=resp(tgoodidx)-mean(resp(tgoodidx));
      r1=linpred(tgoodidx)-mean(linpred(tgoodidx));
      d1=sum(r1.^2);
      if d1>0,
         scf=sum(r0.*r1)./d1;
      else
         scf=1;
      end
      
      % adjust kernel by scaling factor
      strf(nlidx).h=strf(nlidx).h .* scf;
      strf(nlidx).hspace=strf(nlidx).hspace .* scf;
      
      % find optimal threshold, given the new scaling
      linpred=kernpredict(strf(nlidx).h,...
                          (stimbak'-repmat(mstimbak,1,size(stimbak,1))),1,0);
      strf(nlidx).nltype='post-sum thresh';
      strf(nlidx).nlparms=findthresh(linpred(tgoodidx),resp(tgoodidx),0);
   end

   savestrf(length(find(goodbatch))+1,bootidx)=strf(nlidx);
end

goodbatch=cat(2,goodbatch,1);
goodbatchrange=find(goodbatch);
sparams.estbatch=cat(1,sparams.estbatch(:),...
                      sparams.estbatch(sparams.fitbatchidx));
batchcount=batchcount+1;

disp('generating bias-matched spatial kernels...');
mSA2=savestrf(end).sSA2;

for bidx=1:length(savestrf),
   if sparams.estbatch(goodbatchrange(bidx))==sparams.estbatch(end),
      
      % ie, already got the corr bias in it (rev)
      savestrf(bidx).hspacebiased=savestrf(bidx).hspace;
   else
      % extract spatially biased kernel, fully decorrelate then
      % recorrelate according to spatial stats of last kernel and
      % decorrelate that as much as the last kernel was
      % decorrelated
      
      %tt=z{goodbatchrange(bidx)}.mSA2;
      %tt=pinv(tt,0.000001);
      %hspace=z{goodbatchrange(bidx)}.mH(:,1,1);
      %hspace=mSA2*tt*hspace;
      
      hspace=mSA2*savestrf(bidx).hspace;
      
      th=normalizereg(hspace,mSA2,[],params.sfscount,params.sfsstep);
      savestrf(bidx).hspacebiased=th(:,:,savestrf(end).parms.sfsfit);
   end
end

% do rc in time for biased spatial kernels
optspacebidx=length(savestrf);

% convert to format for validation
strf=savestrf;
strfcount=length(strf);

% check how correlated different linear preds are
mod_psth=zeros(size(stimbak,1),strfcount,4);
for bidx=1:strfcount,
   hh=strf(bidx).hspace;
   tt=strf(bidx).tempresp;
   hsb=strf(bidx).hspacebiased;
   
   mod_psth(:,bidx,1)=kernpredict(hh*tt,stimbak',1,0);
   mod_psth(:,bidx,2)=kernpredict((hh.*(hh>0))*tt,stimbak',1,0);
   mod_psth(:,bidx,3)=kernpredict((hh.*(hh<0))*tt,stimbak',1,0);
   mod_psth(:,bidx,4)=kernpredict(hsb*tt,stimbak',1,0);
end

estxc=repmat(eye(strfcount),[1 1 4]);
biasxc=zeros(strfcount,3);
for b1=1:strfcount,
   for b2=b1+1:strfcount,
      for otheridx=1:4,
         ggidx=find(~isnan(mod_psth(:,b1,otheridx)) & ...
                    ~isnan(mod_psth(:,b2,otheridx)));
         if var(mod_psth(ggidx,b1,otheridx))>0 & ...
               var(mod_psth(ggidx,b2,otheridx))>0,
            estxc(b1,b2,otheridx)=xcov(mod_psth(ggidx,b1,otheridx),...
                                       mod_psth(ggidx,b2,otheridx),0,'coeff');
            estxc(b2,b1,otheridx)=estxc(b1,b2,otheridx);
         end
      end
   end
   for otheridx=1:3,
      if var(mod_psth(ggidx,b1,1))>0 & ...
            var(mod_psth(ggidx,b1,otheridx+1))>0,
         biasxc(b1,otheridx)=xcov(mod_psth(ggidx,b1,1),...
                                  mod_psth(ggidx,b1,otheridx+1),0,'coeff');
      end
   end
end

estxc(:,:)
biasxc

clf
linestr={'b-','r-','g-','b--','r--','g--'};
for bidx=1:length(find(goodbatch)),
   plot(0:14:(length(strf(bidx).tempresp)-1)*14,...
        strf(bidx).tempresp,linestr{goodbatchrange(bidx)});
   hold on
   %plot(0:14:(length(strf(bidx).tempresp)-1)*14,...
   %     strf(bidx).tempresp0,[linestr{bidx},'-']);
end
hold off
drawnow

%%
%% VALIDATION SECTION - PREDICT NOVEL RESPONSES
%%

ESTIMATIONPHASE=0;
VALIDATIONPHASE=1;

predbatchcount=length(sparams.predbatch);
if predbatchcount>0,
   
   % multiple validation run classes
   for predidx=1:predbatchcount,
      
      fprintf('Predicting responses to batch %d.\n',...
              sparams.predbatch{predidx});
      
      % figure out pred files for the current batch
      [pcellfiledata,ptimes,pbatchdata]=...
          cellfiletimes(params.cellid,params.predbatch{predidx});
      
      % does this cell have data for batchid=predidx?
      if length(pcellfiledata)>0,
         predparams=params;
         predparams.stimfiles={};
         predparams.respfiles={};
         predparams.stimcrfs=[];
         for ii=1:length(pcellfiledata),
            predparams.stimfiles{ii}=[pcellfiledata(ii).stimpath,...
                    pcellfiledata(ii).stimfile];
            predparams.respfiles{ii}=[pcellfiledata(ii).path,...
                    pcellfiledata(ii).respfile];
            if pcellfiledata(ii).stimfilecrf>0,
               predparams.stimcrfs(ii)=pcellfiledata(ii).stimfilecrf;
            else
               tstimpix=strsep(pcellfiledata(ii).stimiconside,',');
               if length(tstimpix)>0,
                  tstimpix=tstimpix{1};
               end
               sql=['SELECT * FROM gCellMaster WHERE cellid="',...
                    params.cellid,'"'];
               celldata=mysql(sql);
               predparams.stimcrfs(ii)=tstimpix./celldata.rfsize;
            end
         end
         
         tpredstartframe=ptimes(3).start;
         tpredstopframe=ptimes(3).stop;
         tpredfile=ptimes(3).fileidx;
         
         predparams.stimloadcmd=pbatchdata.stimloadcmd;
         predparams.stimloadparms=strsep(pbatchdata.stimloadparms,',');
         predparams.stimfiltercmd=pbatchdata.stimfiltercmd;
         predparams.stimfilterparms=strsep(pbatchdata.stimfilterparms,',');
         predparams.resploadcmd=pbatchdata.resploadcmd;
         predparams.resploadparms=strsep(pbatchdata.resploadparms,',');
         predparams.respfiltercmd=pbatchdata.respfiltercmd;
         predparams.respfilterparms=strsep(pbatchdata.respfilterparms,',');
         predparams.times=ptimes;
         
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
   
elseif params.predfrac>0,
   
   if isempty(BATQUEUEID),
      %keyboard
   end
   
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
gr=repmat(goodbatch(:),[1 bootcount]);
gr=find(gr);

predxc=nan.*zeros(batchcount*bootcount,1,1,predbatchcount);
predxc(gr,:,:,:)=cat(4,predres.predxc);
predinf=nan.*zeros(batchcount*bootcount,1,1,predbatchcount);
predinf(gr,:,:,:)=cat(4,predres.predinf);
predone=nan.*zeros(batchcount*bootcount,1,1,predbatchcount);
predone(gr,:,:,:)=cat(4,predres.predone);
predp=nan.*zeros(batchcount*bootcount,1,1,predbatchcount);
predp(gr,:,:,:)=cat(4,predres.predp);
prederr=nan.*zeros(batchcount*bootcount,1,1,predbatchcount);
prederr(gr,:,:,:)=cat(4,predres.prederr);
predfix=nan.*zeros(batchcount*bootcount,1,1,predbatchcount);
predfix(gr,:,:,:)=cat(4,predres.predfix);
predfixerr=nan.*zeros(batchcount*bootcount,1,1,predbatchcount);
predfixerr(gr,:,:,:)=cat(4,predres.predfixerr);
predmse=nan.*zeros(batchcount*bootcount,1,1,predbatchcount);
predmse(gr,:,:,:)=cat(4,predres.predmse);

if params.predfrac>0 & ~params.fitboot,
   
   predxc=permute(reshape(predxc,batchcount,respcount,attcount,...
                          predbatchcount),[4 1 3 2]);
   predinf=permute(reshape(predinf,batchcount,respcount,attcount,...
                          predbatchcount),[4 1 3 2]);
   predone=permute(reshape(predone,batchcount,respcount,attcount,...
                          predbatchcount),[4 1 3 2]);
   predp=permute(reshape(predp,batchcount,respcount,attcount,...
                         predbatchcount),[4 1 3 2]);
   prederr=permute(reshape(prederr,batchcount,respcount,...
                           attcount,predbatchcount),[4 1 3 2]);
   predmse=permute(reshape(predmse,batchcount,respcount,...
                           attcount,predbatchcount),[4 1 3 2]);
   predfix=permute(reshape(predfix,batchcount,respcount,...
                           attcount,predbatchcount),[4 1 3 2]);
   predfixerr=permute(reshape(predfixerr,batchcount,respcount,...
                              attcount,predbatchcount),[4 1 3 2]);
   
   expxc=bexpxc;
elseif params.predfrac>0,
   
   % predres ( batch x 1) . predxc ( nlcount*bootcount*respcount x 1)
   
   % vpredxc should be batch X nl X boot X resp
   % predxc should be batch X nl X 1 X resp   (averaged over
   % boots!)
   
   vpredxc=permute(reshape(predxc,batchcount,attcount,respcount,...
                           predbatchcount),[4 1 2 3]);
   vpredinf=permute(reshape(predinf,batchcount,attcount,respcount,...
                           predbatchcount),[4 1 2 3]);
   vpredone=permute(reshape(predone,batchcount,attcount,respcount,...
                           predbatchcount),[4 1 2 3]);
   vpredp=permute(reshape(predp,batchcount,attcount,respcount,...
                          predbatchcount),[4 1 2 3]);
   vprederr=permute(reshape(prederr,batchcount,attcount,...
                            respcount,predbatchcount),[4 1 2 3]);
   vpredmse=permute(reshape(predmse,batchcount,attcount,...
                            respcount,predbatchcount),[4 1 2 3]);
   vpredfix=permute(reshape(predfix,batchcount,attcount,...
                            respcount,predbatchcount),[4 1 2 3]);
   vpredfixerr=permute(reshape(predfixerr,batchcount,attcount,...
                               respcount,predbatchcount),[4 1 2 3]);
   
   predxc=mean(vpredxc,3);
   predinf=mean(vpredinf,3);
   predone=mean(vpredone,3);
   predp=mean(vpredp,3);
   prederr=mean(vprederr,3);
   predmse=mean(vpredmse,3);
   predfix=mean(vpredfix,3);
   predfixerr=mean(vpredfixerr,3);
   
   vexpxc=bexpxc;
   expxc=mean(bexpxc,4);
   
else
   predxc=ones([1 size(strf)]).*nan;
   predinf=ones([1 size(strf)]).*nan;
   predone=ones([1 size(strf)]).*nan;
   predp=zeros([1 size(strf)]);
   prederr=ones([1 size(strf)]).*nan;
   predmse=ones([1 size(strf)]).*nan;
   predfix=ones([1 size(strf)]).*nan;
   predfixerr=ones([1 size(strf)]).*nan;
end

for estidx=1:size(predxc,2);
   fprintf('batch %d est: valxc=',sparams.estbatch(estidx));
   for bidx=1:size(predxc,1),
      fprintf('%7.3f',predxc(bidx,estidx));
   end
   fprintf('\n');
end

% save results, zipped if necessary
daterun=date;

params.zipoutfile=1;
params.outfile=[superrundata.respath,superrundata.resfile];

clear z stimbak cdata
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

sres=sprintf('preddata(%d).cellid=''%s'';',superrundata.id,cellid);
sres=sprintf('%s preddata(%d).predxc=%s;',...
             sres,superrundata.id,mat2string(predxc));
sres=sprintf('%s preddata(%d).predp=%s;',...
             sres,superrundata.id,mat2string(predp));
sres=sprintf('%s preddata(%d).prederr=%s;',...
             sres,superrundata.id,mat2string(prederr));
sres=sprintf('%s preddata(%d).predinf=%s;',...
             sres,superrundata.id,mat2string(predinf));
sres=sprintf('%s preddata(%d).predfix=%s;',...
             sres,superrundata.id,mat2string(predfix));
sres=sprintf('%s preddata(%d).predfixerr=%s;',...
             sres,superrundata.id,mat2string(predfixerr));
sres=sprintf('%s preddata(%d).sigfit=%s;',...
             sres,superrundata.id,mat2string(sigfit));
sres=sprintf('%s preddata(%d).sfsfit=%s;',...
             sres,superrundata.id,mat2string(sfsfit));
sres=sprintf('%s preddata(%d).expxc=%s;',...
             sres,superrundata.id,mat2string(expxc));

% determine whether there's already a results record and delete
% it so that it can be replaced.  maybe this should be
% superceded by an archiving scheme someday?
sql=['SELECT * FROM sResults WHERE runid=',num2str(superrundata.id)];
resdata=mysql(sql);
if length(resdata)>0,
   mysql(['DELETE FROM sResults WHERE id=',num2str(resdata(1).id)]);
end

sqlinsert('sResults',...
          'runid',superrundata.id,...
          'batch',superrundata.batch,...
          'matstr',sres);


disp('xcstsep.m: DONE');


