% strf=cellxcdataloaded(stim,resp,params);
%
function strf=cellxcdataloaded(stim,resp,params);

global ESTIMATIONPHASE VALIDATIONPHASE
ESTIMATIONPHASE=1;
VALIDATIONPHASE=0;

disp('cellxcdataloaded.m: INITIALIZING');
if ~exist('params','var'),
   params=[];
   params.respfiles={};
end

% set everything in params structure to default values if they
% haven't been specified
params.stimloadcmd=getparm(params,'stimloadcmd','loadimfile');
params.stimloadparms=getparm(params,'stimloadparms',[]);
params.stimfiltercmd=getparm(params,'stimfiltercmd','');
params.stimfilterparms=getparm(params,'stimfilterparms',[]);
params.resploadcmd=getparm(params,'resploadcmd','respload');
params.resploadparms=getparm(params,'resploadparms',{'',1,1,1});
params.respfiltercmd=getparm(params,'respfiltercmd','');
params.respfilterparms=getparm(params,'respfilterparms',[]);
params.kernfmt=getparm(params,'kernfmt','space');
if isfield(params,'maxlag') & length(params.maxlag)>=2,
   % do nothing, this is a good format for running cellxc
else
   params.maxlag=[getparm(params,'minlag',-6) getparm(params,'maxlag',13)];
end
params.resampcount=getparm(params,'resampcount',20);
params.resampfmt=getparm(params,'resampfmt',1);
params.expfrac=getparm(params,'expfrac',0);
params.fitfrac=getparm(params,'fitfrac',0);
params.predfrac=getparm(params,'predfrac',0.1);
params.decorrspace=getparm(params,'decorrspace',2);
params.decorrtime=getparm(params,'decorrtime',1);
params.sffiltsigma=getparm(params,'sffiltsigma',7);
params.sfscount=getparm(params,'sfscount',30);
params.sfsstep=getparm(params,'sfsstep',4);
params.smoothtime=getparm(params,'smoothtime',2);
params.predsmoothsigma=getparm(params,'predsmoothsigma',0);
params.predtype=getparm(params,'predtype',0);
params.stimwindowcrf=getparm(params,'stimwindowcrf',2);
params.nloutparm=getparm(params,'nloutparm',3);
params.sffiltthresh=getparm(params,'sffiltthresh',0);
params.sffiltsmooth=getparm(params,'sffiltsmooth',0);
params.boundary=getparm(params,'boundary','zero');
params.repexclude=getparm(params,'repexclude',0);
params.docellfit2=getparm(params,'docellfit2',0);
params.stimfmtcode=getparm(params,'stimfmtcode',0);
params.respfmtcode=getparm(params,'respfmtcode',0);
params.nlidxsave=getparm(params,'nlidxsave',2);
params.outfile=getparm(params,'outfile','/tmp/cellxcout.mat');
params.showres=getparm(params,'showres',0);
params.zipoutfile=getparm(params,'zipoutfile',0);
params.cellid=getparm(params,'cellid','CELL');
params.shrinkage=getparm(params,'shrinkage',1);
params.predbatch=getparm(params,'predbatch',{});
params.batch=getparm(params,'batch',[]);
params.fitboot=getparm(params,'fitboot',0);
params.meansub=getparm(params,'meansub',1);
params.tbinms=getparm(params,'tbinms',16);
params.nrandxcov=getparm(params,'nrandxcov',200);
params.sharpspacenorm=getparm(params,'sharpspacenorm',0);
params.keepneg=getparm(params,'keepneg',0);
params.cutextendedsilence=getparm(params,'cutextendedsilence',0);
params.spacetimesep=getparm(params,'spacetimesep',0);
%=getparm(params,'','');

firstseg=1;
rsize=size(resp);
spacecount=size(stim,2);
respcount=size(resp,2);
params.iconside=size(stim');

if params.repexclude,
    % do repeated exclusions
    xcitercore;
elseif isfield(params,'altcore'),
    eval(params.altcore);
else
    % simply do the XC. standard
    xccore;
end

%%
%% FIT SECTION - FIND OPTIMAL REGULARIZATION AND SHRINKAGE PARAMETERS
%%
if isfield(params,'altfit'),
    eval(params.altfit);
else
    xcfit2;
end


return


    

clear sH H ttSA tsSA2 tsSAfull tSR tmR tn tSA1
clear tSA sSA1 sSA2 sSAfull sSA0 SR
clear fdata

%%
%% VALIDATION SECTION - PREDICT NOVEL RESPONSES
%%

ESTIMATIONPHASE=0;
VALIDATIONPHASE=1;

batchcount=length(params.predbatch);

if batchcount==0 && params.fitboot>0 && params.predfrac==0,
   % do nothing, preds already done.
   % now skip down and piece into standard batchcount>0 area
   
elseif batchcount==0 && params.predfrac>0,
   
   %
   % this code usually runs for validation unless you're doing
   % something funky with multiple stimulus types
   %
   
   if isempty(BATQUEUEID),
      %keyboard
   end
   
   % just assume current batch
   if ~isempty(params.batch),
      params.predbatch={params.batch};
   else
      params.predbatch={};  % "not a batch"
   end
   
   nlcount=size(strf,1);
   attcount=size(strf,2);
   respcount=size(strf,3);
   
   error('times not defined');
   
   %if params.fitboot,
   %   predres=xcval(strf(:),params,times(3));
   %else
   %   predres=xcval(reshape(permute(strf,[1 3 2]),nlcount*respcount,...
   %                         attcount),params,times(3));
   %end
elseif batchcount>0,
   
   if params.fitboot,
      rrfull={};
      ppfull={};
      clear predres;
   end
   
   % multiple validation run classes
   for predidx=1:batchcount,
      
      fprintf('Predicting responses to batch %d.\n',params.predbatch{predidx});
      
      % figure out pred files for the current batch
      [pcellfiledata,ptimes,pbatchdata]=...
          cellfiletimes(params.cellid,params.predbatch{predidx});
      
      % does this cell have data for batchid=predidx?
      if length(pcellfiledata)>0,
         
         % predict whole response for datasets with zero predfrac
         if pbatchdata.predfrac==0,
            ptimes(3)=ptimes(1);
         end
         
         predparams=params;
         if params.batch~=params.predbatch{predidx},
            predparams.stimfiles=pbatchdata.stimfiles;
            predparams.respfiles=pbatchdata.respfiles;
            predparams.stimcrfs=pbatchdata.stimcrfs;
            
            predparams.stimloadcmd=pbatchdata.stimloadcmd;
            predparams.stimloadparms=pbatchdata.stimloadparms;
            if strcmp(predparams.stimloadcmd,'loadsiteraster'),
               predparams.stimloadparms{1}.channel=params.stimloadparms{1}.channel;
               predparams.stimloadparms{1}.unit=params.stimloadparms{1}.unit;
            end
            predparams.stimfiltercmd=pbatchdata.stimfiltercmd;
            if ~isempty(pbatchdata.stimfilterparms),
               pbatchdata.stimfilterparms=...
                   strrep(pbatchdata.stimfilterparms,'cellfiledata','pbatchdata');
            end
            predparams.stimfilterparms=strsep(pbatchdata.stimfilterparms,',');
            predparams.resploadcmd=pbatchdata.resploadcmd;
            predparams.resploadparms=pbatchdata.resploadparms;
            predparams.respfiltercmd=pbatchdata.respfiltercmd;
            predparams.respfilterparms=pbatchdata.respfilterparms;
            predparams.times=ptimes;
            
% $$$             % allow for exection of cellid-specific identifiers in strsep
% $$$             if ~isempty(pbatchdata.stimloadparms),
% $$$                pbatchdata.stimloadparms=...
% $$$                    strrep(pbatchdata.stimloadparms,'cellfiledata','pbatchdata');
% $$$             end
% $$$             predparams.stimloadparms=strsep(pbatchdata.stimloadparms,',');
% $$$             predparams.stimfiltercmd=pbatchdata.stimfiltercmd;
% $$$             if ~isempty(pbatchdata.stimfilterparms),
% $$$                pbatchdata.stimfilterparms=...
% $$$                    strrep(pbatchdata.stimfilterparms,'cellfiledata','pbatchdata');
% $$$             end
% $$$             predparams.stimfilterparms=strsep(pbatchdata.stimfilterparms,',');
% $$$             predparams.resploadcmd=pbatchdata.resploadcmd;
% $$$             predparams.resploadparms=strsep(pbatchdata.resploadparms,',');
% $$$             predparams.respfiltercmd=pbatchdata.respfiltercmd;
% $$$             predparams.times=ptimes;
% $$$             predparams.respfilterparms=strsep(pbatchdata.respfilterparms,',');
         end
         
         tpredstartframe=ptimes(3).start;
         tpredstopframe=ptimes(3).stop;
         tpredfile=ptimes(3).fileidx;
         cdata=[];         
         [cdata.stim,cdata.resp,extras]=...
             xcloadstimresp(tpredfile,tpredstartframe,...
                            tpredstopframe,predparams);
         if isfield(extras,'raster'),
            cdata.raster=extras.raster;
         end
      else
         cdata.stim=[];
         cdata.resp=[];
         predparams=params;
         %if predidx>1,
         %   predres(predidx).predxc=[];
         %end
      end
      nlcount=size(strf,1);
      attcount=size(strf,2);
      respcount=size(strf,3);
      
      if params.fitboot,
         for bootidx=1:params.fitboot,
            if length(cdata.resp)>0
               [tresp,tcdata]=xcgetbootset(predparams,bootidx,ptimes,...
                                           cdata.resp,cdata.stim,raster);
               
               predres(predidx,bootidx)=...
                   xcval(strf(:,bootidx),predparams,tcdata);
               
               if bootidx==1,
                  %keyboard
                  rrfull{predidx}=zeros(size(tresp)).*nan;
                  ppfull{predidx}=zeros(length(tresp),...
                          size(predres(predidx,bootidx).mod_psth{1},3)).*nan;
               end
               
               rrfull{predidx}(tcdata.fromidx)=...
                   predres(predidx,bootidx).act_resp{1}(tcdata.toidx);
               ppfull{predidx}(tcdata.fromidx,:)=...
                   predres(predidx,bootidx).mod_psth{1}(tcdata.toidx,:)-...
                   repmat(nanmean(predres(predidx,bootidx).mod_psth{1}(tcdata.toidx,:)),[length(tcdata.toidx),1]);
            else
               predres(predidx,bootidx)=xcval(strf(:,bootidx),predparams,cdata);
            end
         end
         attcount=size(predres(1).predxc,2); 
      else
          %keyboard
         predres(predidx)=xcval(reshape(permute(strf,[1 3 2]),...
                                        nlcount*respcount,attcount),...
                                predparams,cdata);
      end
   end
else
   predres.predxc=[];
   predres.predone=[];
   predres.predinf=[];
   predres.predp=[];
   predres.prederr=[];
   predres.predfix=[];
   predres.predfixerr=[];
   predres.predmse=[];
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

bcount=size(predxc,4);
if params.predfrac>0 & ~params.fitboot,
   predxc=permute(reshape(predxc,nlcount,respcount,size(predxc,3),...
                          bcount),[4 1 3 2]);
   predinf=permute(reshape(predinf,nlcount,respcount,size(predinf,3),...
                          bcount),[4 1 3 2]);
   predone=permute(reshape(predone,nlcount,respcount,size(predone,3),...
                          bcount),[4 1 3 2]);
   predp=permute(reshape(predp,nlcount,respcount,size(predp,3),...
                         bcount),[4 1 3 2]);
   prederr=permute(reshape(prederr,nlcount,respcount,size(prederr,3),...
                           bcount),[4 1 3 2]);
   predmse=permute(reshape(predmse,nlcount,respcount,...
                           size(predmse,3),bcount),[4 1 3 2]);
   predfix=permute(reshape(predfix,nlcount,respcount,...
                           size(predfix,3),bcount),[4 1 3 2]);
   predfixerr=permute(reshape(predfixerr,nlcount,respcount,...
                              size(predfixerr,3),bcount),[4 1 3 2]);
   
elseif params.predfrac>0,
   
   % predres ( batch x 1) . predxc ( nlcount*bootcount*respcount x 1)
   
   % vpredxc should be batch X nl X boot X resp
   % predxc should be batch X nl X 1 X resp   (averaged over boots!)
   vpredxc=permute(reshape(predxc,nlcount,attcount,respcount,...
                           bcount),[4 1 2 3]);
   vpredinf=permute(reshape(predinf,nlcount,attcount,respcount,...
                           bcount),[4 1 2 3]);
   vpredone=permute(reshape(predone,nlcount,attcount,respcount,...
                           bcount),[4 1 2 3]);
   vpredp=permute(reshape(predp,nlcount,attcount,respcount,...
                          bcount),[4 1 2 3]);
   vprederr=permute(reshape(prederr,nlcount,attcount,...
                            respcount,bcount),[4 1 2 3]);
   vpredmse=permute(reshape(predmse,nlcount,attcount,...
                            respcount,bcount),[4 1 2 3]);
   vpredfix=permute(reshape(predfix,nlcount,attcount,...
                            respcount,bcount),[4 1 2 3]);
   vpredfixerr=permute(reshape(predfixerr,nlcount,attcount,...
                               respcount,bcount),[4 1 2 3]);
   
   predxc=mean(vpredxc,3);
   predinf=mean(vpredinf,3);
   predone=mean(vpredone,3);
   predp=mean(vpredp,3);
   prederr=mean(vprederr,3);
   predmse=mean(vpredmse,3);
   predfix=mean(vpredfix,3);
   predfixerr=mean(vpredfixerr,3);
   
   expxc=mean(vexpxc,4);
   
elseif params.fitboot,
   
   % predres ( batch x 1) . predxc ( nlcount*bootcount*respcount x 1)
   
   % vpredxc should be batch X nl X boot X resp
   % predxc should be batch X nl X 1 X resp   (averaged over boots!)
   
   vpredxc=permute(reshape(predxc,nlcount,attcount,respcount,...
                           bcount),[4 1 2 3]);
   vpredinf=permute(reshape(predinf,nlcount,attcount,respcount,...
                           bcount),[4 1 2 3]);
   vpredone=permute(reshape(predone,nlcount,attcount,respcount,...
                           bcount),[4 1 2 3]);
   vpredp=permute(reshape(predp,nlcount,attcount,respcount,...
                          bcount),[4 1 2 3]);
   vprederr=permute(reshape(prederr,nlcount,attcount,...
                            respcount,bcount),[4 1 2 3]);
   vpredmse=permute(reshape(predmse,nlcount,attcount,...
                            respcount,bcount),[4 1 2 3]);
   vpredfix=permute(reshape(predfix,nlcount,attcount,...
                            respcount,bcount),[4 1 2 3]);
   vpredfixerr=permute(reshape(predfixerr,nlcount,attcount,...
                               respcount,bcount),[4 1 2 3]);
   
   if batchcount==0,
      batchcount=1;
      rr=rrfull;
      pp=ppfull;
   else
      rr=rrfull{1};
      pp=ppfull{1};
   end
   predxc=permute(mean(reshape(vpredxc,batchcount,params.fitboot,nlcount),2),[1 3 2]);
   predinf=permute(mean(reshape(vpredinf,batchcount,params.fitboot,nlcount),2),[1 3 2]);
   predone=permute(mean(reshape(vpredone,batchcount,params.fitboot,nlcount),2),[1 3 2]);
   predp=permute(mean(reshape(vpredp,batchcount,params.fitboot,nlcount),2),[1 3 2]);
   prederr=permute(mean(reshape(vprederr,batchcount,params.fitboot,nlcount),2),[1 3 2]);
   predmse=permute(mean(reshape(vpredmse,batchcount,params.fitboot,nlcount),2),[1 3 2]);
   predfix=permute(mean(reshape(vpredfix,batchcount,params.fitboot,nlcount),2),[1 3 2]);
   predfixerr=permute(mean(reshape(vpredfixerr,batchcount,params.fitboot,nlcount),2),[1 3 2]);
   
   expxc=mean(vexpxc,4);

   %disp('untangling prediction data');
   %for predidx=1:length(rrfull),
   %   rr=rrfull{predidx};
   %   pp=ppfull{predidx};
   %   
   %   for ii=1:size(pp,2),
   %      kk=find(~isnan(rr)&~isnan(pp(:,ii)));
   %      [predxc(predidx,ii),prederr(predidx,ii),tt,predp(predidx,ii)]=randxcov(rr(kk),pp(kk,ii),0,100);
   %      vrr=var(rr(kk));
   %      if vrr>0,
   %         predmse(predidx,ii)=var(pp(kk,ii)-rr(kk))./vrr;
   %      end
   %   end
   %end
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

clear tstim tr th tH seplinpred rgoodidx respsave linpred tstrf
clear tmod_psth tgoodidx rpred rbadidx psthval lengthVec 
clear fstim fstim_spike fitidx CSR_ns CS_ns cdata tresp valpred

global SUBXC
if ~isempty(SUBXC),
   subxc=SUBXC;
   subxc(end,:,:)=squeeze(predxc);
end

% save results, zipped if necessary
daterun=date;
if params.zipoutfile,
   skernfile=basename(params.outfile);
   fprintf('cellxcnodb.m: SAVING to %s\n',params.outfile);
   
   % then save everything (inclusive to avoid future screw-ups
   % involving accidentally not saving new important stuff)
   save([tempdir,skernfile]);
   
   if params.showres,
      xcresult([tempdir,skernfile],[],1);
   end
   
   % compress and copy over network if necessary
   unix(['gzip -cf ',tempdir,skernfile,' > ',...
         params.outfile,'.gz']);
   delete([tempdir,skernfile]);
else
   save(params.outfile);
   if ~params.showres
      disp('skipping figure output');
   elseif params.batch>0,
      xcresult(params.outfile,[],1);
   else
      xcresult(params.outfile); % no batch--don't print
   end
end

disp('cellxcnodb.m: DONE');



