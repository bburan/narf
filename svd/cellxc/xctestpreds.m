% function res=xctestpreds(runidx,batch,nlidxfortest,batchfortest);
%
% figure out asymptotic prediction scores based on fits to single/finite
% trial predictions.
%
% the idea is that noise goes as 1/N * sigma(neuron)
% that is, it shrinks linearly as one over the number of
% trials. when 1/N goes to 0 (ie, infinite trials), the noise
% should be zero. any residual error is the quality of the model
% doing the prediction
%
% returns:
% res.asympxc=asympxc;
% res.onetrialxc=onetrialxc;
% res.xc=xc;
% res.asympmse=asympmse;
% res.onetrialmse=onetrialmse;
% res.mse=mse;
% res.trialcount=n(end)
% res.batch=batch;
%
function res=xctestpreds(runidx,batch,nlidxfortest,batchfortest);

disp('xctestpreds.m:');

global BATQUEUEID

persistent pdata

dbopen;
if isnumeric(runidx),
   sql=['SELECT * from sRunData WHERE id=',num2str(runidx)];
   rundata=mysql(sql);
   cellid=rundata.cellid;
else
   cellid=runidx;
   goodbatch=zeros(1,length(batch));
   sql=['SELECT * from sRunData WHERE cellid="',cellid,'"',...
        ' AND batch=',num2str(batch)];
   rundata=mysql(sql);
end

if length(rundata)==0,
   disp('no entries found in db!');
   if nargout>0,
      r=0;
   end
   return
end

outfile=[rundata(1).respath,rundata(1).resfile,'.gz'];

if strcmp(outfile(end-2:end),'.gz'),
   t = tempname;
   [s, w] = unix(sprintf('gunzip <%s >%s', outfile, t));
   if s,
      fprintf('can''t find %s\n', fname);
      z=[];
      return
   end
   
   load(t, '-mat');
   delete(t);
else
   load(outfile);
end

BATQUEUEID=[];
if exist('nlidxfortest','var'),
   nlidx=nlidxfortest;
else
   nlidx=params.nlidxsave;
end
if nlidx>size(strf,1),
   nlidx=1;
end
fprintf('chose to show nlidx=%d\n',nlidx);
strfcount=size(strf,3);
fitcount=params.fitboot;

predcount=strfcount*fitcount;

batch=cat(1,params.predbatch{:});
asympxc=zeros(fitcount,strfcount,length(batch));
onetrialxc=zeros(fitcount,strfcount,length(batch));
xc=zeros(fitcount,strfcount,length(batch));
asympmse=zeros(fitcount,strfcount,length(batch));
onetrialmse=zeros(fitcount,strfcount,length(batch));
mse=zeros(fitcount,strfcount,length(batch));
asymprr=zeros(fitcount,strfcount,length(batch));
onetrialrr=zeros(fitcount,strfcount,length(batch));
ceilxc=zeros(fitcount,strfcount,length(batch));
oneceilxc=zeros(fitcount,strfcount,length(batch));

clear predres

predidx=batchfortest;
pdata.alive=1;

fprintf('Predicting responses to batch %d.\n',params.predbatch{predidx});

% figure out pred files for the current batch
[pcellfiledata,ptimes,pbatchdata]=...
    cellfiletimes(params.cellid,params.predbatch{predidx});

% does this cell have data for batchid=predidx?
if length(pcellfiledata)==0,
   disp('no pred data!');
   res=[];
   return
end

% load response raster and check to see if it's big enough
ro=respload(predparams.respfiles{tpredfile});
ro=ro(tpredstartframe:tpredstopframe,:);
r=compact_raster_matrix3(ro(:,2:end));      

if size(r,1)>2 & size(r,2)>2,
   fprintf('only one trial in cnf data!\n');
   res=[];
   return
end

% stim file and filtering info
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

if isfield(pdata,'respfile') & ...
      strcmp(pdata.respfile,predparams.respfiles{tpredfile}),
   
   % keep current pdata
   disp('keeping previous pdata');
   
else
   
   % load validation stim and resp
   [pdata.stim,pdata.resp]=...
       xcloadstimresp(tpredfile,tpredstartframe,...
                      tpredstopframe,predparams);
   pdata.respfile=predparams.respfiles{tpredfile};
end

% do standard predictions to get mod_psth
tstrf=reshape(strf(nlidx,:,1),fitcount,1);
predres=xcval(tstrf,predparams,pdata);

ro=respload(predparams.respfiles{tpredfile});
ro=ro(tpredstartframe:tpredstopframe,:);
r=compact_raster_matrix3(ro(:,2:end));      

stimstart=max([tpredstartframe-params.maxlag(2) 1]);
lagstart=tpredstartframe-stimstart;
lagend=size(pdata.stim,1)-size(ro,1)-lagstart;
respcount=size(r,2);

r=cat(1,ones(lagstart,respcount)*nan,r,...
      ones(lagend,respcount)*nan);
fprintf('%d repeats of validation data\n',size(r,2));
r2=r;

%r2=ones(size(ro,1),size(r,2)).*nan;
%r2(find(ro(:,1)>-1),:)=r;

pp=predres.mod_psth{1}(:,:,nlidx);
rr=predres.act_resp{1};
okidx=find(~isnan(pp) & ~isnan(rr));

r2=r2(okidx,:);
r2=r2(:,find(sum(isnan(r2),1)==0));
pp=pp(okidx);
rr=rr(okidx);

if params.predsmoothsigma > 0 & params.respfmtcode==0,
   %fprintf('Smoothing actual response with optimal filter...\n');
   %pfilt=[1/9 2/9 1/3 2/9 1/9]';
   
   fprintf('(rr: smsig=%.1f)',params.predsmoothsigma);
   tt=(-10:10)';
   pfilt=exp(-tt.^2/(2*params.predsmoothsigma.^2))./...
         (sqrt(2*pi)*params.predsmoothsigma);
   pfilt=pfilt(find(pfilt > 0.1*max(pfilt)));
   fr=r2;
   fr(isnan(fr))=0;
   fr=conv2(fr,pfilt,'same');
   r2(~isnan(r2))=fr(~isnan(r2));
end
fprintf('\n');

[ccr,ccrM]=frednoise(r2);
M=size(r2,2);

ncount=round((respcount+1)/2);
n=zeros(ncount,fitcount);
s=zeros(ncount,fitcount);
m=zeros(ncount,fitcount);
cxy=zeros(ncount,respcount,fitcount);

fprintf('fitidx:');
strfidx=1;
for fitidx=1:fitcount,
   fprintf(' %d',fitidx);
   
   pp=predres.mod_psth{1}(okidx,:,fitidx);
   [m(:,fitidx),s(:,fitidx),n(:,fitidx),cxy(:,:,fitidx)]=predtrialcurve(pp,r2);
end
fprintf('\n');


% save important results
clear res
res.n=n;
res.m=m;
res.s=s;
res.cxy=cxy;
res.ccr=ccr;
res.ccrM=ccrM;

return

% display xc results
cla
mn=n(:,1);
mm=mean(m,2);
me=mean(s,2);
errorbar(mn,mm,me);

if 0
   useidx=floor(length(n)./2):length(n);
   
   if sum(n(useidx)==0)>0 | sum(m(useidx)==0)>0,
      disp('ZEEERO predicted response!');
      xcparms=[0 0];
      mfit=m;
   else
      
      % y axis is predxc
      xcparms=polyfit(1./n(useidx),1./m(useidx),1);
      mfit=1./(1./n .* xcparms(1) + xcparms(2));
      asympxc(fitidx,strfidx,predidx)=1./xcparms(2);
      xc(fitidx,strfidx,predidx)=m(end);
      onetrialxc(fitidx,strfidx,predidx)=m(1);
   end
end


