


cellid='modelhis';
batchid=24;


dbopen;
sql=['SELECT * from sRunData WHERE cellid="',cellid,'"',...
     ' AND batch=',num2str(batchid)];
rundata=mysql(sql);


[cellfiledata,times,batchdata]=cellfiletimes(cellid,batchid);
% fill up params structure for passing to cellxcnodb
params=dbget('sBatch','',batchid);
params.times=times;
params.cellid=cellid;

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
params.batch=batchid;

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
params.repexclude=0;
params.meansub=1;

params.zipoutfile=1;
params.outfile=[rundata.respath,rundata.resfile];

if length(batchdata.parmstring)>0,
   eval(char(batchdata.parmstring));
end

if isfield(params,'maxlag') & length(params.maxlag)>=2,
   % do nothing, this is a good format for running cellxc
else
   params.maxlag=[getparm(params,'minlag',-6) getparm(params,'maxlag',13)];
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

xcloadfiles;

maxlag=params.maxlag;

ract=load(params.respfiles{1});
htrue=reshape(ract.STRF.Phase0(:,:,1:(maxlag(2)+1)),256,maxlag(2)+1);

%
% BELOW THIS IS RANDOM CRAP THAT CAN BE MESSED WITH
%

%keyboard

hspacetrue=htrue(:,4);
htimetrue=htrue'*hspacetrue./(hspacetrue'*hspacetrue);
%stim=stim*hspacetrue;

resp=kernpredict(htrue,stim',1,0,1);


if 0,
   % only want natrev component for the time being
   resp=resp(1:7228);
   stim=stim(1:7228,:);
   
   
   htrue(:,6:end)=0;
   
   % fix stim to be gaussian white noise
   stim=randn(size(stim));
   resp=kernpredict(htrue,stim',1,0,1);
end

rgoodidx=find(~isnan(resp));

% number of valid time bins
n=length(rgoodidx);

% calculate mean response and stimulus
mR=sum(resp(rgoodidx))' ./ n;
mS=sum(stim(rgoodidx,:))' ./ n;

% calculate stimulus-response cross correlation
fprintf('Computing stimulus-response cross correlation');
SR=zeros(spacecount,diff(maxlag)+1);
for tt=maxlag(1):maxlag(2),
   fprintf('.');
   trg=rgoodidx;
   trg=trg(find(trg-tt>0 & trg-tt<size(stim,1)));
   SR(:,tt-maxlag(1)+1)=(resp(trg)'*stim(trg-tt,:))' ./ n - mS*mR;
end
fprintf('\n');

% stimulus 2nd order spatial autocorrelation
disp('Computing spatial autocorrelation...');
sSA2=stim(rgoodidx,:)'*stim(rgoodidx,:) - mS*mS';

% first order temporal autocorrelation
disp('Computing temporal autocorrelation...');
tSA=zeros(diff(maxlag)*2+1,1);
for xx=1:spacecount,
   tstim=stim(rgoodidx,xx);
   tSA=tSA + xcorr(tstim,diff(maxlag),'biased') ./ spacecount;
end
tSA=tSA - mS'*mS./spacecount;

tolval=[5e-6 1e-6 5e-7 1e-7 1e-8 1e-9 1e-10];
H=normalize(SR,sSA2,tSA,tolval);
sI=sSA2^-1;
tc=zeros(size(SR,1),size(SR,2));
for hh=1:size(SR,2),
   tc(:,hh)=sI*SR(:,hh);
end

H=cat(3,SR,tc,H,cat(2,zeros(size(H,1),-maxlag(1)),htrue(:,1:(maxlag(2)+1))));
Hcc=zeros(size(H,3),1);
tstim=stim'-repmat(mS,[1 size(stim,1)]);

for pc=1:length(Hcc),
   
   tH=H(:,-maxlag(1)+1:end,pc);
   
   predresp=kernpredict(tH,tstim,1,0);  % 1 phase, no rect
   vidx=find(~isnan(resp) & ~isnan(predresp));
   cc=corrcoef(predresp(vidx),resp(vidx));
   Hcc(pc)=cc(1,2);
   
end

showkern(H(:,-maxlag(1)+1:end,:),'space')
Hcc

