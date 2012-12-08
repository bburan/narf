% challdemo.m
%
% test prediction challenge data with cellxcnodb.m
%
% path requirements:
%   ... /cellxc/
%   ... /mutils/
%   ... /gen/
%   ... /db/
%
% created SVD 6/21/05 -- ripped off of cellxcdemo
%
function challdemo(cellid,kernfmt,dataset)

% make sure parameters are being fully reset
clear params

% no boosting
BOOSTING=1;

if ~exist('cellid','var'),
   cellid='r0221a';
end
if ~exist('dataset','var'),
   dataset='review';
end

params.cellid=cellid;
if strcmpi(dataset,'review'),
   datapath='/auto/b0/data/npc/data/01_v1_nvm/';
   params.stimfiles={sprintf('%s%s_data.mat',datapath,cellid),...
                     sprintf('%sval_data/%s_val.mat',datapath,cellid)};
   params.respfiles={sprintf('%s%s_data.mat',datapath,cellid),...
                     sprintf('%sval_data/%s_val.mat',datapath,cellid)};
else
   datapath='/auto/data/nsl/users/svd/data/v1_nis_data/';
   params.stimfiles={sprintf('%s%s_data.mat',datapath,cellid)};
   params.respfiles={sprintf('%s%s_data.mat',datapath,cellid)};
end

% estimate the model kernel
params.resploadcmd='respload';
params.resploadparms={'resp',0,1,1};
params.stimloadcmd='loadchallstim';
params.stimloadparms={};

if ~exist('kernfmt','var'),
   params.kernfmt='pfft';
else
   params.kernfmt=kernfmt;
end

if strcmp(params.kernfmt,'space'),
   params.stimfiltercmd='';
   params.stimfilterparms={};
elseif strcmp(params.kernfmt,'pfft'),
   params.stimfiltercmd='movpower';
   params.stimfilterparms={0,0,1,1,50};
else
   params.stimfiltercmd='movphasesep';
   params.stimfilterparms={0,0,1,1,0};
end

% this runs slower but produces less biased strf estimates
params.altfit='xcfit2';

if ~BOOSTING,
   % no boosting
   params.outfile=['/auto/b0/svd/npc/chall-',cellid,'-',...
                   params.kernfmt,'-',dataset,'.mat'];
   params.sfscount=10;
else
   % boosting
   params.minlag=0;
   params.maxlag=10;
   params.altcore='cdcore';
   params.outfile=['/auto/b0/svd/npc/chall-',cellid,'-bs-',params.kernfmt,'-',dataset,'.mat'];
   params.sfscount=15;
end

params.showres=1;
params.sfsstep=7;
params.sffiltsigma=5;
params.smoothtime=0;
params.fitfrac=0;
params.predfrac=0.075;

cellxcnodb(params);

disp('generating npc-compatible prediction');
res=load(params.outfile);

data=[];
load(res.params.stimfiles{1},'vstim');
data.stim=vstim;
data.resp=zeros(size(data.stim,3),1);

if ~isempty(res.params.stimfiltercmd),
   data.stim=feval(res.params.stimfiltercmd,data.stim,...
                   res.params.stimfilterparms{:})';
else
   sds=size(data.stim);
   data.stim=reshape(data.stim,sds(1)*sds(2),sds(3))';
end

predres=xcval(res.strf(2),res.params,data);

pp=predres.mod_psth{1}(:,1,1);
cheatpp=predres.act_resp{1}(:,1);

pp(find(isnan(pp)))=0;
cheatpp(find(isnan(cheatpp)))=0;

if BOOSTING,
   npcfile=['/auto/b0/svd/npc/realpreds-bs-',params.kernfmt,...
            '-',dataset,'.mat'];
else
   npcfile=['/auto/b0/svd/npc/realpreds-',params.kernfmt,...
            '-',dataset,'.mat'];
end
if exist(npcfile,'file'),
   load(npcfile);
else
   prediction=[];
end

matched=0;
for nn=1:length(prediction),
   if strcmp(prediction(nn).cellid,cellid),
      matched=1;
      prediction(nn).response=pp;
   end
end
if ~matched,
   nn=length(prediction)+1;
   prediction(nn).cellid=cellid;
   prediction(nn).response=pp;
end
save(npcfile,'prediction');

if 0
   % only need to have run this once
   load /auto/b0/svd/npc/cheatpreds.mat
   matched=0;
   for nn=1:length(prediction),
      if strcmp(prediction(nn).cellid,cellid),
         matched=1;
         prediction(nn).response=pp;
      end
   end
   if ~matched,
      nn=length(prediction)+1;
      prediction(nn).cellid=cellid;
      prediction(nn).response=cheatpp;
   end
   save /auto/b0/svd/npc/cheatpreds.mat prediction
end

return

cd /auto/b0/data/npc/data/01_v1_nvm
cd /auto/data/nsl/users/svd/data/v1_nis_data/
cellinfo
for ii=1:length(celldata),
   dbaddqueuemaster(['challdemo(''',celldata(ii).cellid,''',''pfft'',''natrev'');'],...
                    [celldata(ii).cellid,'/pfft']);
end


cd /auto/b0/data/npc/data/01_v1_nvm
cellinfo
for ii=1:length(celldata),
   cellid=celldata(ii).cellid;
   outfile=['/auto/b0/svd/npc/chall-',cellid,'-bs-pfft.mat'];
   res=load(params.outfile);
   
   data=load(res.params.stimfiles{2});
   if ~isempty(res.params.stimfiltercmd),
      data.stim=feval(res.params.stimfiltercmd,data.stim,...
                      res.params.stimfilterparms{:})';
   else
      sds=size(data.stim);
      data.stim=reshape(data.stim,sds(1)*sds(2),sds(3))';
   end
   
   predres=xcval(res.strf(2),res.params,data);
   
   pp=predres.mod_psth{1}(:,1,1);
   cheatpp=predres.act_resp{1}(:,1);
   
   pp(find(isnan(pp)))=0;
   cheatpp(find(isnan(cheatpp)))=0;
   
   npcfile=['/auto/b0/svd/npc/realpreds-bs-',params.kernfmt,'.mat'];
   load(npcfile);
   matched=0;
   for nn=1:length(prediction),
      if strcmp(prediction(nn).cellid,cellid),
         matched=1;
         prediction(nn).response=pp;
      end
   end
   if ~matched,
      nn=length(prediction)+1;
      prediction(nn).cellid=cellid;
      prediction(nn).response=pp;
   end
   save(npcfile,'prediction');
   
end







%
%  video demo
%

% estimate the model kernel
params.stimfiles={'/auto/k5/hayden/van_data/2003-11-23/v0098.video4.007.imsm'};
params.respfiles={'/auto/k5/hayden/van_data/2003-11-23/v0098.video4.007.resp.mat'};
params.resploadcmd='respload';
params.resploadparms={'',0,1,1};
params.stimloadcmd='loadimfile';
params.stimloadparms={32,0,16};

params.stimfiltercmd='movpower';
params.stimfilterparms={0,0,1,1,0};
params.kernfmt='pfft';

% this may help but runs slower
params.altfit='xcfit2';

params.outfile=['/tmp/chall.',cellid,'.mat'];
params.showres=1;
params.sfsstep=8;
params.sfscount=20;
params.sffiltsigma=5;
params.smoothtime=0;
params.fitfrac=0;
params.predfrac=0.1;

cellxcnodb(params);

