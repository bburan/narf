% function [stim,resp,extras]=xcloadstimresp(fileidx,starttimes,stoptimes,
%                                     params,attidx);
%
% load stimulus/response data sets
%
% fileidx= vector of indiced into respfile/stimfile list to load
% starttimes,stoptimes - for each file in fileidx, the first and
%                        last time bins to load
% params - structure with a bunch of useful parameters:
%    these fields are required (no defaults provided):
%          .stimfiles={'/path/to/stim1.imsm','/path/to/stim2.imsm',...}
%          .respfiles={'/path/to/resp1.mat','/path/to/resp2.mat',...}
%          .resploadcmd    [='respload']
%          .resploadparms  [={'',1,1,1}]
%          .respfiltercmd  [='']
%          .respfilterparms [={}];
%          .stimloadcmd    [='loadimfile']
%          .stimloadparms  [={'',1,1,1}]
%          .stimfiltercmd  [='']
%          .stimfilterparms [={}];
%    optional:
%          .maxlag [= [-6 13] ] -- puts a buffer of nan equal to
%                                  13+6 between each file
%
% returns:
%   stim T x N matrix
%   resp T x 1 vector
%   extras.origstimsize = [X Y ...] vector of original stim
%                         dimensions in each frame (where XxY = N)
%
function [stim,resp,extras]=xcloadstimresp(fileidx,starttimes,stoptimes,...
                                           params,attidx);

global BATQUEUEID

if ~exist('attidx','var'),
   attidx=1;
end

stim=[];
resp=[];
curlen=0;
uselens=stoptimes-starttimes;
fileidx=fileidx(uselens>1);

% maxlag is used for padding with nans between stim sets.
if ~isfield(params,'maxlag'),
   params.maxlag=[0 0];
end
if params.maxlag(1)>0,
   params.maxlag(1)=0;
end
if params.maxlag(2)<0,
   params.maxlag(2)=0;
end

for fidx=1:length(fileidx(:)),
   
   % load next stim file
   tstimloadparms=params.stimloadparms;
   if (strcmp(params.stimloadcmd,'loadimfile') |...
       strcmp(params.stimloadcmd,'loadimscaled')) & ...
         length(params.stimloadparms)>0 & params.stimloadparms{1}==0,
      tstimloadparms{1}=round(params.stimloadparms{3} .* ...
                params.stimcrfs(fileidx(fidx))./params.stimwindowcrf);
   end
   if (strcmp(params.stimloadcmd,'loadimfix') & ...
       length(params.stimloadparms)>0 & params.stimloadparms{2}==0),
      tstimloadparms{2}=round(params.stimloadparms{4} .* ...
                              params.stimcrfs(fidx)./params.stimwindowcrf);
   end
   stimstart=max([starttimes(fidx)-params.maxlag(2) 1]);
   % display start and stop frames
   %[stimstart stoptimes(fidx)-params.maxlag(1)]
   tstim=feval(params.stimloadcmd,params.stimfiles{fileidx(fidx)},stimstart,...
               stoptimes(fidx)-params.maxlag(1),tstimloadparms{:});
   if ~isempty(params.stimfiltercmd),
      tstim=feval(params.stimfiltercmd,tstim,params.stimfilterparms{:});
   end
   if isfield(params,'stimdumpfile'),
      [framecount1,iconside1,stimfile,path,filetype1]=imfileinfo(params.stimfiles{fileidx(fidx)});
      framestep=1000;
      n=stimstart-1;
      frames1=min([framecount1 stoptimes(fidx)-params.maxlag(1)-stimstart+1]);
      for clumpidx=1:ceil(frames1/framestep),
         nstart=n+1;
         nstop=min([n+framestep stoptimes(fidx)-params.maxlag(1)]);
         mov=loadimfile(params.stimfiles{fileidx(fidx)},nstart,nstop);
         n=n+size(mov,3);
         appendimfile(mov,params.stimdumpfile);
      end
      clear mov
   end
   
   if issparse(tstim),
      tstim=full(tstim);
   end
   
   % reshape to space X time if necessary
   iconside=size(tstim);
   bslen=iconside(end);
   iconside=iconside(1:(end-1));
   if length(iconside)>1,
      % reshape all spatial dims into one
      tstim=reshape(tstim,prod(iconside),bslen);
   end
   
   clear fun
   
   % take transpose to put time in rows and append to stim
   stim=cat(1,stim,tstim');
   
   % load next respfile
   disp(sprintf('Loading response: %s...',params.respfiles{fileidx(fidx)}));
   
   % resp=resploadatt(respfile,respvarname);
   % resp is time X space(phase) X attcode
   if strcmp(params.resploadcmd,'loadspikeraster') && ...
         isstruct(params.resploadparms{1}) && ...
         params.resploadparms{1}.psthonly==1,
      disp('saving raster in extras');
      params.resploadparms{1}.psthonly=0;
      savingraster=1;
   elseif strcmp(params.resploadcmd,'loadspikeraster') && ...
         length(params.resploadparms)==5 && ...
         params.resploadparms{5}==1,
      disp('saving raster in extras');
      params.resploadparms{5}=0;
      savingraster=1;
   else
      savingraster=0;
   end
   
   tresp=feval(params.resploadcmd,params.respfiles{fileidx(fidx)},...
               params.resploadparms{:});
   
   % filter response (not used yet)
   if ~isempty(params.respfiltercmd),
      tresp=feval(params.respfiltercmd,tresp,params.respfilterparms{:});
   end
   
   % isolate just a single attentional set per run... reduces load
   tresp=tresp(:,:,attidx);
   if size(tresp,1)==1,
      tresp=tresp';
   end
   
   % trim response to appropriate time range
   tresp=tresp(starttimes(fidx):stoptimes(fidx),:);
   
   if ((isfield(params,'smoothexp') & params.smoothexp)) & ...
         params.respfmtcode==0,
      
      % if it's single trial data, smooth with optimal linear filter
      runique=unique(tresp(~isnan(tresp(:,1)),1));
      if length(runique)<(size(tresp,1)/20),
         smoothsig=params.smoothexp;
         fprintf('single trial smoothing. sigma=%.1f\n',smoothsig);
         tt=(-10:10)';
         pfilt=exp(-tt.^2/(2*smoothsig.^2))./(sqrt(2*pi)*smoothsig);
         pfilt=pfilt(find(pfilt > 0.1*max(pfilt)));
         %pfilt=[0.25; .5; 0.25]';
         fr=tresp;
         fr(isnan(fr))=0;
         fr=conv2(fr,pfilt,'same');
         tresp(~isnan(tresp))=fr(~isnan(tresp));
      end
   end
   
   rsize=size(tresp);
   tresp=tresp(:,:); % reshape to 2D matrix
   tresp(1:params.maxlag(2),:)=nan;
   lagstart=starttimes(fidx)-stimstart;
   lagend=size(tstim,2)-size(tresp,1)-lagstart;
   
   
   if fidx>1 & size(resp,2)>size(tresp,2),
      tresp(:,size(tresp,2)+1:size(resp,2))=nan;
   elseif fidx>1,
      resp(:,size(resp,2)+1:size(tresp,2))=nan;
   end
   respcount=size(tresp,2);
   resp=cat(1,resp,ones(lagstart,respcount)*nan,tresp,...
            ones(lagend,respcount)*nan);
   
   % done, keep track
   curlen=curlen+stoptimes(fidx)-starttimes(fidx)+1;
   
   % update queue if active
   if (exist('dbsetqueue')==2),
       dbsetqueue;
   end  
end

if savingraster,
   extras.raster=resp;
   resp=nanmean(resp,2);
end

extras.origstimsize=iconside;


