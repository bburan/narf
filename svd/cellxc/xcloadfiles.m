% xcloadfiles.m
%
% script to load data from multiple stimulus/response file sets and
% preprocess into standard rc format. ie, outputs are: 
%
%    stim X x T
%    resp T x M
%
% created SVD 1/23/03 - ripped off of cellxc.m
%

disp('xcloadfiles.m:');

if ~exist('attidx','var'),
   attidx=1;
end

% movstep: size of movie segments to send to movxc -- want this to
% be as large as possible without straining computer's memory
movstep=200000;

fidx=1;
firstseg=1;

% initialize data matrices
curlen=0;
stim=[];
resp=[];

% while there are still data files available and we want stim/resp
% data from the current file and we haven't gone over the maximum
% length, iterate through the available data

filecount=length(params.respfiles);
stimstart=zeros(filecount,1);
stimstop=zeros(filecount,1);

tresploadparms=params.resploadparms;
if strcmp(params.resploadcmd,'loadspikeraster') && ...
      isstruct(tresploadparms{1}) && ...
      tresploadparms{1}.psthonly==1,
   disp('saving raster in extras');
   tresploadparms{1}.psthonly=0;
   savingraster=1;
elseif strcmp(params.resploadcmd,'loadspikeraster') && ...
      length(tresploadparms)==5 && ...
      tresploadparms{5}==1,
   disp('saving raster in extras');
   tresploadparms{5}=0;
   savingraster=1;
else
   savingraster=0;
end

while fidx<=filecount & ...
            (stoptimes(fidx)-starttimes(fidx)>0 | isnan(stoptimes(fidx))) & ...
            (curlen==0 | stoptimes(fidx)-starttimes(fidx)+1+curlen<movstep),
   
   clear tstim
   % load next stim file
   tstimloadparms=params.stimloadparms;
   if ((strcmp(params.stimloadcmd,'loadimfile') |...
        strcmp(params.stimloadcmd,'loadimscaled')) & ...
       length(params.stimloadparms)>0 & params.stimloadparms{1}==0) ...
      tstimloadparms{1}=round(params.stimloadparms{3} .* ...
                              params.stimcrfs(fidx)./params.stimwindowcrf);
   end
   if (strcmp(params.stimloadcmd,'loadimfix') & ...
       length(params.stimloadparms)>0 & params.stimloadparms{2}==0),
      tstimloadparms{2}=round(params.stimloadparms{4} .* ...
                              params.stimcrfs(fidx)./params.stimwindowcrf);
   end
   
   if ~isnan(stoptimes(fidx,attidx)) & ...
         (strcmp(params.stimloadcmd,'loadimfile') |...
          strcmp(params.stimloadcmd,'loadimscaled')),
      stimstart(fidx)=max([starttimes(fidx,attidx)-params.maxlag(2) 1]);
      
      framecount=imfileinfo(params.stimfiles{fidx},1);
      stimstop(fidx)=min([stoptimes(fidx,attidx)-params.maxlag(1) framecount]);
   elseif ~isnan(stoptimes(fidx,attidx)),
      stimstart(fidx)=max([starttimes(fidx,attidx)-params.maxlag(2) 1]);
      stimstop(fidx)=stoptimes(fidx,attidx)-params.maxlag(1);
   else
      stimstop(fidx)=0;
   end
   
   fprintf('%d: %d-->%d now %d-->%d\n',fidx,starttimes(fidx,attidx),...
           stoptimes(fidx,attidx),stimstart(fidx),stimstop(fidx));
   if strcmp(params.stimloadcmd,'loadstimfrombaphy'),
      [tstim,stimparam]=...
          feval(params.stimloadcmd,params.stimfiles{fidx},stimstart(fidx),...
                stimstop(fidx),tstimloadparms{:});
      params.stimparam=stimparam;
   else
      tstim=feval(params.stimloadcmd,params.stimfiles{fidx},...
                  stimstart(fidx),stimstop(fidx),tstimloadparms{:});
   end
   
   % apply stimulus filter
   if ~isempty(params.stimfiltercmd),
      tstim=feval(params.stimfiltercmd,tstim,params.stimfilterparms{:});
   end
   if issparse(tstim),
      tstim=full(tstim);
   end
   
   % reshape to space X time if necessary
   iconside=size(tstim);
   bslen=iconside(end);
   if bslen<stimstop(fidx)-stimstart(fidx)+1,
      stimstop(fidx)=bslen+stimstart(fidx)-1;
   end
   
   iconside=iconside(1:(end-1));
   if length(iconside)>1,
      % reshape all spatial dims into one
      tstim=reshape(tstim,prod(iconside),bslen);
   end
   
   % save for use in other stuff;
   params.iconside=iconside;
   
   % take transpose to put time in rows and append to stim
   stim=cat(1,stim,tstim');
   
   % load next respfile
   disp(sprintf('Loading response: %s...',params.respfiles{fidx}));
   
   
   % resp=resploadatt(respfile,respvarname);
   % resp is time X space(phase) X attcode
   tresp=feval(params.resploadcmd,params.respfiles{fidx},...
               tresploadparms{:});

   if size(tresp,1)==1,
      tresp=tresp';
   end
   % filter response (not used yet)
   if ~isempty(params.respfiltercmd),
      tresp=feval(params.respfiltercmd,tresp,params.respfilterparms{:});
   end
   rsize=size(tresp);
   
   % trim response to appropriate time range
   if ~isnan(stoptimes(fidx,attidx)),
      if rsize(1)<stimstop(fidx),
         tresp=cat(1,tresp(stimstart(fidx):end,:,:),...
                   ones([stimstop(fidx)-rsize(1) rsize(2:end)]));
      else
         tresp=tresp(stimstart(fidx):stimstop(fidx),:,:);
      end
   end
   
   tresp=tresp(:,:); % reshape to 2D matrix
   respcount=size(tresp,2);
   tresp(1:params.maxlag(2),:)=nan;
   tresp(length(tresp)+params.maxlag(1):end,:)=nan;
   
   if ((isfield(params,'smoothexp') & params.smoothexp)) & ...
         params.respfmtcode==0,
      
      % if it's single trial data, smooth with optimal linear filter
      runique=unique(tresp(~isnan(tresp(:,1)),1));
      if length(runique)<(size(tresp,1)/20),
         disp('single trial smoothing');
         smoothsig=1;
         tt=(-10:10)';
         pfilt=exp(-tt.^2/(2*smoothsig.^2))./(sqrt(2*pi)*smoothsig);
         pfilt=pfilt(find(pfilt > 0.1*max(pfilt)));
         pfilt=[0.25; .5; 0.25];
         fr=tresp;
         fr(isnan(fr))=0;
         fr=conv2(fr,pfilt,'same');
         tresp(~isnan(tresp))=fr(~isnan(tresp));
      end
   end
   
   % match number of reps/att conditions by padding with nans in
   % data with fewer reps.
   if size(resp,2)>size(tresp,2),
      tresp(:,size(tresp,2)+1:size(resp,2))=nan;
   elseif size(resp,2)<size(tresp,2) & size(resp,1)>0,
      resp(:,size(resp,2)+1:size(tresp,2))=nan;
   end
   
   resp=cat(1,resp,tresp);
   
   % done, keep track
   curlen=curlen+stimstop(fidx)-stimstart(fidx)+1;
   
   % update queue if active
   if exist('BATQUEUEID','var') & BATQUEUEID>0,
      dbsetqueue(BATQUEUEID,1);
   end
   
   fidx=fidx+1;
end

fprintf('stim is %d x %d. resp is %d x %d\n',size(stim,1),size(stim,2),...
        size(resp,1),size(resp,2));

spacecount=size(stim,2);
if savingraster,
   raster=resp;
   resp=nanmean(resp,2);
   respcount=size(resp,2);
else
   raster=[];
end
rsize=size(resp);

if ~params.meansub,
   disp('adding constant channel to stim');
   stim(:,size(stim,2))=1;
   
   %disp('subtracting mean resp!!!');
   %for ii=1:size(resp,2),
   %   resp(:,ii)=resp(:,ii)-nanmean(resp(:,ii));
   %end
end

clear tstim



