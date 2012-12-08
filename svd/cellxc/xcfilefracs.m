% function [times,params]=xcfilefracs(params);
%
% you need to provide some params (* required, [default]):
% params.stimfiles={'/path/to/stim1.imsm','/path/to/stim2.imsm',...} *
%       .respfiles={'/path/to/resp1.mat','/path/to/resp2.mat',...} *
%       .resploadcmd    [='respload']
%       .resploadparms  [={'',1,1,1}]
%       .respfiltercmd  [='']
%       .respfilterparms [={}];
%       .predfrac [=0]  % frac of data for validation
%       .fitfrac [=0]   % fraction for intermediate fit
%       .expfrac [=1-predfrac-fitfrac]   % fraction for main fit
%
% returns times() structure
%   times(1) : describes data for main fit
%   times(2) : describes data for supplementary (meta) fit
%   times(3) : describes data for validataion
%   times( ).fileidx : indexes entries in params.respfiles/.stimfiles
%           .start   : time bin to start in each file
%           .stop    : time bin to stop
%
% created SVD 3/03 - ripped off cellfiletimes
% modified SVD 9/04 - removed multiple att cond functionality
%                   - fixed val data selection to choose stimulus
%                     that was repeated many times
%                   - merge repeated files into a single "file"
% 
function [times,params]=xcfilefracs(params);

MAXPREDLEN=70000; % plenty data for fitting regularization/threshold parms?

%dbopen;

% check params and set default values where applicable
respfiles=getparm(params,'respfiles','ERROR');
stimfiles=getparm(params,'stimfiles','ERROR');
expfrac=getparm(params,'expfrac',0);
fitfrac=getparm(params,'fitfrac',0);
predfrac=getparm(params,'predfrac',0.1);
decorrspace=getparm(params,'decorrspace',2);
resploadcmd=getparm(params,'resploadcmd','respload');
resploadparms=getparm(params,'resploadparms',{'',1,1,1});
respfiltercmd=getparm(params,'respfiltercmd','');
respfilterparms=getparm(params,'respfilterparms',[]);
resampcount=getparm(params,'resampcount',20);

if 0,
    dbopen;
   sql=['SELECT concat(stimpath,stimfile) as stim,',...
        ' concat(path,respfile) as resp, stimfilecrf',...
        ' FROM sCellFile where cellid="e0107" and runclassid=2'];
   filedata=mysql(sql);
   filecount=length(filedata);
   for fidx=1:filecount,
      params.respfiles{fidx}=filedata(fidx).resp;
      params.stimfiles{fidx}=filedata(fidx).stim;
      params.stimcrfs(fidx)=filedata(fidx).stimfilecrf;
   end
   params.resploadcmd='respload';
   params.resploadparms={'',1,1,1};
   params.fitfrac=0.1;
   params.predfrac=0.1;
   params.respfiltercmd='';
   params.respfilterparms={};
   params.resampcount=20;
end

if length(respfiles)==0 | length(stimfiles)==0,
   times=[];
   return
end

filecount=length(respfiles);
keepidx=ones(filecount,1);
for fidx=filecount:-1:2,
   for f2idx=1:fidx-1,
      if strcmp(stimfiles{fidx},stimfiles{f2idx}) & ...
            keepidx(fidx) & keepidx(f2idx),
         respfiles{f2idx}=[respfiles{f2idx},'+',respfiles{fidx}];
         respfiles{fidx}='';
         keepidx(fidx)=0;
         if isfield(params,'repcount'),
            params.repcount(f2idx)=...
                params.repcount(fidx)+params.repcount(f2idx);
         end
      end
   end
end

keepidx=find(keepidx);
respfiles={respfiles{keepidx}};
stimfiles={stimfiles{keepidx}};
if isfield(params,'revstimresp') && params.revstimresp,
   params.respfiles=stimfiles;
   params.stimfiles=respfiles;
else
   params.respfiles=respfiles;
   params.stimfiles=stimfiles;
end
if isfield(params,'resplen'),
   params.resplen=params.resplen(keepidx);
end
if isfield(params,'repcount'),
   params.repcount=params.repcount(keepidx);
end
if isfield(params,'stimcrfs'),
   params.stimcrfs=params.stimcrfs(keepidx);
end
filecount=length(params.respfiles);

% figure out separate start and stop times for each attentional state!
resplens=zeros(filecount,1);
attcount=1;
rvalid={};

for fidx=1:filecount,
   resp=feval(resploadcmd,params.respfiles{fidx},resploadparms{:});
   % filter response (eg resample, pick attentional state, etc)
   if ~isempty(respfiltercmd),
      resp=feval(respfiltercmd,resp,respfilterparms{:});
   end
   
   rsize=size(resp);
   if rsize(1)==1,
      resp=resp(:);
      rsize=size(resp);
   end
   if rsize(2)>1,
      resp=nanmean(resp')';
   end
   
   rvalid{fidx}=find(~isnan(resp(:,1)));
   resplens(fidx)=length(rvalid{fidx});
end

starttimes=ones(filecount,1);
stoptimes=zeros(filecount,1);
cumresplens=cumsum(resplens);
totlen=sum(resplens);

if predfrac>0,
   predstartframe=round(totlen*(1-predfrac));
   % min pred is 200 frames
   if totlen*predfrac<200 & totlen>400,
      predstartframe=totlen-199;
   end
   
   predfile=sum(cumresplens<predstartframe,1)+1;
   % make sure that last file is getting used in pred file
   if predfile<filecount,
      predfile=filecount;
      predstartframe=sum(resplens(1:(filecount-1)))+1;
   end
else
   predstartframe=totlen;
   predstopframe=totlen;  
   predfile=filecount;
end

predstartframe=predstartframe - sum(resplens(1:(predfile-1)));
predstopframe=min([predstartframe+MAXPREDLEN resplens(predfile)]);

% start fit region offset from pred region
tpredstart=predstartframe+sum(resplens(1:predfile-1));
fitstartframe=round(tpredstart-totlen*fitfrac);
fitfile=sum(cumresplens<fitstartframe,1)+1;
if fitfile>filecount,
   fitfile=filecount;
end
if fitfile<filecount-1,
   fitfile=filecount-1;
   fitstartframe=cumresplens(filecount-2)+1;
elseif cumresplens(fitfile)-fitstartframe < floor(totlen*fitfrac/2) & ...
      fitfile<filecount,
   fitfile=fitfile+1;
   fitstartframe=cumresplens(fitfile-1)+1;
end

if expfrac<=0,
   texpfrac=1-fitfrac-(predstopframe-predstartframe)./totlen;
   expstopframe=round(totlen*texpfrac)-1;
else
   expstopframe=round(totlen*expfrac)-1;
end
expstopfile=sum(cumresplens<expstopframe,1)+1;
expstopframe=expstopframe - sum(resplens(1:(expstopfile-1)));
if expstopframe<resampcount & expstopfile>1,
   expstopfile=expstopfile(attidx)-1;
   expstopframe=resplens(expstopfile);
end
fitstartframe=fitstartframe - sum(resplens(1:(fitfile-1)));

for fidx=1:fitfile,
   if resplens(fidx)>0,
      starttimes(fidx)=rvalid{fidx}(1);
      if fidx<expstopfile,
         stoptimes(fidx)=rvalid{fidx}(resplens(fidx));
      elseif fidx==expstopfile,
         stoptimes(fidx)=rvalid{fidx}(expstopframe);
      else
         stoptimes(fidx)=0;
      end
   end
end

% figure out file and start/stop idexs
% for the time being, assume fitfrac > 0 
if fitfile==predfile,
   fitstopframe=min([predstartframe-1 fitstartframe+MAXPREDLEN]);
else
   predstartframe=1;
   fitstopframe=min([fitstartframe+MAXPREDLEN resplens(fitfile)]);
end

predstartframe=rvalid{predfile}(predstartframe);
predstopframe=rvalid{predfile}(predstopframe);
fitstartframe=rvalid{fitfile}(fitstartframe);
if fitstopframe>0,
   fitstopframe=rvalid{fitfile}(fitstopframe);
end

if fitfrac<=0,
   % fit to first MAXPREDLEN frames of exploratory
   % data... usually that means all of it for review.
   fprintf('Using all (up to %d frames) exploratory data for fit.\n',...
           MAXPREDLEN);
   
   expcumlen=cumsum(stoptimes-starttimes+1);
   tfitstopframe=min([expcumlen(end) MAXPREDLEN]);
   fitfilemax=min(find(tfitstopframe<=expcumlen));
   
   fitlen=0;
   for fitidx=1:fitfilemax,
      fitfile(fitidx)=fitidx;
      fitstartframe(fitidx)=starttimes(fitidx);
      fitstopframe(fitidx)=stoptimes(fitidx);
      fitlen=fitlen+stoptimes(fitidx)-starttimes(fitidx)+1;
      if fitlen>tfitstopframe,
         fitstopframe(fitidx)=fitstopframe(fitidx)-fitlen+tfitstopframe;
      end
   end
end

%maxexpfidx=max(find(stoptimes-starttimes>0));
maxexpfidx=length(starttimes);
times(1).fileidx=repmat((1:maxexpfidx)',[1 attcount]);
times(1).start=starttimes;
times(1).stop=stoptimes;
times(1).cat='exp';

times(2).fileidx=fitfile;
times(2).start=fitstartframe;
times(2).stop=fitstopframe;
times(2).cat='fit';

times(3).fileidx=predfile;
times(3).start=predstartframe;
times(3).stop=predstopframe;
times(3).cat='pred';


