% function times=xcfilefracs(params);
%
% fields for params (ie, params.<fieldname>)
% "ERROR" means required
%
% created SVD 3/03 - ripped of cellfiletimes
% 
function times=xcfilefracsold(params);

MAXPREDLEN=7000; % plenty data for fitting regularization/threshold parms?

dbopen;

respfiles=getparm(params,'respfiles','ERROR');
expfrac=getparm(params,'expfrac',0);
fitfrac=getparm(params,'fitfrac',0);
predfrac=getparm(params,'predfrac',0.1);
decorrspace=getparm(params,'decorrspace',2);
resploadcmd=getparm(params,'resploadcmd','respload');
resploadparms=getparm(params,'resploadparms',{'',1,1,1});
respfiltercmd=getparm(params,'respfiltercmd','');
respfilterparms=getparm(params,'respfilterparms',[]);
resampcount=getparm(params,'resampcount',20);

if length(respfiles)==0,
   times=[];
   return
end

% figure out separate start and stop times for each attentional state!
rvalid={};
filecount=length(respfiles);
for fidx=1:filecount,
   resp=feval(resploadcmd,respfiles{fidx},resploadparms{:});
   % filter response (eg resample, pick attentional state, etc)
   if ~isempty(respfiltercmd),
      resp=feval(respfiltercmd,resp,respfilterparms{:});
   end
   rsize=size(resp);
   
   if rsize(1)==1,
      resp=resp(:);
      rsize=size(resp);
   end
   
   % adjust attcount to actual number of attention states
   if fidx==1,
      attcount=size(resp,3);  % dim 3 is number attentional states.
      resplens=zeros(filecount,attcount);
   end
   for attidx=1:attcount,
      rvalid{fidx,attidx}=find(~isnan(squeeze(resp(:,1,attidx))));
      resplens(fidx,attidx)=length(rvalid{fidx,attidx});
   end
end

starttimes=ones(filecount,attcount);
stoptimes=zeros(filecount,attcount);
cumresplens=cumsum(resplens,1);
totlen=sum(resplens,1);

if predfrac>0,
   predstartframe=round(totlen*(1-predfrac));
   % min pred is 200 frames
   if totlen*predfrac<200 & totlen>400,
      predstartframe=totlen-199;
   end
   
   predfile=sum(cumresplens<predstartframe,1)+1;
   % make sure that last file is getting used in pred file
   for attidx=1:attcount,
      if predfile(attidx)<filecount,
         predfile(attidx)=filecount;
         predstartframe(attidx)=sum(resplens(1:(filecount-1)))+1;
      end
   end
else
   predstartframe=totlen;
   predstopframe=totlen;  
   predfile=filecount;
end

predstopframe=zeros(1,attcount);
for attidx=1:attcount,
   predstartframe(attidx)=predstartframe(attidx) - ...
       sum(resplens(1:(predfile(attidx)-1),attidx));
   predstopframe(attidx)=min([predstartframe(attidx)+MAXPREDLEN ...
                    resplens(predfile(attidx),attidx)]);
end


expstopfile=zeros(1,attcount);
expstopframe=zeros(1,attcount);
fitfile=zeros(1,attcount);
fitstartframe=zeros(1,attcount);
fitstopframe=zeros(1,attcount);
for attidx=1:attcount,
   tpredstart=predstartframe(attidx)+...
       sum(resplens(1:predfile(attidx)-1,attidx));
   fitstartframe(attidx)=round(tpredstart-totlen(attidx)*fitfrac);
   fitfile(attidx)=sum(cumresplens(:,attidx)<fitstartframe(attidx),1)+1;
   if fitfile>filecount,
      fitfile=filecount;
   end
   
   if expfrac<=0,
      texpfrac=1-fitfrac-(predstopframe(attidx)- ...
                          predstartframe(attidx))./totlen(attidx);
      expstopframe(attidx)=round(totlen(attidx)*texpfrac)-1;
   else
      expstopframe(attidx)=round(totlen(attidx)*expfrac)-1;
   end
   expstopfile(attidx)=sum(cumresplens(:,attidx)<expstopframe(attidx),1)+1;
   if fitfile(attidx)<filecount-1,
      fitfile(attidx)=filecount-1;
      fitstartframe(attidx)=cumresplens(filecount-2,attidx)+1;
   elseif cumresplens(fitfile(attidx),attidx)-fitstartframe(attidx) < ...
         floor(totlen(attidx)*fitfrac/2) & fitfile(attidx)<filecount,
      fitfile(attidx)=fitfile(attidx)+1;
      fitstartframe(attidx)=cumresplens(fitfile(attidx)-1)+1;
   end
   
   expstopframe(attidx)=expstopframe(attidx) - ...
       sum(resplens(1:(expstopfile(attidx)-1),attidx));
   if expstopframe(attidx)<resampcount & expstopfile(attidx)>1,
      expstopfile(attidx)=expstopfile(attidx)-1;
      expstopframe(attidx)=resplens(expstopfile(attidx),attidx);
   end
   fitstartframe(attidx)=fitstartframe(attidx) - ...
       sum(resplens(1:(fitfile(attidx)-1),attidx));
   
   for fidx=1:fitfile(attidx),
      %keyboard
      if resplens(fidx,attidx)>0,
         starttimes(fidx,attidx)=rvalid{fidx,attidx}(1);
         if fidx<expstopfile(attidx),
            stoptimes(fidx,attidx)=rvalid{fidx,attidx}(resplens(fidx,attidx));
         elseif fidx==expstopfile(attidx),
            stoptimes(fidx,attidx)=rvalid{fidx,attidx}(expstopframe(attidx));
         else
            stoptimes(fidx,attidx)=0;
         end
      end
   end
   
   % figure out file and start/stop idexs
   % for the time being, assume fitfrac > 0 
   %  fitfrac=fitfrac;
   if fitfile(attidx)==predfile(attidx),
      fitstopframe(attidx)=min([predstartframe(attidx)-1 ...
                    fitstartframe(attidx)+MAXPREDLEN]);
   else
      predstartframe(attidx)=1;
      fitstopframe(attidx)=min([fitstartframe(attidx)+MAXPREDLEN ...
                    resplens(fitfile,attidx)]);
   end
   
   predstartframe(attidx)=rvalid{predfile(attidx),attidx}(predstartframe(attidx));
   predstopframe(attidx)=rvalid{predfile(attidx),attidx}(predstopframe(attidx));
   fitstartframe(attidx)=rvalid{fitfile(attidx),attidx}(fitstartframe(attidx));
   if fitstopframe(attidx)>0,
      fitstopframe(attidx)=rvalid{fitfile(attidx),attidx}(fitstopframe(attidx));
   
   end
   
   if fitfrac<=0,
      
      % 0 DISABLES NEWER (BETTER? OVERFIT?) EXP
      if 1,
      % fit to first MAXPREDLEN frames of exploratory
      % data... usually that means all of it for review.
      fprintf('Using all (up to %d frames) exploratory data for fit.\n',...
              MAXPREDLEN);
      
      expcumlen=cumsum(stoptimes(:,attidx)-starttimes(:,attidx)+1);
      tfitstopframe=min([expcumlen(end) MAXPREDLEN]);
      fitfilemax=min(find(tfitstopframe<=expcumlen));
      
      fitlen=0;
      for fitidx=1:fitfilemax,
         fitfile(fitidx,attidx)=fitidx;
         fitstartframe(fitidx,attidx)=starttimes(fitidx,attidx);
         fitstopframe(fitidx,attidx)=stoptimes(fitidx,attidx);
         fitlen=fitlen+stoptimes(fitidx,attidx)-starttimes(fitidx,attidx)+1;
         if fitlen>tfitstopframe,
            fitstopframe(fitidx,attidx)=fitstopframe(fitidx,attidx)-...
                fitlen+tfitstopframe;
         end
      end
      else
      
      % old crap sort of hacky for picking exploratory data .. but
      % does it actually work better??? Or is it a case of
      % overfitting???
      % find longest exploratory file for doing fitting.
      disp('Using longest exploratory file for fit.');
      
      fitfile(1,attidx)=1;
      maxexplen=0;
      for fidx=1:filecount,
         if stoptimes(fidx,attidx)-starttimes(fidx,attidx)>maxexplen,
            fitfile(1,attcount)=fidx;
            maxexplen=stoptimes(fidx,attidx)-starttimes(fidx,attidx);
         end
      end
      fprintf('Max exp run len=%d (%d).\n',maxexplen,fitfile(1,attcount));
      
      fitstartframe(attidx)=starttimes(fitfile(attidx),attidx);
      fitstopframe(attidx)=stoptimes(fitfile(attidx),attidx);
      if fitstopframe(attidx)-fitstartframe(attidx)>=MAXPREDLEN,
         fitstopframe(attidx)=fitstartframe(attidx)+MAXPREDLEN-1;
      end
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


