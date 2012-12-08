% lnrcore.m
%
% dumb linear reconstruction using a single time bin (smoothed by
% params.smoothparm)
%
% created SVD 2011-11-04 -- ripped off mlrcore.
%
disp([mfilename,':']);

% silly kludge where "stim" and "resp" are reversed from standard
% STRF framework.  reverse back with new variables.
s0=resp;
s0(isnan(s0))=0;

r0=stim;
cellcount=size(r0,2);
resplen=size(r0,1);

% smooth response to generate appropriate integration time
fprintf('Response window=%d bins\n',params.smoothwin);
smfilt=[zeros(params.smoothwin-1,1);ones(params.smoothwin,1)];
r=conv2(r0,smfilt,'same');

% remove silent periods if selected
EXCLUDESILENCE=1;
if EXCLUDESILENCE,
   ff=double(s0==0);
   smfilt=ones(params.maxlag(2)-params.maxlag(1)+1,1);
   ff=conv2(ff,smfilt,'same');
   ff=find(ff);
   r(ff,:)=nan;
end

if isfield(params,'adddelay') && params.adddelay,
   % two responses, one offset by a smoothbin
   fprintf('adding 2nd resp for each cell shifted back by %d bins\n',...
           params.smoothwin);
   r=permute(r,[1 3 2]);
   rshift=shift(r,params.smoothwin);
   rshift(1:params.smoothwin,:,:)=nan;
   r=cat(2,r,rshift);
   r=reshape(r,size(r,1),size(r,2).*size(r,3));
   cellcount=size(r,2);
end

ff=find(sum(isnan(r),2)==0);

mS=mean(s0(ff));
s0=s0-mS;

plagcount=params.maxlag(2)-params.maxlag(1)+1;
s=ones(resplen,plagcount).*nan;
for ii=1:plagcount,
   os=ii-params.maxlag(2)-1;
   if os<0
      s(1-os:resplen,ii)=s0(1:resplen+os);
   else
      s(1:resplen-os,ii)=s0(os+1:resplen);
   end
end

ff=find(~isnan(r(:,1)) & ~isnan(r(:,end)));

[H,BRAC]=revrecCore(s(ff,:),r(ff,:));

if 0,
   siteid=strsep(params.cellid,'-');
   siteid=siteid{1};
   
   figure;
   subplot(6,1,1);
   rasterfs=params.stimloadparms{1}.rasterfs;
   plot((1:600)./rasterfs,s0(1:600),'k','LineWidth',2);
   title(sprintf('%s stimulus envelope',siteid));
   
   subplot(6,1,2);
   rasterfs=params.stimloadparms{1}.rasterfs;
   plot((1:600)./rasterfs,r0(1:600,:),'LineWidth',2);
   xlabel('time (s)');
   ylabel('spike count');
   if length(params.stimloadparms{1}.channel)>1,
      title(sprintf('%s cells (%d,%d) (%d,%d)',...
                    siteid,params.stimloadparms{1}.channel(1),...
                    params.stimloadparms{1}.unit(1),...
                    params.stimloadparms{1}.channel(2),...
                    params.stimloadparms{1}.unit(2)));
   else
      title(sprintf('Cell %s (lag 1, lag 2)',params.cellid));
   end
   
   dotstr={'b.','g.'};
   linstr={'b--','g--'};
   for cellidx=1:cellcount
      uff=shuffle(ff);
      uff=ff(1:1250);
      
      subplot(3,2,2+cellidx);
      plot(r(uff,cellidx),s0(uff-1),dotstr{cellidx});
      hold on
      p=polyfit(r(ff,cellidx),s0(ff-1),1);
      x1=max(r(:,cellidx));
      plot([0 x1],[0 x1].*p(1)+p(2),linstr{cellidx},'LineWidth',2);
      hold off
      xlabel(sprintf('cell %d response',cellidx));
      ylabel('stimulus');
      
      subplot(3,2,4+cellidx);
      lags=(-params.maxlag(2)):(-params.maxlag(1));
      plot(lags,H(cellidx,:),linstr{cellidx}(1:2),'LineWidth',2);
      title('linear recon filter');
      axis([lags(1) lags(end) ...
            min(H(:))-std(H(:))./2 max(H(:))+std(H(:))./2]);
   end

   fullpage portrait
   
   keyboard
end



% don't need a separate fit routine. just make an "STRF" here
disp('generating STRF');
strf=struct;
strf.h=H;
strf.mS=zeros(size(r0,2),1);
strf.name='lnR';
strf.architecture='lnR';
strf.parms.adddelay=params.adddelay;
strf.parms.smoothwin=params.smoothwin;
strf.parms.maxlag=params.maxlag;
strf.nlparms=[];
strf.nltype='none';

clear s0 r r0 tresp ff


% strf structure:
% strf( ).h    :  linear kernel
% strf( ).mS   :  mean stimulus used in fit for subtracting before pred
% strf( ).name :  info about what this is
% strf( ).parms.sfsidx   : best fit idx
% strf( ).parms.sigidx   : best fit idx
% strf( ).parms.kernfmt  : kernfmt for showkern
% strf( ).parms.iconside : stim size info for showkern
% strf( ).nlparms : parms for output nl
% strf( ).nltype  : string desc of nl
% strf( ).xc      : results of fit xc test.

