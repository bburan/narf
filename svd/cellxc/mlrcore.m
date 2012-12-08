% mlrcore.m
%
% direct probability reconstruction using a single time bin
%
% created SVD 2011-10-28
%
disp([mfilename,':']);

% silly kludge where "stim" and "resp" are reversed from standard
% STRF framework.  reverse back with new variables.
s0=resp;
s0(isnan(s0))=0;

r0=stim;
rmin=nanmin(r0);
cellcount=size(r0,2);
for ii=1:cellcount,
   r0(:,ii)=r0(:,ii)-rmin(ii);
end
resplen=size(r0,1);

% smooth response to generate appropriate integration time
fprintf('Response window=%d bins\n',params.smoothwin);
smfilt=[zeros(params.smoothwin-1,1);ones(params.smoothwin,1)];
r=conv2(r0,smfilt,'same');

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

   
% remove silent periods if selected
EXCLUDESILENCE=1;
if EXCLUDESILENCE,
   ff=double(s0==0);
   smfilt=ones(params.maxlag(2)-params.maxlag(1)+1,1);
   ff=conv2(ff,smfilt,'same');
   ff=find(ff);
   r(ff,:)=nan;
end

disp('binning response');
if cellcount==2,
    maxr=5;
elseif cellcount==3,
    maxr=4;
end
rbins=zeros(maxr+1,cellcount);
ff=find(~isnan(r(:,1)));
rbinned=zeros(size(r));
for cc=1:cellcount,
   ru=unique(r(ff,cc));
   n=hist(r(ff,cc),ru);
   ridx=1;
   rstep=sum(n)./(maxr);
   for ri=1:maxr,
      rcount=min(find(cumsum(n(ridx:end))>=rstep))+ridx-1;
      if ~isempty(rcount),
         if rcount==length(ru),
            rbins(ri+1,cc)=ru(end);
         else
            rbins(ri+1,cc)=ru(rcount+1);
         end
         ridx=rcount+1;
         rstep=sum(n(ridx:end))./(maxr-ri+1);
      end
   end
   rbins(end,cc)=inf;
   for ii=1:maxr,
      rbinned(ff(find(r(ff,cc)>=rbins(ii,cc) & ...
                      r(ff,cc)<rbins(ii+1,cc))),cc)=ii;
   end
end

if min(min(rbinned(ff,:)))==0,
   disp('0 binned value!');
   keyboard
end


% figure out how to bin the stimulus
disp('binning stimulus');
if params.mlindep,
    sbincount=24;
else
    sbincount=20;
end
sbins=linspace(min(s0),max(s0)+nanstd(s0)./100,sbincount+1)';
sbincenters=sbins(1:(end-1))+diff(sbins)./2;

sbinned=zeros(size(s0));
for ii=1:sbincount,
   sbinned(find(s0>=sbins(ii) & s0<sbins(ii+1)))=ii;
end

plagcount=params.maxlag(2)-params.maxlag(1)+1;
if params.mlindep,
   disp('calculating ML stim according to independent spiking probabilities');
   shist0=zeros(sbincount,plagcount,maxr,cellcount);
   rcount=zeros(maxr,cellcount);
   for c1=1:cellcount,
      for r1=1:maxr,
         fprintf('cell %d rbinned %d\n',c1,r1);
         rii=find(rbinned(:,c1)==r1);
         rcount(r1,c1)=length(rii);
         for ii=1:length(rii),
            sidx=rii(ii)+(-params.maxlag(2):-params.maxlag(1));
            for jj=1:length(sidx),
               shist0(sbinned(sidx(jj)),jj,r1,c1)=...
                   shist0(sbinned(sidx(jj)),jj,r1,c1)+1;
            end
         end
         if rcount(r1,c1)==0 && r1>1,
             shist0(:,:,r1,c1)=shist0(:,:,r1-1,c1);
         end
             
      end
   end
   
   % p(s|ri)= count(s=sj,r=ri)./count(r=ri)
   shist0=shist0./repmat(sum(shist0,1),[sbincount 1 1 1]);
   shist0(isnan(shist0))=0;
   
   % if indep then p(s|r1,r2)=p(s|r1)*p(s|r2)
   smax=zeros(plagcount,maxr,maxr);
   shist=zeros(sbincount,plagcount,maxr,maxr);
   for r1=1:maxr,
      for r2=1:maxr,
          
          sh1=gsmooth(shist0(:,:,r1,1),[1.2 0.25]);
          sh2=gsmooth(shist0(:,:,r2,2),[1.2 0.25]);
          sr1r2=sh1.*sh2;
          shist(:,:,r1,r2)=sr1r2;
          %sr1r2=shist(:,:,r1,1).*shist(:,:,r2,2);
          for jj=1:size(shist0,2),
              smax(jj,r1,r2)=...
                  sbincenters(min(find(sr1r2(:,jj)==max(sr1r2(:,jj)))));
          end
      end
   end
   
elseif cellcount==2,
   
   disp('calculating ML stim according to joint spiking probabilities');
   shist=zeros(sbincount,plagcount,maxr,maxr);
   smax=zeros(plagcount,maxr,maxr);
   rcount=zeros(maxr,maxr);
   for r1=1:maxr,
      for r2=1:maxr,
         rii=find(rbinned(:,1)==r1 & rbinned(:,2)==r2);
         rcount(r1,r2)=length(rii);
         for ii=1:length(rii),
            sidx=rii(ii)+(-params.maxlag(2):-params.maxlag(1));
            for jj=1:length(sidx),
               shist(sbinned(sidx(jj)),jj,r1,r2)=...
                shist(sbinned(sidx(jj)),jj,r1,r2)+1;
            end
         end
         
         shist(:,:,r1,r2)=gsmooth(shist(:,:,r1,r2),[1.2 0.25]);
         sh=shist(:,:,r1,r2);
         %figure(1);clf;
         %subplot(2,1,1);imagesc(shist(:,:,r1,r2));
         %subplot(2,1,2);imagesc(sh);
         %[r1 r2],pause(0.1);
         
         for jj=1:size(shist,2),
             %smax(jj,r1,r2)=sbincenters(min(find(shist(:,jj,r1,r2)==...
             %                                   max(shist(:,jj,r1,r2)))));
            smax(jj,r1,r2)=sbincenters(min(find(sh(:,jj)==max(sh(:,jj)))));
         end
      end
   end
   
   %keyboard
   
elseif cellcount==3,
   disp('calculating ML stim according to joint spiking probabilities');
   shist=zeros(sbincount,plagcount,maxr,maxr);
   smax=zeros(plagcount,maxr,maxr,maxr);
   rcount=zeros(maxr,maxr,maxr);
   for r1=1:maxr,
      for r2=1:maxr,
          for r3=1:maxr,
              shist=zeros(sbincount,plagcount);
              rii=find(rbinned(:,1)==r1 & rbinned(:,2)==r2 & rbinned(:,3)==r3);
              rcount(r1,r2,r3)=length(rii);
              for ii=1:length(rii),
                  sidx=rii(ii)+(-params.maxlag(2):-params.maxlag(1));
                  for jj=1:length(sidx),
                      shist(sbinned(sidx(jj)),jj)=...
                          shist(sbinned(sidx(jj)),jj)+1;
                  end
              end
              for jj=1:size(shist,2),
                  smax(jj,r1,r2)=sbincenters(min(find(shist(:,jj)==...
                                                      max(shist(:,jj)))));
              end
              figure(1);clf;imagesc(shist);[r1 r2 r3],pause(0.1);
          end
      end
   end
   disp('triplet ML recon not fully debugged');
   keyboard
end

% don't need a separate fit routine. just make an "STRF" here
disp('generating STRF');
strf=struct;
strf.h=smax;
strf.architecture='ML';
strf.name='ML';
strf.mS=rmin';
strf.parms.adddelay=params.adddelay;
strf.parms.sbins=sbins;
strf.parms.rbins=rbins;
strf.parms.smoothwin=params.smoothwin;
strf.parms.maxlag=params.maxlag;
strf.nlparms=[];
strf.nltype='none';

% strf structure (linear model):
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
      uff=[];
      for r1=1:maxr,
         sff=ff(find(rbinned(ff,cellidx)==r1));
         sff=shuffle(ff);
         sff=sff(1:min(250,length(sff)));
         uff=[uff;sff];
      end
      
      subplot(3,2,2+cellidx);
      plot(r(uff,cellidx),s0(uff-1),dotstr{cellidx});
      hold on
      p=polyfit(r(ff,cellidx),s0(ff-1),1);
      x1=max(r(:,cellidx));
      plot([0 x1],[0 x1].*p(1)+p(2),linstr{cellidx});
      hold off
      xlabel(sprintf('cell %d response',cellidx));
      ylabel('stimulus');
      
   end
   subplot(3,2,5)
   imagesc(rbins(1:maxr,2),rbins(1:maxr,1),rcount');
   axis xy
   axis square
   xlabel('cell 1 spike count');
   ylabel('cell 2 spike count');
   colorbar
   colormap(1-gray);
   title('response distribution');

   fullpage portrait
   
   
   
   figure;
   sbincenters=sbins(1:(end-1))+(sbins(2)-sbins(1))./2;
   lags=(-params.maxlag(2)):(-params.maxlag(1));
   for r1=1:maxr,
      for r2=1:maxr,
         subplot(maxr,maxr,r1+(maxr-r2).*maxr);
         imagesc(lags,sbincenters,...
                 shist(:,:,r1,r2));
         hold on
         plot(lags,smax(:,r1,r2),'k','LineWidth',2);
         hold off
         axis xy;
         %axis([params.maxlag 0 max(sbins)]);
         title(sprintf('(r1,r2)=(%d,%d)',rbins(r1,1),rbins(r2,2)));
      end
   end
   fullpage landscape
   
   keyboard
end


clear s0 r r0 tresp ff


