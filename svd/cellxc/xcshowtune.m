% function [tunedata,fitdata]=xcshowtune(kernfiles/cellid,batch)
%
function [tunedata,fitdata]=xcshowtune(cellid,batch)

dbopen;

if nargout>0,
   tunedata=[];
   fitdata=[];
end

if iscell(cellid),
   kernfiles=cellid;
   batchcount=length(resfiles);
elseif strcmp(cellid(end-3:end),'.mat') | ...
      strcmp(cellid(end-2:end),'.gz'),
   kernfiles={cellid};
   batchcount=1;
else
   goodbatch=zeros(1,length(batch));
   batchcount=0;
   resfiles={};
   for ii=1:length(batch),
      sql=['SELECT * from sRunData WHERE cellid="',cellid,'"',...
           ' AND batch=',num2str(batch(ii))];
      trd=mysql(sql);
      if ~isempty(trd),
         goodbatch(ii)=1;
         batchcount=batchcount+1;
         kernfiles{batchcount}=[trd.respath,trd.kernfile];
         rundata(batchcount)=trd;
      end
   end
   batch=batch(find(goodbatch));
   batchcount=length(batch);
end

fprintf('xcshowtune.m:  %d files:\n',batchcount);

for batchidx=1:batchcount
   fprintf(' %s\n',kernfiles{batchidx});
   if strcmp(kernfiles{batchidx}(end-2:end),'.gz'),
      r=zload(kernfiles{batchidx});
   else
      r=load(kernfiles{batchidx});
   end
   nlidx=2;
   if nlidx>length(r.strf),
      nlidx=1;
   end
   
   r.fitdata.batch=batch(batchidx);
   if batchidx==1,
      tunedata=r.tunedata;
      fitdata=r.fitdata;
      strf=r.strf(nlidx,:,1);
   else
      tunedata=cat(2,tunedata,r.tunedata);
      fitdata=cat(2,fitdata,r.fitdata);
      strf=cat(1,strf,r.strf(nlidx,:,1));
      
   end
end

obincount=length(tunedata(1).obins);

% plot tuning curves and associated fits
figure
for batchidx=1:batchcount
   
   bootcount=size(tunedata,1);
   
   subplot(4,batchcount,batchidx);
   mor=cat(3,tunedata(:,batchidx).por);
   eor=std(mor(:,1,:),1,3) .* sqrt(bootcount-1);
   mor=mean(mor(:,1,:),3);
   errorbar(tunedata(1,batchidx).obins,mor,eor);
   
   hold on
   plot(tunedata(1,batchidx).obins,...
        flatcos(fitdata(batchidx).orfit,tunedata(1,batchidx).obins),'r--');
   hold off
   
   title(sprintf('%s bat %d or',cellid,batch(batchidx)));
   
   subplot(4,batchcount,batchidx+batchcount);
   msf=cat(3,tunedata(:,batchidx).psf) ;
   esf=std(msf(:,1,:),1,3) .* sqrt(bootcount-1);
   msf=mean(msf(:,1,:),3);
   errorbar(tunedata(1,batchidx).sfbins,msf,esf);
   title(sprintf('%s bat %d sf',cellid,batch(batchidx)));
   
   hold on
   plot(tunedata(1,batchidx).sfbins,...
        gauss1(fitdata(batchidx).sffit,log2(tunedata(1,batchidx).sfbins)),'r--');
   hold off
   
   subplot(4,batchcount,batchidx+batchcount*2);
   mtime=cat(3,tunedata(:,batchidx).ptime) ;
   etime=std(mtime(:,1,:),1,3) .* sqrt(bootcount-1);
   mtime=mean(mtime(:,1,:),3);
   errorbar(tunedata(1,batchidx).tbins,mtime,etime);
   title(sprintf('%s bat %d time',cellid,batch(batchidx)));
   
   hold on
   plot(tunedata(1,batchidx).tbins,...
        expdecaydiff(fitdata(batchidx).timefit,tunedata(1,batchidx).tbins),'r--');
   hold off
end


% plot surfaces representing tuning curves for each jackknife
figure;
for batchidx=1:batchcount
   
   subplot(4,batchcount,batchidx);
   mor=cat(3,tunedata(:,batchidx).por);
   imagesc(squeeze(mor(:,1,:)));
   title(sprintf('%s bat %d or',cellid,batch(batchidx)));
   
   subplot(4,batchcount,batchidx+batchcount);
   msf=cat(3,tunedata(:,batchidx).psf) ;
   imagesc(squeeze(msf(:,1,:)));
   title(sprintf('%s bat %d sf',cellid,batch(batchidx)));
   
   subplot(4,batchcount,batchidx+batchcount*2);
   mtime=cat(3,tunedata(:,batchidx).ptime) ;
   imagesc(squeeze(mtime(:,1,:)));
   title(sprintf('%s bat %d time',cellid,batch(batchidx)));
   
end


%plot the new statistics
figure

spacecount=size(strf(1).h,1);
tcount=size(strf(1).h,2);

ttimepos=zeros(tcount,bootcount,batchcount);
ttimeneg=zeros(tcount,bootcount,batchcount);
ttimesep=zeros(tcount,bootcount,batchcount);
tspacemaxpos=zeros(spacecount,bootcount,batchcount);
tspacemaxneg=zeros(spacecount,bootcount,batchcount);
tspaceall=zeros(spacecount,bootcount,batchcount);
tspacesep=zeros(spacecount,bootcount,batchcount);
teigvals=zeros(tcount,bootcount,batchcount);

for batchidx=1:batchcount,
   for bootidx=1:bootcount,
      h=strf(batchidx,bootidx).h;
      
      ttimepos(:,bootidx,batchidx)=sum(h.*(h>0),1)';
      ttimeneg(:,bootidx,batchidx)=sum(h.*(h<0),1)';
      
      tmax=min(find(ttimepos(:,bootidx,batchidx)==...
                    max(ttimepos(:,bootidx,batchidx))));
      tspacemaxpos(:,bootidx,batchidx)=h(:,tmax);
      tmin=min(find(ttimeneg(:,bootidx,batchidx)==...
                    min(ttimeneg(:,bootidx,batchidx))));
      tspacemaxneg(:,bootidx,batchidx)=h(:,tmin);
      tspaceall(:,bootidx,batchidx)=sum(h,2);
      
      [u,s,v]=svd(h);
      if tspaceall(:,bootidx,batchidx)'*u(:,1) > 0;
         ttimesep(:,bootidx,batchidx)=v(:,1);
         tspacesep(:,bootidx,batchidx)=u(:,1);
      else
         ttimesep(:,bootidx,batchidx)=-v(:,1);
         tspacesep(:,bootidx,batchidx)=-u(:,1);
      end
      teigvals(:,bootidx,batchidx)=diag(s);
   end
end

% get average spatial/temporal kernel for each class
timepos=squeeze(mean(ttimepos,2));
etimepos=squeeze(std(ttimepos,1,2)) .* sqrt(bootcount-1);
timeneg=squeeze(mean(ttimeneg,2));
etimeneg=squeeze(std(ttimeneg,1,2)) .* sqrt(bootcount-1);
spacemaxpos=squeeze(mean(tspacemaxpos,2));
espacemaxpos=squeeze(std(tspacemaxpos,1,2)) .* sqrt(bootcount-1);
spacemaxneg=squeeze(mean(tspacemaxneg,2));
espacemaxneg=squeeze(std(tspacemaxneg,1,2)) .* sqrt(bootcount-1);
spaceall=squeeze(mean(tspaceall,2));
espaceall=squeeze(std(tspaceall,1,2)) .* sqrt(bootcount-1);

timesep=squeeze(mean(ttimesep,2));
etimesep=squeeze(std(ttimesep,1,2)) .* sqrt(bootcount-1);
spacesep=squeeze(mean(tspacesep,2));
espacesep=squeeze(std(tspacesep,1,2)) .* sqrt(bootcount-1);

% compute s-t separability and std err
tsepfrac=teigvals(1,:,:).^2./sum(teigvals.^2,1);
sepfrac=squeeze(mean(tsepfrac,2));
esepfrac=squeeze(std(tsepfrac,1,2)) .* sqrt(bootcount-1);

% computer corr between kernels
tpc=zeros(bootcount,batchcount,batchcount);
tnc=zeros(bootcount,batchcount,batchcount);
tsc=zeros(bootcount,batchcount,batchcount);
spc=zeros(bootcount,batchcount,batchcount);
spposc=zeros(bootcount,batchcount,batchcount);
spnegc=zeros(bootcount,batchcount,batchcount);
snc=zeros(bootcount,batchcount,batchcount);
ssc=zeros(bootcount,batchcount,batchcount);
ssposc=zeros(bootcount,batchcount,batchcount);
ssnegc=zeros(bootcount,batchcount,batchcount);
sac=zeros(bootcount,batchcount,batchcount);

for b1=1:batchcount,
   for b2=b1:batchcount,
      for bootidx=1:25,
         
         % corr for temp stuff
         tpc(bootidx,b1,b2)=xcorr(ttimepos(:,bootidx,b1),...
                                  ttimepos(:,bootidx,b2),0,'coeff');
         tpc(bootidx,b2,b1)=tpc(bootidx,b1,b2);
         tnc(bootidx,b1,b2)=xcorr(ttimeneg(:,bootidx,b1),...
                                  ttimeneg(:,bootidx,b2),0,'coeff');
         tnc(bootidx,b2,b1)=tnc(bootidx,b1,b2);
         tsc(bootidx,b1,b2)=xcorr(ttimesep(:,bootidx,b1),...
                                  ttimesep(:,bootidx,b2),0,'coeff');
         tsc(bootidx,b2,b1)=tsc(bootidx,b1,b2);
         
         % corr for space stuff
         ts1=tspacemaxpos(:,bootidx,b1);
         ts2=tspacemaxpos(:,bootidx,b2);
         if sum(abs(ts1)) & sum(abs(ts2)),
            spc(bootidx,b1,b2)=xcorr(ts1,ts2,0,'coeff');
            spc(bootidx,b2,b1)=spc(bootidx,b1,b2);
         end
         if sum(abs(ts1.*(ts1>0))) & sum(abs(ts2.*(ts2>0))),
            spposc(bootidx,b1,b2)=xcorr(ts1.*(ts1>0),ts2.*(ts2>0),0,'coeff');
            spposc(bootidx,b2,b1)=spposc(bootidx,b1,b2);
         end
         if sum(abs(ts1.*(ts1<0))) & sum(abs(ts2.*(ts2<0))),
            spnegc(bootidx,b1,b2)=xcorr(ts1.*(ts1<0),ts2.*(ts2<0),0,'coeff');
            spnegc(bootidx,b2,b1)=spnegc(bootidx,b1,b2);
         end
         
         
         
         if sum(abs(tspacemaxneg(:,bootidx,b1))) & ...
               sum(abs(tspacemaxneg(:,bootidx,b2))),
            snc(bootidx,b1,b2)=xcorr(tspacemaxneg(:,bootidx,b1),...
                                     tspacemaxneg(:,bootidx,b2),0,'coeff');
            snc(bootidx,b2,b1)=snc(bootidx,b1,b2);
         end
         
         if sum(abs(tspaceall(:,bootidx,b1))) & ...
               sum(abs(tspaceall(:,bootidx,b2))),
            sac(bootidx,b1,b2)=xcorr(tspaceall(:,bootidx,b1),...
                                     tspaceall(:,bootidx,b2),0,'coeff');
            sac(bootidx,b2,b1)=sac(bootidx,b1,b2);
         end
         
         ts1=tspacesep(:,bootidx,b1);
         ts2=tspacesep(:,bootidx,b2);
         if sum(abs(ts1)) & sum(abs(ts2)),
            ssc(bootidx,b1,b2)=xcorr(ts1,ts2,0,'coeff');
            ssc(bootidx,b2,b1)=ssc(bootidx,b1,b2);
         end
         if sum(abs(ts1.*(ts1>0))) & sum(abs(ts2.*(ts2>0))),
            ssposc(bootidx,b1,b2)=xcorr(ts1.*(ts1>0),ts2.*(ts2>0),0,'coeff');
            ssposc(bootidx,b2,b1)=ssposc(bootidx,b1,b2);
         end
         if sum(abs(ts1.*(ts1<0))) & sum(abs(ts2.*(ts2<0))),
            ssnegc(bootidx,b1,b2)=xcorr(ts1.*(ts1<0),ts2.*(ts2<0),0,'coeff');
            ssnegc(bootidx,b2,b1)=ssnegc(bootidx,b1,b2);
         end
         
      end
   end
end
etpc=squeeze(std(tpc,1,1)) .* sqrt(bootcount-1);
tpc=squeeze(mean(tpc,1));
etnc=squeeze(std(tnc,1,1)) .* sqrt(bootcount-1);
tnc=squeeze(mean(tnc,1));
etsc=squeeze(std(tsc,1,1)) .* sqrt(bootcount-1);
tsc=squeeze(mean(tsc,1));

espc=squeeze(std(spc,1,1)) .* sqrt(bootcount-1);
spc=squeeze(mean(spc,1));
espposc=squeeze(std(spposc,1,1)) .* sqrt(bootcount-1);
spposc=squeeze(mean(spposc,1));
espnegc=squeeze(std(spnegc,1,1)) .* sqrt(bootcount-1);
spnegc=squeeze(mean(spnegc,1));
esnc=squeeze(std(snc,1,1)) .* sqrt(bootcount-1);
snc=squeeze(mean(snc,1));
esac=squeeze(std(sac,1,1)) .* sqrt(bootcount-1);
sac=squeeze(mean(sac,1));
essc=squeeze(std(ssc,1,1)) .* sqrt(bootcount-1);
ssc=squeeze(mean(ssc,1));
essposc=squeeze(std(ssposc,1,1)) .* sqrt(bootcount-1);
ssposc=squeeze(mean(ssposc,1));
essnegc=squeeze(std(ssnegc,1,1)) .* sqrt(bootcount-1);
ssnegc=squeeze(mean(ssnegc,1));

labels={'tpc','tnc','tsc','spc','spc+','spc-',...
        'snc','sac','ssc','ssc+','ssc-'};

for b1=1:batchcount,
   for b2=b1+1:batchcount,
      subplot(batchcount-1,batchcount-1,b2-1+(b1-1)*(batchcount-1));
      
      m=[tpc(b1,b2) tnc(b1,b2) tsc(b1,b2) ...
         spc(b1,b2) spposc(b1,b2) spnegc(b1,b2) ...
         snc(b1,b2) sac(b1,b2) ...
         ssc(b1,b2) ssposc(b1,b2) ssnegc(b1,b2)];
      bar(m);
      title(sprintf('%d v %d',batch(b1),batch(b2)))
      
      axis([0 length(m)+1 -0.5 1]);
      
      xticks(1:length(m),labels);
   end
end

keyboard


fprintf('TUNING SUMMARY:\n');
fprintf('%5s %5s %5s  %5s %5s  %5s %5s\n',...
        'BATCH','OR','ORBW','SF','SFBW','LAT','ARAT');
for batchidx=1:batchcount,
   fprintf('%5d ',batch(batchidx));
   fprintf('%5.1f %5.1f  %5.2f %5.2f  %5.1f %5.2f\n',...
           fitdata(batchidx).orfit(1:2),fitdata(batchidx).sffit(1:2),...
           fitdata(batchidx).timefit(1),...
           (fitdata(batchidx).timefit(3)-fitdata(batchidx).timefit(6)) ./ ...
           fitdata(batchidx).timefit(3));
end











