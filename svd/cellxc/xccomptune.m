% function r=xccomptune(kernfiles/cellid,batch)
%
function r=xccomptune(cellid,batch)

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

fprintf('xccomptune.m:  %d files:\n',batchcount);

% load kernel files
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

bootcount=size(strf,2);
tcount=size(strf(1).h,2);

kernfmt=strf(1).parms.kernfmt;
if strcmp(kernfmt,'fft') | strcmp(kernfmt,'lin2'),
   phasecount=4;
elseif length(kernfmt)>4 & ...
      (kernfmt(end-1)=='+' | kernfmt(end-1)=='-'),
   phasecount=str2num(kernfmt(end));
else
   phasecount=1;
end

%phasecount=1;
spacecount=size(strf(1).h,1)/phasecount;

obincount=15;
sfbincount=8;

ttime=zeros(tcount,bootcount,batchcount);
ttimepos=zeros(tcount,bootcount,batchcount);
ttimeneg=zeros(tcount,bootcount,batchcount);
tctimepos=zeros(tcount,bootcount,batchcount);
tctimeneg=zeros(tcount,bootcount,batchcount);
ttimesep=zeros(tcount,bootcount,batchcount);
tspacemaxpos=zeros(spacecount,bootcount,batchcount);
tspacemaxneg=zeros(spacecount,bootcount,batchcount);
tspaceall=zeros(spacecount,bootcount,batchcount);
tspacesep=zeros(spacecount,bootcount,batchcount);
torall=zeros(obincount,bootcount,batchcount);
tsfall=zeros(sfbincount,bootcount,batchcount);
teigvals=zeros(tcount,bootcount,batchcount);

% get temporal response functions
for bootidx=1:bootcount,
   for batchidx=1:batchcount,
      h=strf(batchidx,bootidx).h;
      
      % normalize according to review bias correction passband
      if isfield(strf(batchidx,bootidx),'powunbiased'),
         if batchidx==1,
            commonpowunbiased=strf(batchidx,bootidx).powunbiased;
            if max(commonpowunbiased)>1,
               commonpowunbiased=commonpowunbiased./max(commonpowunbiased);
            end
         else
            powunbiased=strf(batchidx,bootidx).powunbiased;
            if max(powunbiased)>1,
               powunbiased=powunbiased./max(powunbiased);
            end
            
            ratub=commonpowunbiased./powunbiased;
            ratub(find(ratub>1))=1;
            ratub(find(ratub<0))=0;
            
            h=h .* repmat(ratub,1,tcount);
         end
      end
      if phasecount>1,
         h=squeeze(sum(reshape(h,spacecount,phasecount,tcount),2));
      end
      
      % positive and negative temporal responses
      ttime(:,bootidx,batchidx)=sum(h,1)';
      ttimepos(:,bootidx,batchidx)=sum(h.*(h>0),1)';
      ttimeneg(:,bootidx,batchidx)=sum(h.*(h<0),1)';
      
      ch=cumsum(h,2);
      tctimepos(:,bootidx,batchidx)=sum(ch.*(ch>0),1)';
      tctimeneg(:,bootidx,batchidx)=sum(ch.*(ch<0),1)';
      
      % space-time separable division of kernel
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

clear r

% get average spatial/temporal kernel and associated error for each
% class by averaging over jackknifes
timepos=squeeze(mean(ttimepos,2));
etimepos=squeeze(std(ttimepos,1,2)) .* sqrt(bootcount-1);
timeneg=squeeze(mean(ttimeneg,2));
etimeneg=squeeze(std(ttimeneg,1,2)) .* sqrt(bootcount-1);
timesep=squeeze(mean(ttimesep,2));
etimesep=squeeze(std(ttimesep,1,2)) .* sqrt(bootcount-1);
spacesep=squeeze(mean(tspacesep,2));
espacesep=squeeze(std(tspacesep,1,2)) .* sqrt(bootcount-1);

r.tempresp=squeeze(mean(ttime,2));
r.etempresp=squeeze(std(ttime,1,2)) .* sqrt(bootcount-1);

tpostot=squeeze(sum(ttime.*(ttime>0)));
tnegtot=squeeze(sum(ttime.*(ttime<0)));

% pick time of first sig time lag, last sig time lag, peak
% exitation, peak suppression
t0=zeros(batchcount,1);
t1=zeros(batchcount,1);
tmin=zeros(batchcount,1);
tmax=zeros(batchcount,1);
for batchidx=1:batchcount,
   
   trstd=std(r.tempresp(:,batchidx));
   t0(batchidx)=min(find(abs(r.tempresp(:,batchidx))>trstd));
   t1(batchidx)=max(find(abs(r.tempresp(:,batchidx))>trstd));
   
   tmin(batchidx)=min(find(timeneg(:,batchidx)==min(timeneg(:,batchidx))));
   tmax(batchidx)=min(find(timepos(:,batchidx)==max(timepos(:,batchidx))));

   if t1(batchidx)<t0(batchidx)
      t1(batchidx)=t0(batchidx);
   end
   if tmin(batchidx)<t0(batchidx)
      tmin(batchidx)=t0(batchidx);
   end
   if tmax(batchidx)<t0(batchidx)
      tmax(batchidx)=t0(batchidx);
   end
end

[t0 t1 tmax tmin]

% get spatial response functions
for bootidx=1:bootcount,
   for batchidx=1:batchcount,
      h=strf(batchidx,bootidx).h;
      
      % normalize according to review bias correction passband
      if isfield(strf(batchidx,bootidx),'powunbiased'),
         if batchidx==1,
            commonpowunbiased=strf(batchidx,bootidx).powunbiased;
            if max(commonpowunbiased)>1,
               commonpowunbiased=commonpowunbiased./max(commonpowunbiased);
            end
         else
            powunbiased=strf(batchidx,bootidx).powunbiased;
            if max(powunbiased)>1,
               powunbiased=powunbiased./max(powunbiased);
            end
            
            ratub=commonpowunbiased./powunbiased;
            ratub(find(ratub>1))=1;
            ratub(find(ratub<0))=0;
            
            h=h .* repmat(ratub,1,tcount);
         end
      end
      
      if phasecount>1,
         h=squeeze(sum(reshape(h,spacecount,phasecount,tcount),2));
      end
      
      % get spatial profile at max,min,all time lag
      buidx=batchidx;
      tspacemaxpos(:,bootidx,batchidx)=sum(h(:,t0(buidx):tmax(buidx)),2);
      tspacemaxneg(:,bootidx,batchidx)=sum(h(:,t0(buidx):tmin(buidx)),2);
      tspaceall(:,bootidx,batchidx)=sum(h(:,t0(buidx):t1(buidx)),2);
      
      Xmax=sqrt(spacecount*2);
      [cfilt,cfiltconj,cmask,cmaskconj]=gencfilt(Xmax,Xmax);
      tsf=zeros(Xmax*Xmax,1);
      tsf(cfilt,:)=tspaceall(:,bootidx,batchidx);
      tsf(cfiltconj,:)=tspaceall(:,bootidx,batchidx);
      tsf=reshape(tsf,Xmax,Xmax);
      tgr=sf2gr(tsf,obincount,sfbincount);
      
      [tor,s,tsf]=svd(tgr);
      if mean(tgr,1)*tsf(:,1) > 0;
         torall(:,bootidx,batchidx)=tor(:,1);
         tsfall(:,bootidx,batchidx)=tsf(:,1);
      else
         torall(:,bootidx,batchidx)=-tor(:,1);
         tsfall(:,bootidx,batchidx)=-tsf(:,1);
      end
   end
end

spacemaxpos=squeeze(mean(tspacemaxpos,2));
espacemaxpos=squeeze(std(tspacemaxpos,1,2)) .* sqrt(bootcount-1);
spacemaxneg=squeeze(mean(tspacemaxneg,2));
espacemaxneg=squeeze(std(tspacemaxneg,1,2)) .* sqrt(bootcount-1);
spaceall=squeeze(mean(tspaceall,2));
espaceall=squeeze(std(tspaceall,1,2)) .* sqrt(bootcount-1);
orall=squeeze(mean(torall,2));
eorall=squeeze(std(torall,1,2)) .* sqrt(bootcount-1);
sfall=squeeze(mean(tsfall,2));
esfall=squeeze(std(tsfall,1,2)) .* sqrt(bootcount-1);


% compute s-t separability and std err
tsepfrac=teigvals(1,:,:).^2./sum(teigvals(1:3,:,:).^2,1);
r.sepidx=squeeze(mean(tsepfrac,2))';
r.esepidx=squeeze(std(tsepfrac,1,2))' .* sqrt(bootcount-1);


tadaptidx=zeros(size(tpostot));
tsuppidx=zeros(tcount,bootcount,batchcount);

for bootidx=1:bootcount
   for b1=1:batchcount,
      
      spsum=ttimepos(:,bootidx,b1);
      snsum=ttimeneg(:,bootidx,b1);
      
      spsum(find(spsum==0 & snsum==0))=1;
      tsuppidx(:,bootidx,b1)=-snsum./(spsum-snsum);
      
      spvec=tspacemaxpos(:,bootidx,b1);
      spvec=spvec(find(spvec~=0));
      spvec=spvec-mean(spvec);
      
      spsum=sum(spvec.^2.*(spvec>0));
      snsum=sum(spvec.^2.*(spvec<0));
      if spsum==0 & snsum==0,
         tsuppidx(1,bootidx,b1)=0;
      else
         tsuppidx(1,bootidx,b1)=sqrt(snsum./(spsum+snsum));
      end
      
      %spvec=tspacemaxneg(:,bootidx,b1);
      %spvec=tspacesep(:,bootidx,b1);
      spvec=tspaceall(:,bootidx,b1);
      spvec=spvec(find(spvec~=0));
      spvec=spvec-mean(spvec);
      
      spsum=sum(spvec.^2.*(spvec>0));
      snsum=sum(spvec.^2.*(spvec<0));
      
      if spsum==0 & snsum==0,
         tsuppidx(2,bootidx,b1)=0;
      else
         tsuppidx(2,bootidx,b1)=sqrt(snsum./(spsum+snsum));
      end
      
      %tvec=ttimesep(:,bootidx,b1);
      tvec=ttime(:,bootidx,b1);
      tvec=tvec-mean(tvec([1 end]));
      
      tpsum=sum(tvec.*(tvec>0));
      tnsum=-sum(tvec.*(tvec<0));
      
      if tpsum==0 & tnsum==0,
         tadaptidx(bootidx,b1)=0;
      else
         tadaptidx(bootidx,b1)=(tnsum./(tpsum+tnsum));
      end
   end
end

r.suppidx=squeeze(mean(tsuppidx,2));
r.esuppidx=squeeze(std(tsuppidx,1,2)) .* sqrt(bootcount-1);
r.adaptidx=squeeze(mean(tadaptidx));
r.eadaptidx=squeeze(std(tadaptidx,1)) .* sqrt(bootcount-1);


if 0,
   %keyboard
   
   for b1=1:batchcount,
      for b2=1:batchcount,
         
         r.tpc(b1,b2)=dist2(timepos(:,b1),timepos(:,b2),...
                            etimepos(:,b1),etimepos(:,b2));
         r.etpc(b1,b2)=0;
         r.tnc(b1,b2)=dist2(timeneg(:,b1),timeneg(:,b2),...
                            etimeneg(:,b1),etimeneg(:,b2));
         r.etnc(b1,b2)=0;
         r.tsc(b1,b2)=dist2(timesep(:,b1),timesep(:,b2),...
                            etimesep(:,b1),etimesep(:,b2));
         r.etsc(b1,b2)=0;
         
         r.spc(b1,b2)=dist2(spacemaxpos(:,b1),spacemaxpos(:,b2),...
                            espacemaxpos(:,b1),espacemaxpos(:,b2),0);
         r.espc(b1,b2)=0;
         r.spposc(b1,b2)=dist2(spacemaxpos(:,b1),spacemaxpos(:,b2),...
                            espacemaxpos(:,b1),espacemaxpos(:,b2),1);
         r.espposc(b1,b2)=0;
         r.spnegc(b1,b2)=dist2(spacemaxpos(:,b1),spacemaxpos(:,b2),...
                            espacemaxpos(:,b1),espacemaxpos(:,b2),-1);
         r.espnegc(b1,b2)=0;
         
         r.snc(b1,b2)=dist2(spacemaxneg(:,b1),spacemaxneg(:,b2),...
                            espacemaxneg(:,b1),espacemaxneg(:,b2));
         r.esnc(b1,b2)=0;
         r.sac(b1,b2)=dist2(spaceall(:,b1),spaceall(:,b2),...
                            espaceall(:,b1),espaceall(:,b2),0);
         r.esac(b1,b2)=0;
         r.saposc(b1,b2)=dist2(spaceall(:,b1),spaceall(:,b2),...
                            espaceall(:,b1),espaceall(:,b2),1);
         r.esaposc(b1,b2)=0;
         r.sanegc(b1,b2)=dist2(spaceall(:,b1),spaceall(:,b2),...
                            espaceall(:,b1),espaceall(:,b2),-1);
         r.esanegc(b1,b2)=0;
         
         r.ssc(b1,b2)=dist2(spacesep(:,b1),spacesep(:,b2),...
                            espacesep(:,b1),espacesep(:,b2),0);
         r.essc(b1,b2)=0;
         r.ssposc(b1,b2)=dist2(spacesep(:,b1),spacesep(:,b2),...
                               espacesep(:,b1),espacesep(:,b2),1);
         r.essposc(b1,b2)=0;
         r.ssnegc(b1,b2)=dist2(spacesep(:,b1),spacesep(:,b2),...
                               espacesep(:,b1),espacesep(:,b2),-1);
         r.essnegc(b1,b2)=0;
      end
   end
   
   r.batch=batch;
   r.cellid=cellid;
   r.labels={'tpc','tnc','tsc','spc','spc+','spc-',...
             'snc','sac','ssc','ssc+','ssc-'};
   
   return
   
end



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
orc=zeros(bootcount,batchcount,batchcount);
sfc=zeros(bootcount,batchcount,batchcount);

for b1=1:batchcount,
   for b2=b1:batchcount,
      for bootidx=1:bootcount,
         
         % corr for temp stuff
         if sum(abs(ttimepos(:,bootidx,b1))) & ...
               sum(abs(ttimepos(:,bootidx,b2))),
            tpc(bootidx,b1,b2)=dist1(ttimepos(:,bootidx,b1),...
                                     ttimepos(:,bootidx,b2),0);
            tpc(bootidx,b2,b1)=tpc(bootidx,b1,b2);
         end
         if sum(abs(ttimeneg(:,bootidx,b1))) & ...
               sum(abs(ttimeneg(:,bootidx,b2))),
            tnc(bootidx,b1,b2)=dist1(ttimeneg(:,bootidx,b1),...
                                     ttimeneg(:,bootidx,b2),0);
            tnc(bootidx,b2,b1)=tnc(bootidx,b1,b2);
         end
         if sum(abs(ttimesep(:,bootidx,b1))) & ...
               sum(abs(ttimesep(:,bootidx,b2))),
            tsc(bootidx,b1,b2)=dist1(ttimesep(:,bootidx,b1),...
                                     ttimesep(:,bootidx,b2),0);
            tsc(bootidx,b2,b1)=tsc(bootidx,b1,b2);
         end
         
         % corr for space stuff
         ts1=tspacemaxpos(:,bootidx,b1);
         ts2=tspacemaxpos(:,bootidx,b2);
         if sum(abs(ts1)) & sum(abs(ts2)),
            spc(bootidx,b1,b2)=dist1(ts1,ts2,0);
            spc(bootidx,b2,b1)=spc(bootidx,b1,b2);
         end
         if sum(abs(ts1.*(ts1>0))) & sum(abs(ts2.*(ts2>0))),
            spposc(bootidx,b1,b2)=dist1(ts1,ts2,1);
            spposc(bootidx,b2,b1)=spposc(bootidx,b1,b2);
         end
         if sum(abs(ts1.*(ts1<0))) & sum(abs(ts2.*(ts2<0))),
            spnegc(bootidx,b1,b2)=dist1(ts1,ts2,-1);
            spnegc(bootidx,b2,b1)=spnegc(bootidx,b1,b2);
         end
         
         if sum(abs(tspacemaxneg(:,bootidx,b1))) & ...
               sum(abs(tspacemaxneg(:,bootidx,b2))),
            snc(bootidx,b1,b2)=dist1(tspacemaxneg(:,bootidx,b1),...
                                     tspacemaxneg(:,bootidx,b2));
            snc(bootidx,b2,b1)=snc(bootidx,b1,b2);
         end
         
         ts1=torall(:,bootidx,b1);
         ts2=torall(:,bootidx,b2);
         if sum(abs(ts1)) & sum(abs(ts2)),
            orc(bootidx,b1,b2)=dist1(ts1,ts2,0);
            orc(bootidx,b2,b1)=orc(bootidx,b1,b2);
         end
         ts1=tsfall(:,bootidx,b1);
         ts2=tsfall(:,bootidx,b2);
         if sum(abs(ts1)) & sum(abs(ts2)),
            sfc(bootidx,b1,b2)=dist1(ts1,ts2,0);
            sfc(bootidx,b2,b1)=sfc(bootidx,b1,b2);
         end
         
         ts1=tspaceall(:,bootidx,b1);
         ts2=tspaceall(:,bootidx,b2);
         if sum(abs(ts1)) & sum(abs(ts2)),
            sac(bootidx,b1,b2)=dist1(ts1,ts2,0);
            sac(bootidx,b2,b1)=sac(bootidx,b1,b2);
         end
         if sum(abs(ts1.*(ts1>0))) & sum(abs(ts2.*(ts2>0))),
            saposc(bootidx,b1,b2)=dist1(ts1,ts2,1);
            saposc(bootidx,b2,b1)=saposc(bootidx,b1,b2);
         end
         if sum(abs(ts1.*(ts1<0))) & sum(abs(ts2.*(ts2<0))),
            sanegc(bootidx,b1,b2)=dist1(ts1,ts2,-1);
            sanegc(bootidx,b2,b1)=sanegc(bootidx,b1,b2);
         end
         %if sum(abs(tspaceall(:,bootidx,b1))) & ...
         %      sum(abs(tspaceall(:,bootidx,b2))),
         %   sac(bootidx,b1,b2)=dist1(tspaceall(:,bootidx,b1),...
         %                               tspaceall(:,bootidx,b2));
         %   sac(bootidx,b2,b1)=sac(bootidx,b1,b2);
         %end
         
         ts1=tspacesep(:,bootidx,b1);
         ts2=tspacesep(:,bootidx,b2);
         if sum(abs(ts1)) & sum(abs(ts2)),
            ssc(bootidx,b1,b2)=dist1(ts1,ts2,0);
            ssc(bootidx,b2,b1)=ssc(bootidx,b1,b2);
         end
         if sum(abs(ts1.*(ts1>0))) & sum(abs(ts2.*(ts2>0))),
            ssposc(bootidx,b1,b2)=dist1(ts1,ts2,1);
            ssposc(bootidx,b2,b1)=ssposc(bootidx,b1,b2);
         end
         if sum(abs(ts1.*(ts1<0))) & sum(abs(ts2.*(ts2<0))),
            ssnegc(bootidx,b1,b2)=dist1(ts1,ts2,-1);
            ssnegc(bootidx,b2,b1)=ssnegc(bootidx,b1,b2);
         end
         
      end
   end
end

r.etpc=squeeze(std(tpc,1,1)) .* sqrt(bootcount-1);
r.tpc=squeeze(mean(tpc,1));
r.etnc=squeeze(std(tnc,1,1)) .* sqrt(bootcount-1);
r.tnc=squeeze(mean(tnc,1));
r.etsc=squeeze(std(tsc,1,1)) .* sqrt(bootcount-1);
r.tsc=squeeze(mean(tsc,1));

r.espc=squeeze(std(spc,1,1)) .* sqrt(bootcount-1);
r.spc=squeeze(mean(spc,1));
r.espposc=squeeze(std(spposc,1,1)) .* sqrt(bootcount-1);
r.spposc=squeeze(mean(spposc,1));
r.espnegc=squeeze(std(spnegc,1,1)) .* sqrt(bootcount-1);
r.spnegc=squeeze(mean(spnegc,1));
r.esnc=squeeze(std(snc,1,1)) .* sqrt(bootcount-1);
r.snc=squeeze(mean(snc,1));
r.esac=squeeze(std(sac,1,1)) .* sqrt(bootcount-1);
r.sac=squeeze(mean(sac,1));
r.esaposc=squeeze(std(saposc,1,1)) .* sqrt(bootcount-1);
r.saposc=squeeze(mean(saposc,1));
r.esanegc=squeeze(std(sanegc,1,1)) .* sqrt(bootcount-1);
r.sanegc=squeeze(mean(sanegc,1));
r.essc=squeeze(std(ssc,1,1)) .* sqrt(bootcount-1);
r.ssc=squeeze(mean(ssc,1));
r.essposc=squeeze(std(ssposc,1,1)) .* sqrt(bootcount-1);
r.ssposc=squeeze(mean(ssposc,1));
r.essnegc=squeeze(std(ssnegc,1,1)) .* sqrt(bootcount-1);
r.ssnegc=squeeze(mean(ssnegc,1));
r.eorc=squeeze(std(orc,1,1)) .* sqrt(bootcount-1);
r.orc=squeeze(mean(orc,1));
r.esfc=squeeze(std(sfc,1,1)) .* sqrt(bootcount-1);
r.sfc=squeeze(mean(sfc,1));
r.batch=batch;
r.cellid=cellid;
r.labels={'tpc','tnc','tsc','spc','spc+','spc-',...
          'snc','sac','ssc','ssc+','ssc-'};

return

%plot the comp statistics
figure

for b1=1:batchcount,
   for b2=b1+1:batchcount,
      subplot(batchcount-1,batchcount-1,b2-1+(b1-1)*(batchcount-1));
      
      m=[r.tpc(b1,b2) r.tnc(b1,b2) r.tsc(b1,b2) ...
         r.spc(b1,b2) r.spposc(b1,b2) r.spnegc(b1,b2) ...
         r.snc(b1,b2) r.sac(b1,b2) ...
         r.ssc(b1,b2) r.ssposc(b1,b2) r.ssnegc(b1,b2)];
      bar(m);
      title(sprintf('%d v %d',batch(b1),batch(b2)))
      
      axis([0 length(m)+1 -0.5 1]);
      
      xticks(1:length(m),r.labels);
   end
end


return


function r=dist1(a,b,sign);

if ~exist('sign','var') | ~sign,
   % do nothing
elseif sign==-1,
   a=a.*(a<0);
   b=b.*(b<0);
else
   a=a.*(a>0);
   b=b.*(b>0);
end

useidx=find(a-mean(a)~=0 | b-mean(b)~=0);

if 0 & length(useidx)>0,
   r=sqrt(mean((a(useidx)-b(useidx)).^2))./ ...
     sqrt(mean(a(useidx).^2)+mean(b(useidx).^2));
   
elseif length(useidx)>0,
   r=xcorr(a(useidx),b(useidx),0,'coeff');
   %r=xcov(a(useidx),b(useidx),0,'coeff');
else
   r=0;
end


function r=dist2(a,b,aerr,berr,sign);

if ~exist('sign','var') | ~sign,
   % do nothing
elseif sign==-1,
   a=a.*(a<0);
   b=b.*(a<0);
else
   a=a.*(a>0);
   b=b.*(a>0);
end

aerr(find(aerr==0))=1;
berr(find(berr==0))=1;

useidx=find((a~=0 | b~=0) & ~isnan(a) & ~isnan(b));

if length(useidx)>0,
   a=a(useidx);
   b=b(useidx);
   aerr=aerr(useidx);
   berr=berr(useidx);
   
   r=sqrt(mean((a-b).^2./(aerr.^2+berr.^2)))./ ...
     sqrt(mean(a.^2./aerr.^2)+mean(b.^2./berr.^2));
else
   r=1;
end

return


a=a./aerr;
b=b./berr;

%r=mean((a-b).^2)./sqrt(mean(a.^2).*mean(b.^2));
%return

SEC=0.0;

a=(a).*(abs(a)>SEC);
b=(b).*(abs(b)>SEC);

useidx=find((a~=0 | b~=0) & ~isnan(a) & ~isnan(b));

if length(useidx)>0 & sum(abs(a))>0 & sum(abs(b))>0,
   r=xcorr(a(useidx),b(useidx),0,'coeff');
else
   r=0;
end
