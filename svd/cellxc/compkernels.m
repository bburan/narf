% function r=compkernels(strf);
%
%
function r=compkernels(strf);

strfcount=size(strf,1);
bootcount=size(strf,2);
tcount=size(strf(1).h,2);

fprintf('compkernels.m:  %d strfs X %d jackknifes:\n',...
        strfcount,bootcount);

kernfmt=strf(1).parms.kernfmt;
if strcmp(kernfmt,'fft') | strcmp(kernfmt,'lin2'),
   phasecount=4;
elseif length(kernfmt)>4 & ...
      (kernfmt(end-1)=='+' | kernfmt(end-1)=='-'),
   phasecount=str2num(kernfmt(end));
else
   phasecount=1;
end

phasecount=1;
spacecount=size(strf(1).h,1)/phasecount;

obincount=15;
sfbincount=8;

ttime=zeros(tcount,bootcount,strfcount);
tspace=zeros(spacecount,bootcount,strfcount);
tspace0=zeros(spacecount,bootcount,strfcount);

ttimesep=zeros(tcount,bootcount,strfcount);
tspacesep=zeros(spacecount,bootcount,strfcount);

ttimepos=zeros(tcount,bootcount,strfcount);
ttimeneg=zeros(tcount,bootcount,strfcount);
tctimepos=zeros(tcount,bootcount,strfcount);
tctimeneg=zeros(tcount,bootcount,strfcount);
tspacemaxpos=zeros(spacecount,bootcount,strfcount);
tspacemaxneg=zeros(spacecount,bootcount,strfcount);
torall=zeros(obincount,bootcount,strfcount);
tsfall=zeros(sfbincount,bootcount,strfcount);

teigvals=zeros(tcount,bootcount,strfcount);
tadaptidx=zeros(bootcount,strfcount);
tsuppidx=zeros(bootcount,strfcount);

% get temporal response functions
for bootidx=1:bootcount,
   for batchidx=1:strfcount,
      h=strf(batchidx,bootidx).h;
      
      tspaceall=sum(h,2);
      if isfield(strf(batchidx,bootidx),'hspace'),
         hspace=strf(batchidx,bootidx).hspace;
      else
         hspace=tspaceall;
      end
      
      if isfield(strf(batchidx).parms,'sfscount');
         if batchidx==1,
            mSA2=strf(batchidx).sSA2;
            hsb=hspace;
            
         else
            
            if size(mSA2,2)==size(hspace,1),
               hspace=mSA2*hspace;
               th=normalizereg(hspace,mSA2,[],strf(1).parms.sfscount,....
                            strf(1).parms.sfsstep);
               hsb=th(:,:,strf(1).parms.sfsfit);
            else
               hsb=hspace;
            end
      
         end
         % normalize according to review bias correction passband
      else
         if isfield(strf(batchidx,bootidx),'powunbiased'),
            if batchidx==1,
               commonpowunbiased=strf(batchidx,bootidx).powunbiased;
               if max(commonpowunbiased)>1,
                  commonpowunbiased=commonpowunbiased./max(commonpowunbiased);
               end
            else
               powunbiased=strf(batchidx,bootidx).powunbiased;
               
               if min(powunbiased)<0,
                  spowunbiased=sort(powunbiased);
                  powunbiased=powunbiased - mean(spowunbiased(1:5));
                  powunbiased(find(powunbiased<0.0001))=0.0001;
               end
               if max(powunbiased)>1,
                  spowunbiased=sort(powunbiased);
                  powunbiased=powunbiased./mean(spowunbiased(end-5:end));
                  powunbiased(find(powunbiased>1))=1;
               end
               
               ratub=commonpowunbiased./powunbiased;
               ratub(find(ratub>1))=1;
               ratub(find(ratub<0))=0;
               
               h=h .* repmat(ratub,1,tcount);
               hspace=hspace .* ratub;
            end
         end
         
         % residual bias from review added
         if isfield(strf(batchidx,bootidx),'hspacebiased'),
            hsb=strf(batchidx,bootidx).hspacebiased;
         else
            hsb=hspace;
         end
      end
      
      if phasecount>1,
         h=squeeze(sum(reshape(h,spacecount,phasecount,tcount),2));
         hspace=squeeze(sum(reshape(hspace,spacecount,phasecount),2));
         hsb=squeeze(sum(reshape(hsb,spacecount,phasecount),2));
      end
      
      [u,s,v]=svd(h);
      if hspace'*u(:,1) > 0;
         ttimesep(:,bootidx,batchidx)=v(:,1);
         tspacesep(:,bootidx,batchidx)=u(:,1);
      else
         ttimesep(:,bootidx,batchidx)=-v(:,1);
         tspacesep(:,bootidx,batchidx)=-u(:,1);
      end
      teigvals(:,bootidx,batchidx)=diag(s);
      
      if isfield(strf(batchidx,bootidx),'tempresp0'),
         ttime(:,bootidx,batchidx)=strf(batchidx,bootidx).tempresp0;
      else
         ttime(:,bootidx,batchidx)=ttimesep(:,bootidx,batchidx);
      end
      
      tspace(:,bootidx,batchidx)=hsb;
      tspace0(:,bootidx,batchidx)=hspace;
      
      tvec=ttime(:,bootidx,batchidx);
      %tvec=tvec-mean(tvec([1 end]));
      
      tpsum=sum(tvec.^2.*(tvec>0));
      tnsum=sum(tvec.^2.*(tvec<0));
      %tpsum=sum(tvec.*(tvec>0));
      %tnsum=-sum(tvec.*(tvec<0));
      
      if tpsum==0 & tnsum==0,
         tadaptidx(bootidx,batchidx)=0;
      elseif 1,
         tadaptidx(bootidx,batchidx)=(tnsum./(tpsum+tnsum));
      elseif tnsum<tpsum,
         tadaptidx(bootidx,batchidx)=sqrt(2.*(tnsum./(tpsum+tnsum)));
      else
         tadaptidx(bootidx,batchidx)=sqrt(2.*(tpsum./(tpsum+tnsum)));
      end
      
      svec=tspace(:,bootidx,batchidx);
      %spsum=sum(svec.^2.*(svec>0));
      %snsum=sum(svec.^2.*(svec<0));
      spsum=sum(svec.*(svec>0));
      snsum=-sum(svec.*(svec<0));
      
      if spsum==0 & snsum==0,
         tsuppidx(bootidx,batchidx)=0;
      elseif 1,
         tsuppidx(bootidx,batchidx)=sqrt(snsum./(spsum+snsum));
      elseif snsum<spsum,
         tsuppidx(bootidx,batchidx)=sqrt(2.*(snsum./(spsum+snsum)));
      else
         tsuppidx(bootidx,batchidx)=sqrt(2.*(spsum./(spsum+snsum)));
      end
   end
end

clear r

r.tempresp=squeeze(mean(ttime,2));
r.etempresp=squeeze(std(ttime,1,2)) .* sqrt(bootcount-1);
r.hspace=squeeze(mean(tspace,2));
r.ehspace=squeeze(std(tspace,1,2)) .* sqrt(bootcount-1);

% find average adaptation index
r.adaptidx=squeeze(mean(tadaptidx,1));
r.eadaptidx=squeeze(std(tadaptidx,1,1)) .* sqrt(bootcount-1);

% find average suppression index
r.suppidx=squeeze(mean(tsuppidx,1));
r.esuppidx=squeeze(std(tsuppidx,1,1)) .* sqrt(bootcount-1);

% compute s-t separability and std err
tsepfrac=teigvals(1,:,:).^2./(sum(teigvals(:,:,:).^2,1) + ...
                              (sum(teigvals(:,:,:).^2,1)==0)) ;
r.sepidx=squeeze(mean(tsepfrac,2))';
r.esepidx=squeeze(std(tsepfrac,1,2))' .* sqrt(bootcount-1);


% compute corr between kernels
tc=zeros(bootcount,strfcount,strfcount);
sc=zeros(bootcount,strfcount,strfcount);
scp=zeros(bootcount,strfcount,strfcount);
scn=zeros(bootcount,strfcount,strfcount);

for b1=1:strfcount,
   for b2=b1:strfcount,
      for bootidx=1:bootcount,
         
         % corr for temp stuff
         tc(bootidx,b1,b2)=dist1(ttime(:,bootidx,b1),...
                                 ttime(:,bootidx,b2),0);
         tc(bootidx,b2,b1)=tc(bootidx,b1,b2);
         
         % corr for spatial stuff
         ts1=tspace(:,bootidx,b1);
         ts2=tspace(:,bootidx,b2);

         sc(bootidx,b1,b2)=dist1(ts1,ts2,0);
         sc(bootidx,b2,b1)=sc(bootidx,b1,b2);
         scp(bootidx,b1,b2)=dist1(ts1,ts2,1);
         scp(bootidx,b2,b1)=scp(bootidx,b1,b2);
         scn(bootidx,b1,b2)=dist1(ts1,ts2,-1);
         scn(bootidx,b2,b1)=scn(bootidx,b1,b2);
         
         ts1=tspace0(:,bootidx,b1);
         ts2=tspace0(:,bootidx,b2);
         sc0(bootidx,b1,b2)=dist1(ts1,ts2,0);
         sc0(bootidx,b2,b1)=sc0(bootidx,b1,b2);
         scp0(bootidx,b1,b2)=dist1(ts1,ts2,1);
         scp0(bootidx,b2,b1)=scp0(bootidx,b1,b2);
         scn0(bootidx,b1,b2)=dist1(ts1,ts2,-1);
         scn0(bootidx,b2,b1)=scn0(bootidx,b1,b2);
         
      end
   end
end

r.tc=squeeze(mean(tc,1));
r.etc=squeeze(std(tc,1,1)) .* sqrt(bootcount-1);
r.sc=squeeze(mean(sc,1));
r.esc=squeeze(std(sc,1,1)) .* sqrt(bootcount-1);
r.scp=squeeze(mean(scp,1));
r.scn=squeeze(mean(scn,1));
r.sc0=squeeze(mean(sc0,1));
r.scp0=squeeze(mean(scp0,1));
r.scn0=squeeze(mean(scn0,1));


%keyboard

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

if var(a(useidx))>0 & var(b(useidx))>0,
   r=xcorr(a(useidx),b(useidx),0,'coeff');
   %r=xcov(a(useidx),b(useidx),0,'coeff');
else
   r=0;
end



return


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



