%function kvatunesum(batch,forcereload[=0],recalc[=0])
%
% good generic code source for batch-wide summaries
%
% created SVD 10/2004 - ripped off of rescalesum.m
%
function kvatunesum(batch,forcereload,RECALC)

if ~exist('batch','var'),
   batch=99
end
if ~exist('forcereload','var'),
   forcereload=0;
end
if ~exist('RECALC','var'),
   RECALC=0;
end

% this is to enable the two-peak fitter
addpath ~/code/kendrick

fprintf('%s.m: batch=%d\n',mfilename,batch);
resfile=sprintf('/auto/k5/david/tmp/kvaparms/kvatunesum.batch%d.mat',batch);

if ~forcereload & exist(resfile,'file');
   RELOAD=0;
   
   % don't need to re-run all the high-level analyses
   fprintf('loading saved results: %s\n',resfile);
   load(resfile);
else
   RELOAD=1;
   
   % find all sRunData entries for this batch
   rundata=dbgetrundata(batch);
   
   clear res
   for ii=1:length(rundata);
      
      % load and process cell-specific results
      tres=kvatune(rundata(ii).cellid,batch);
      
      % stupid code to deal with stupid matlab
      if tres.bad,
         res(ii).bad=1;
      else
         res(ii)=tres;  
      end
   end
   
   clear RECALC
   save(resfile);
   RECALC=1;
end

if RECALC==1,
   disp('(re)calculating tuning properties.');
   
   cellids={res.cellid};

   % remove all cells that don't have results files
   ggidx=find(~cat(1,res.bad));
   res=res(ggidx);
   rundata=rundata(ggidx);
   
   cellcount=length(res);
   
   cc=cat(3,res.cc);
   cc0=squeeze(cc(1,1,:));
   ccrand=cat(1,res.randxc);
   cc=[squeeze(cc(2,3:end,:))' ones(cellcount,1)];
   
   if 0,
      % old attempt at measuring cart/ncart selectivity .. rather
      % ugly distribution, eh?
      figure(1);
      clf
      
      p1=[1 -.5 -.5]./norm([1 -.5 -.5]);
      p2=[0 1 -1]./norm([0 -1 1]);
      p3=[1 1 1]./norm([1 1 1]);
      P=[p1' p2' p3'];
      
      d2=dmax*P;
      
      scatter(d2(:,1),d2(:,2));
      A=5.*eye(3)*P;
      hold on
      for ii=1:size(A,2);
         plot([0 A(ii,1)],[0 A(ii,2)],'r-');
      end
      hold off
      
      keyboard
   end
   
   batchdata=dbget('sBatch',batch);
   batchdata.stimloadparms=strsep(batchdata.stimloadparms,',');
   batchdata.stimfilterparms=strsep(batchdata.stimfilterparms,',');
   Xmax=batchdata.stimloadparms{3};
   
   % relative amplitude of two excitatory gabors
   %a=squeeze(beta2(:,5,1:2));
   %a=a./repmat(sum(a,2),[1 2]);
   a=cc(:,1:2);
   a=a./repmat(a(:,2),[1 2]);
   
   % which gabor is bigger
   biga=(a(:,2)>a(:,1))+1;
   
   sf2=zeros(cellcount,2);
   owid2=zeros(cellcount,2);
   catidx=zeros(cellcount,1);
   
   beta1=cat(3,res.beta1);
   
   figure(1);
   if size(beta1,1)>0,
      beta1=permute(beta1(:,1,:),[3 1 2]);
      beta1=reshape(beta1,cellcount,5,size(beta1,2)/5);
      
      beta2=cat(3,res.beta2);
      beta2=permute(beta2(:,1,:),[3 1 2]);
      beta2=reshape(beta2,cellcount,5,size(beta2,2)/5);
      
      orpeak=mod(beta1(:,1,1).*180/pi,180);
      sfpeak=beta1(:,2,1);
      sfbw=beta1(:,4,1);
      sfnpeak=beta1(:,2,2);
      owid=beta1(:,3,1).*180/pi;
      odiff=mod((beta2(:,1,1)-beta2(:,1,2))*180/pi,180);
      odiff=90-abs(odiff-90);
      sfdiff=abs(beta2(:,2,1)-beta2(:,2,2));
      
      for ii=1:cellcount,
         sf2(ii,1)=beta2(ii,2,biga(ii));
         sf2(ii,2)=beta2(ii,2,3-biga(ii));
         owid2(ii,1)=beta2(ii,3,biga(ii))*180/pi;
         owid2(ii,2)=beta2(ii,3,3-biga(ii))*180/pi;
         catidx(ii)=(a(ii,biga(ii))>0.85)+1;
      end
   end
   
   % tuning crap
   bootcount=size(res(1).H,2);
   if strcmp(batchdata.kernfmt,'space'),
      obincount=sqrt(size(res(ii).H(:,bootidx,1),1));
      sfbincount=obincount;
   else
      obincount=16;
      sfbincount=8;
   end
   obins=linspace(0,180,obincount+1);
   obins=obins(1:end-1);
   sfbins=linspace(1,round((Xmax+1)/2),sfbincount+1);
   sfbins=sfbins(1:end-1);
   ortuning=zeros(obincount,cellcount,bootcount);
   sftuning=zeros(sfbincount,cellcount,bootcount);
   cm=zeros(cellcount,bootcount);
   cstd=zeros(cellcount,bootcount);
   mcstd=zeros(cellcount,1);
   mcm=zeros(cellcount,1);
   portuning=zeros(obincount,cellcount,bootcount);
   psftuning=zeros(sfbincount,cellcount,bootcount);
   
   orpeak=zeros(cellcount,1);
   porpeak=zeros(cellcount,1);
   pcm=zeros(cellcount,bootcount);
   pcstd=zeros(cellcount,bootcount);
   mpcstd=zeros(cellcount,1);
   seprat=zeros(cellcount,bootcount);
   pseprat=zeros(cellcount,bootcount);
   peakdist0=zeros(cellcount,bootcount);
   twopeakscore0=zeros(cellcount,bootcount);
   twopeakscore=zeros(cellcount,1);
   sfpeak0=zeros(cellcount,1);
   sfbw0=zeros(cellcount,1);
   sfpeak1=zeros(cellcount,1);
   sfbw1=zeros(cellcount,1);
   sfbw1good=zeros(cellcount,1);
   seprat2=zeros(cellcount,bootcount);
   pseprat2=zeros(cellcount,bootcount);
   
   figure(1);
   drawnow;
   
   for ii=1:cellcount
      cellids{ii}
      
      for bootidx=1:size(res(ii).H,2),
         seprat2(ii,bootidx)=pfftsep(res(ii).H(:,bootidx,1));
         if (sum((res(ii).H(:,bootidx,1)>0)))>0,
            pseprat2(ii,bootidx)=pfftsep(res(ii).H(:,bootidx,1).*...
                                         (res(ii).H(:,bootidx,1)>0));
         end
         %pause(0.1);
      end
      
      % [nanmedian(seprat2(ii,:)) nanmedian(pseprat2(ii,:))]
      
      for bootidx=1:size(res(ii).H,2),
         if strcmp(batchdata.kernfmt,'space'),
            tsfgr=reshape(res(ii).H(:,bootidx,1),obincount,sfbincount);
         else
            tsfIR=pfft2sf(res(ii).H(:,bootidx,1));
            tsfgr=sf2gr(tsfIR,obincount,sfbincount);
         end
         
         psfgr=tsfgr.*(tsfgr>0);
         
         if bootidx==1,
            if strcmp(batchdata.kernfmt,'space'),
               tsf0=reshape(mean(res(ii).H(:,:,1),2),obincount,sfbincount);
            else
               tsf0=pfft2sf(mean(res(ii).H(:,:,1),2));
               tsf0=sf2gr(tsf0,obincount,sfbincount);
            end
            mor=mean(tsf0,2);
            msf=mean(tsf0,1)';
            mpsf=mean(psfgr,1)';
         end
         
         [u,s,v]=svd(tsfgr);
         if u(:,1)'*mor >0,
            por=u(:,1); %.*s(1);
            psf=v(:,1); %.*s(1);
         else
            por=-u(:,1); %.*s(1);
            psf=-v(:,1); %.*s(1);
         end
         mor=mor+por;
         if sum(diag(s))>0,
            seprat(ii,bootidx)=s(1)./sum(diag(s));
         end
         [u,s,v]=svd(psfgr);
         if v(:,1)'*mpsf >0,
            ppor=u(:,1); %.*s(1);
            ppsf=v(:,1); %.*s(1);
         else
            ppor=-u(:,1); %.*s(1);
            ppsf=-v(:,1); %.*s(1);
         end
         if sum(diag(s))>0,
            pseprat(ii,bootidx)=s(1)./sum(diag(s));
         end
         
         % use sum, not separable decomposition to get or and sf
         % tuning curves
         if 0,
            por=mean(tsfgr,2);
            psf=mean(tsfgr',2);
         end
         
         ortuning(:,ii,bootidx)=por(:,1);
         sftuning(:,ii,bootidx)=psf(:,1);
         portuning(:,ii,bootidx)=ppor(:,1);
         psftuning(:,ii,bootidx)=ppsf(:,1);
         
         [cm(ii,bootidx),cstd(ii,bootidx)]=circstats(por(:,1));
         [pcm(ii,bootidx),pcstd(ii,bootidx)]=circstats(ppor(:,1));
         
         % adjust to deal with extra pi half of circle and then
         % convert to degrees
         cstd(ii,bootidx)=cstd(ii,bootidx)./2 .* 180/pi; 
         pcstd(ii,bootidx)=pcstd(ii,bootidx)./2 .* 180/pi;
         
         % now convert from std to width at half-height (assuming
         % that or tuning is a gaussian)
         cstd(ii,bootidx)=cstd(ii,bootidx).* sqrt(2.*log(2)).*2; 
         pcstd(ii,bootidx)=pcstd(ii,bootidx).* sqrt(2.*log(2)).*2;
         
         % can't be greater than 180 deg
         if cstd(ii,bootidx)>180,
            cstd(ii,bootidx)=180;
         end
         if pcstd(ii,bootidx)>180,
            pcstd(ii,bootidx)=180;
         end
      end
      
      if 1,
         % find cstd of ortuning to the (1./0.4)-th power to
         % compensate for exponent on movpower to get correct orbw
         % to correspond to width at half height of response as
         % the SAME stim is rotated from the peak.
         mort=mean(ortuning(:,ii,:),3);
         if sum(abs(mort))>0,
            [aa,mcstd(ii)]=circstats(abs(mort).^(1./0.4).*sign(mort));
            orpeak(ii)=mod(aa./2.*180/pi,180);
            mcstd(ii)=mcstd(ii)./2 .* 180/pi;
            mcstd(ii)=mcstd(ii).* sqrt(2.*log(2)).*2;
         else
            disp('flat or tuning curve');
            mcstd(ii)=180;
         end
         
         mort=mean(portuning(:,ii,:),3);
         if sum(abs(mort))>0,
            [aa,mpcstd(ii)]=circstats(abs(mort).^(1./0.4).*sign(mort));
            porpeak(ii)=mod(aa./2.*180/pi,180);
            mpcstd(ii)=mpcstd(ii)./2 .* 180/pi;
            mpcstd(ii)=mpcstd(ii).* sqrt(2.*log(2)).*2;
         else
            disp('flat por tuning curve');
            mpcstd(ii)=180;
         end
         
         % can't be greater than 180 deg
         if mcstd(ii)>180,
            mcstd(ii)=180;
         end
         if mpcstd(ii)>180,
            mpcstd(ii)=180;
         end
      
      else
         % alternatively fit gaussian and take params from fit
         if sum(abs(mean(ortuning(:,ii,:),3)))>0,
            obeta=fitgauss1d(obins,mean(ortuning(:,ii,:),3),[],1);
            orpeak(ii)=obeta(1);
            mcstd(ii)=(obeta(2)+obeta(2)) .* sqrt(2.*log(2^.4));
          else
            disp('flat or tuning curve');
            mcstd(ii)=180;
         end
         
         if sum(mean(portuning(:,ii,:),3))>0,
            obeta=fitgauss1d(obins,mean(portuning(:,ii,:),3),[],1);
            orpeak(ii)=obeta(1);
            mpcstd(ii)=(obeta(2)+obeta(2)) .* sqrt(2.*log(2^.4));
          else
            disp('flat por tuning curve');
            mpcstd(ii)=180;
         end
         % can't be greater than 180 deg
         if mcstd(ii)>180,
            mcstd(ii)=180;
         end
         if mpcstd(ii)>180,
            mpcstd(ii)=180;
         end
      end
      
      %sfpower=sqrt(mean(sftuning(:,ii,:).^2,3));
      sfpower=(mean(sftuning(:,ii,:),3));
      
      % sf measurements V1 : fit gaussian and pull out stats
      beta=fitgauss1d(log2(sfbins),sfpower);
      %beta=fitgauss1d(log(sfbins),sfpower);
      
      % peak sf tuning is mean of sf tuning curve
      sfpeak0(ii)=2.^beta(1);
      
      % sf bw is width at 1/2 height divided by peak sf
      %sfbw0(ii)=(2.^(beta(2) .* sqrt(2.*log(2))).*2)./sfpeak0(ii);
      %sfbw0(ii)=(2.^(beta(1)+beta(2))-2.^(beta(1)-beta(2)))./2.^beta(1);
      
      % sf bw defined by de valois is (log2 H2 - log2 L2)
      % where H2 and L2 are sfs at half height.
      %sfbw0(ii)=(log(2.^(beta(1)+beta(2)))-log(2.^(beta(1)-beta(2))))./log(2);
      %sfbw0(ii)=(beta(2)+beta(2))./log(2);
      sfbw0(ii)=(beta(2)+beta(2)) .* sqrt(2.*log(2));
      
      % sf measurements V2 : use center-of-mass
      sfnorm=sfpower;
      
      sfnorm=sfnorm - min(sfnorm .*(sfnorm<0));
      sfnorm=sfnorm./(sum(sfnorm) + (sum(sfnorm)==0));
      sfpeak1(ii)=2.^(log2(sfbins)*sfnorm / sum(sfnorm));      
      
      sfsum=sfnorm.*log2(sfbins(:));
      sfsum(sfsum==0)=min(sfsum(sfsum>0))./1000;
      sfsum=sfsum./(sum(sfsum) + (sum(sfsum)==0));
      sfsum=cumsum(sfsum);
      
      sfpoints=2.^interp1(sfsum,log2(sfbins),[.125 .5 .875]);
      
      % use real center of mass calc (above)
      %sfpeak1(ii)=sfpoints(2);
      if isnan(sfpoints(1)) & isnan(sfpoints(3)),
         sfbw1(ii)=nan;
      elseif isnan(sfpoints(1)),
         sfbw1(ii)=2.*sfpoints(3)./sfpoints(2);
      elseif isnan(sfpoints(3)),
         sfbw1(ii)=-2.*sfpoints(1)./sfpoints(2);
      else
         sfbw1(ii)=(sfpoints(3)-sfpoints(1))./sfpoints(2);
         sfbw1good(ii)=1;
      end
      
      if sfnorm(end)>max(sfnorm)/2,
         sfbw1good(ii)=-2;      % high-pass
      elseif sfnorm(1)>max(sfnorm)/2,
         sfbw1good(ii)=-1;      % low-pass
      end
      
      clf
      subplot(2,1,1);
      plot(obins,mean(ortuning(:,ii,:),3),'r--x');
      % not fitting gauss for orientation tuning
      %hold on
      %plot(linspace(0,180,100),...
      %     gauss1(obeta,linspace(0,180,100),1),'b-');
      %plot(rbins*180/pi,gauss1(obeta,rbins,1));
      %hold off
      title(sprintf('%s orpeak=%.0f orw=%.1f porw=%.1f',res(ii).cellid,...
                    orpeak(ii),mcstd(ii),mpcstd(ii)));
      
      subplot(2,1,2);
      %plot(log2(sfbins),gauss1(beta,log2(sfbins)),'b-');
      plot(linspace(log2(sfbins(1)),log2(sfbins(end))),...
           gauss1(beta,linspace(log2(sfbins(1)),log2(sfbins(end)))),'b-');
      hold on
      plot(log2(sfbins),sfpower,'r--x');
      hold off
      title(sprintf('sfpk=%.1f sfbw=%.1f',sfpeak0(ii),sfbw0(ii)));
      
      %res(ii).cellid
      [sfpeak0(ii) sfbw0(ii) sfpeak1(ii) sfbw1(ii) sfbw1good(ii)]
      
      if 0,
         por=mean(ortuning(:,ii,:),3);
         eor=std(por).*ones(size(por));
         if sum(por)==0,
            f=zeros(size(por));
         else
            f = findpeaks(por,eor,20,1,0.5,0,1);
         end
         if sum(f)>2,
            f = findpeaks(por,eor,20,1,0.4,0,1);
         end
         if sum(f)>2,
            f = findpeaks(por,eor,20,1,0.25,0,1);
         end
         if sum(f)>2,
            f = findpeaks(por,eor,20,1,0.1,0,1);
         end
         if sum(f)>2,
            f = findpeaks(por,eor,20,1,0.05,0,1);
         end
         
         peakidx=find(f);
         if length(peakidx)==2,
            min1=min(por(peakidx(1):peakidx(2)));
            min2=min(por([1:peakidx(1) peakidx(2):end]));
            
            minmax=min(por(peakidx));
            maxmin=max([min1 min2]);
            twopeakscore(ii)=(minmax-maxmin)./...
                (max([por(:);0])-min([por(:);0]));
         else
            twopeakscore(ii)=0;
         end
         %keyboard
         
         if sum(f)>2,
            disp('more than two peaks');
            twopeakscore0(ii,:)=nan;
         elseif sum(f)==2,
            peakidx=find(f);
            peakdist0(ii,:)=diff(obins(peakidx));
            if peakdist0(ii,1)>90,
               peakdist0(ii,bootidx)=peakdist0(ii,bootidx)-90;
            end
            
            min1=min(por(peakidx(1):peakidx(2)));
            min1idx=find(por==min1);
            min2=min(por([1:peakidx(1) peakidx(2):end]));
            min2idx=find(por==min2);
            if por(min1idx)>por(min2idx),
               mmm=min1idx;
               min1idx=min2idx;
               min2idx=mmm;
            end
            if por(peakidx(1))<por(peakidx(2));
               peakidx=flipud(peakidx);
            end
            
            for bootidx=1:size(res(ii).H,2),
               tpor=ortuning(:,ii,bootidx);

               %twopeakscore0(ii,bootidx)=(tpor(peakidx(2))-tpor(min2idx))./...
               %    (max([tpor(peakidx(1)) 0])-min([tpor(min1idx) 0]));
               twopeakscore0(ii,bootidx)=(tpor(peakidx(2))-tpor(min2idx))./...
                   (max([tpor(peakidx(1)) 0])-min([tpor(min1idx) 0]));
               if twopeakscore0(ii,bootidx)>1,
                  %twopeakscore0(ii,bootidx)=1;
               end
               if twopeakscore0(ii,bootidx)<0,
                  %twopeakscore0(ii,bootidx)=0;
               end
               
            end
         else
            peakdist0(ii,:)=0;
            twopeakscore0(ii,:)=0;
         end
      else
         
         por=mean(ortuning(:,ii,:),3);
         eor=std(por).*ones(size(por));
         if sum(por)==0,
            f=zeros(size(por));
         else
            f = findpeaks(por,eor,20,1,0.5,0,1);
         end
         if sum(f)>2,
            f = findpeaks(por,eor,20,1,0.4,0,1);
         end
         if sum(f)>2,
            f = findpeaks(por,eor,20,1,0.25,0,1);
         end
         if sum(f)>2,
            f = findpeaks(por,eor,20,1,0.1,0,1);
         end
         if sum(f)>2,
            f = findpeaks(por,eor,20,1,0.05,0,1);
         end
         
         peakidx=find(f);
         if length(peakidx)==2,
            min1=min(por(peakidx(1):peakidx(2)));
            min2=min(por([1:peakidx(1) peakidx(2):end]));
            
            minmax=min(por(peakidx));
            maxmin=max([min1 min2]);
            twopeakscore(ii)=(minmax-maxmin)./...
                (max([por(:);0])-min([por(:);0]));
         else
            twopeakscore(ii)=0;
         end
         
         for bootidx=1:size(res(ii).H,2),
            por=ortuning(:,ii,bootidx);
            eor=std(por).*ones(size(por));
            if sum(por)==0,
               f=zeros(size(por));
            else
               f = findpeaks(por,eor,20,1,0.5,0,1);
            end
            if sum(f)>2,
               f = findpeaks(por,eor,20,1,0.4,0,1);
            end
            if sum(f)>2,
               f = findpeaks(por,eor,20,1,0.25,0,1);
            end
            if sum(f)>2,
               f = findpeaks(por,eor,20,1,0.1,0,1);
            end
            if sum(f)>2,
               f = findpeaks(por,eor,20,1,0.05,0,1);
            end
            if sum(f)>2,
               disp('more than two peaks');
               twopeakscore0(ii,bootidx)=nan;
            elseif sum(f)==2,
               peakidx=find(f);
               peakdist0(ii,bootidx)=diff(obins(peakidx));
               if peakdist0(ii,bootidx)>90,
                  peakdist0(ii,bootidx)=peakdist0(ii,bootidx)-90;
               end
               
               min1=min(por(peakidx(1):peakidx(2)));
               min2=min(por([1:peakidx(1) peakidx(2):end]));
               
               minmax=min(por(peakidx));
               maxmin=max([min1 min2]);
               
               %twopeakscore0(ii,bootidx)=(minmax-maxmin)./max([por(:);0]);
               %twopeakscore0(ii,bootidx)=(minmax-maxmin)./mean(abs(por));
               %twopeakscore0(ii,bootidx)=(minmax-maxmin)./(max(por)-min(por));
               %twopeakscore0(ii,bootidx)=(minmax-maxmin)./...
               %    (max([por(:)])-min([por(:)]));
               twopeakscore0(ii,bootidx)=(minmax-maxmin)./...
                   (max([por(:);0])-min([por(:);0]));
               
            else
               peakdist0(ii,bootidx)=0;
               twopeakscore0(ii,bootidx)=0;
            end
         end
         
      end
      fprintf('2m: %.2f 2s: %.2f\n',mean(twopeakscore(ii,:)),...
              std(twopeakscore0(ii,:)).*sqrt(size(twopeakscore0,2)-1));
      drawnow
      
   end
   
   if bootcount>1,
      %twopeakscore=nanmedian(twopeakscore0')';
      peakdist=nanmedian(peakdist0')';
   else
      twopeakscore=twopeakscore0;
      peakdist=peakdist0;
   end
end

%mcstd=mean(cstd,2);
%mpcstd=mean(pcstd,2);
mseprat=median(seprat,2);
mpseprat=median(pseprat,2);

%mseprat=sqrt(mseprat);
%mpseprat=sqrt(mseprat);

NATFRACTHRESH=1./(384);

sparseness=cat(1,res.sparseness);
sparseplus=cat(1,res.sparseplus);
entropy=cat(1,res.entropy);

gmax=[cat(1,res.cartmax) cat(1,res.polmax) ...
      cat(1,res.polmax) cat(1,res.hypmax)];

if size(gmax,2)>4,
   %gmax=gmax(:,1:2:7); % use Hz values
   gmax=gmax(:,2:2:8); % use norm response values
end

% don't separate polar and concentric
%gmax(:,2)=0;

% do separate polar and concentric
polparms=cat(1,res.polparms);
gmax(find(polparms(:,3)<=2),2)=0;
gmax(find(polparms(:,3)>2),3)=0;

bestcat=zeros(cellcount,1);
for ii=1:cellcount,
   bestcat(ii)=min(find(gmax(ii,:)==max(gmax(ii,:))));
end
natfracbetter=cat(1,res.natfracbetter);
natovercart=cat(1,res.natovercart);

natmax=cat(1,res.natmax);
natmaxerr=cat(1,res.natmaxerr);
cartvs=cat(1,res.cartvscc);
if size(cartvs,2)>4,
   disp('temp cartvs fix!!!!')
   cartvs=[cartvs(:,1) cartvs(:,3:end)];
end
cartvserr=cat(1,res.cartvsccerr);
synthvs=cat(1,res.synthvs);
synthvserr=cat(1,res.synthvserr);

% cart vs ncart
bestcat2=ones(size(bestcat));
bestcat2(cartvs(:,2)>cartvs(:,1))=2;

bestcat2e=-bestcat2;
bestcat2e(cartvs(:,2)+cartvserr(:,2)<cartvs(:,1)-cartvserr(:,1))=1;
bestcat2e(cartvs(:,2)-cartvserr(:,2)>cartvs(:,1)+cartvserr(:,1))=2;

% cart vs nat
bestcat3=ones(size(bestcat));
bestcat3(cartvs(:,3)>cartvs(:,1))=2;

bestcat3e=-bestcat3;
bestcat3e(cartvs(:,3)+cartvserr(:,3)<cartvs(:,1)-cartvserr(:,1))=1;
bestcat3e(cartvs(:,3)-cartvserr(:,3)>cartvs(:,1)+cartvserr(:,1))=2;

% cart vs ncart vs nat
bestcat4=ones(size(bestcat));
bestcat4(cartvs(:,2)>cartvs(:,1) & cartvs(:,2)>cartvs(:,3))=2;
bestcat4(cartvs(:,3)>cartvs(:,1) & cartvs(:,3)>cartvs(:,2))=3;

% cart vs ncart vs nat
bestcat4e=-bestcat4;
xxx=sqrt(2);
bestcat4e(cartvs(:,1)-cartvserr(:,1)./xxx>cartvs(:,2)+cartvserr(:,2)./xxx & ...
          cartvs(:,1)-cartvserr(:,1)./xxx>cartvs(:,3)+cartvserr(:,3)./xxx)=1;
bestcat4e(cartvs(:,2)-cartvserr(:,2)./xxx>cartvs(:,1)+cartvserr(:,1)./xxx & ...
          cartvs(:,2)-cartvserr(:,2)./xxx>cartvs(:,3)+cartvserr(:,3)./xxx)=2;
bestcat4e(cartvs(:,3)-cartvserr(:,3)./xxx>cartvs(:,1)+cartvserr(:,1)./xxx & ...
          cartvs(:,3)-cartvserr(:,3)./xxx>cartvs(:,2)+cartvserr(:,2)./xxx)=3;

%bestcat4=bestcat;
%bestcat4(natmax(:,2)>max(gmax,[],2))=5;

%bestcat3=ones(size(bestcat));
%bestcat3(natmax(:,1)>gmax(:,1))=2;

% summarize after excluding noisy cells
fprintf('best3==1: %.3f best3==2: %.3f best5==5: %.3f best1>1: %.3f\n',...
        [sum(bestcat3==1 & cc0>ccrand) sum(bestcat3==2 & cc0>ccrand) ...
         sum(bestcat4==5 & cc0>ccrand) sum(bestcat>1 & cc0>ccrand) ]./...
        sum(cc0>ccrand));

% without excluding noisy cells
%[sum(bestcat3==1) sum(bestcat3==2) ...
% sum(bestcat4==5) sum(bestcat>1) ]./length(cc0)

catcount=zeros(max(bestcat),1);
g=zeros(max(bestcat),1);
ge=zeros(max(bestcat),1);
g2=zeros(max(bestcat),1);
g2e=zeros(max(bestcat),1);
g3=zeros(max(bestcat),1);
g3e=zeros(max(bestcat),1);
g4=zeros(max(bestcat),1);
g4e=zeros(max(bestcat),1);
Ncat=8;
catcstd=zeros(max(bestcat),Ncat);
for ii=1:max(bestcat),
   useidx=find(cc0>=ccrand & bestcat==ii);
   %useidx=find(bestcat==ii);
   if length(useidx)>0,
      g(ii)=median(mcstd(useidx));
      ge(ii)=std(mcstd(useidx)) ./ sqrt(length(useidx));
      
      g2(ii)=median(twopeakscore(useidx));
      g2e(ii)=std(twopeakscore(useidx)) ./ sqrt(length(useidx));
      
      g3(ii)=median(mseprat(useidx));
      g3e(ii)=std(mseprat(useidx)) ./ sqrt(length(useidx));

      g4(ii)=median(mseprat(useidx));
      g4e(ii)=std(mseprat(useidx)) ./ sqrt(length(useidx));
      
      g5(ii)=sqrt(mean(cc0(useidx).^2));
      g5e(ii)=std(cc0(useidx)) ./ sqrt(length(useidx));
      
      catcount(ii)=length(find(cc0>=ccrand & bestcat==ii));
      
      [catcstd(ii,:),hbins]=hist(mcstd(useidx),linspace(15,180,Ncat));
   end
end

%catcount2=hist(bestcat,1:max(bestcat));
catcount2=hist(bestcat2,1:max(bestcat2));
catcount2e=hist(bestcat2e,-2:2);
catcount2e=catcount2e([4 5; 2 1]);
catcount3=hist(bestcat3,1:max(bestcat3));
catcount3e=hist(bestcat3e,-2:2);
catcount3e=catcount3e([4 5; 2 1]);
catcount4=hist(bestcat4,1:max(bestcat4));
%catcount4=[catcount4(1) sum(catcount4(2:4)) catcount4(5)];
catcount4e=hist(bestcat4e,-3:3);
catcount4e=catcount4e([5 6 7; 3 2 1]);

catstr={'cart','rad','conc','hyp'};
catstr2={'cart','ncart','nat'};
catstr3={'cart','nat'};
catfmt={'bs','r.','ko','g^'};
catfmt2={'bs','ko','ko','ko','g^'};

cellids={res.cellid};

if RELOAD | RECALC,
   %done with stuff that gets saved to the processed data
   %file. rest of script is mostly display
   
   % save res to file specific to this batch
   clear RECALC
   save(resfile);
end


if 1,
   figure(2);
   clf
   
   acells=find(cc0>=ccrand);
   bcells=find(cc0<ccrand);
   %acells=find(bestcat3==1);
   %bcells=find(bestcat3==2);
   
   subplot(3,2,1);
   n1=hist(mpcstd(acells),10:20:170)';
   n2=hist(mpcstd(bcells),10:20:170)';
   bar(10:20:170,[n1 n2],'stacked');
   title(sprintf('batch %d mpcstd m=%.2f',batch,median(mpcstd)));
   
   subplot(3,2,2);
   n1=hist(twopeakscore(acells),0.05:0.1:0.95)';
   n2=hist(twopeakscore(bcells),0.05:0.1:0.95)';
   bar(0.05:0.1:0.95,[n1 n2],'stacked');
   title(sprintf('batch %d twopeakscore m=%.2f',batch,nanmedian(twopeakscore)));
   
   subplot(3,2,3);
   n1=hist(sfpeak0(acells),sfbins)';
   n2=hist(sfpeak0(bcells),sfbins)';
   bar(sfbins,[n1 n2],'stacked');
   title(sprintf('batch %d sfpeak0 m=%.2f',batch,median(sfpeak0)));
   
   subplot(3,2,4);
   n1=hist(sfbw0(acells),0.5:0.33:2.5)';
   n2=hist(sfbw0(bcells),0.5:0.33:2.5)';
   bar(0.5:0.33:2.5,[n1 n2],'stacked');
   title(sprintf('batch %d sfbw0  m=%.2f',batch,median(sfbw0)));
   
   if 0,
      subplot(3,2,5);
      n1=hist(sfpeak1(acells),sfbins)';
      n2=hist(sfpeak1(bcells),sfbins)';
      bar(sfbins,[n1 n2],'stacked');
      title(sprintf('batch %d sfpeak1  m=%.2f',batch,median(sfpeak1)));
      
      subplot(3,2,6);
      n1=hist(sfbw1(acells),0.5:0.33:2.5)';
      n2=hist(sfbw1(bcells),0.5:0.33:2.5)';
      bar(0.5:0.33:2.5,[n1 n2],'stacked');
      title(sprintf('batch %d sfbw1  m=%.2f',batch,median(sfbw1)));
     
   else
      subplot(3,2,5);
      sepbins=linspace(0.05,0.95,10);
      n1=hist(mseprat(acells),sepbins)';
      n2=hist(mseprat(bcells),sepbins)';
      bar(sepbins,[n1 n2],'stacked');
      title(sprintf('sf-or sep (med=%.2f)',nanmedian(mseprat(acells))));
      
      subplot(3,2,6);
      n1=hist(mpseprat(acells),sepbins)';
      n2=hist(mpseprat(bcells),sepbins)';
      bar(sepbins,[n1 n2],'stacked');
      title(sprintf('sf-or psep (med=%.2f)',nanmedian(mpseprat(acells))));
   end
   drawnow
end


if 1,
   figure(1);
   clf
   
   subplot(4,3,1);
   ccbins=linspace(0.05,0.95,10);
   n1=hist(cc0(acells),ccbins)';
   n2=hist(cc0(bcells),ccbins)';
   bar(ccbins,[n1 n2],'stacked');
   title(sprintf('predcc (median=%.2f)',sqrt(median(cc0.^2))));
   
   acells=find(bestcat4==1);
   bcells=find(bestcat4==2);
   ccells=find(bestcat4==3);
   
   gg=zeros(3,1);
   gge=zeros(3,1);
   
   subplot(4,3,4);
   
   gg(1)=median(mcstd(acells));
   gg(2)=median(mcstd(bcells));
   gg(3)=median(mcstd(ccells));
   gge(1)=std(mcstd(acells)) ./ sqrt(length(acells));
   gge(2)=std(mcstd(bcells)) ./ sqrt(length(bcells));
   gge(3)=std(mcstd(ccells)) ./ sqrt(length(ccells));
   
   errorbar(1:length(gg),gg,gge,'k+');
   hold on
   bar(1:length(gg),gg)
   hold off
   title('median circ pstd per nc class');
   xticks(1:length(gg),catstr2)
      
   subplot(4,3,7);
   
   gg(1)=median(twopeakscore(acells));
   gg(2)=median(twopeakscore(bcells));
   gg(3)=median(twopeakscore(ccells));
   gge(1)=std(twopeakscore(acells)) ./ sqrt(length(acells));
   gge(2)=std(twopeakscore(bcells)) ./ sqrt(length(bcells));
   gge(3)=std(twopeakscore(ccells)) ./ sqrt(length(ccells));
   
   errorbar(1:length(gg),gg,gge,'k+');
   hold on
   bar(1:length(gg),gg)
   hold off
   title('median twopeakscore per nc class');
   xticks(1:length(gg),catstr2)
   
   subplot(4,3,10);
   
   if 1,
      gg(1)=median(sfbw0(acells));
      gg(2)=median(sfbw0(bcells));
      gg(3)=median(sfbw0(ccells));
      gge(1)=std(sfbw0(acells)) ./ sqrt(length(acells));
      gge(2)=std(sfbw0(bcells)) ./ sqrt(length(bcells));
      gge(3)=std(sfbw0(ccells)) ./ sqrt(length(ccells));
      svar='sfbw0';
   else
      gg(1)=median(sfpeak0(acells));
      gg(2)=median(sfpeak0(bcells));
      gg(3)=median(sfpeak0(ccells));
      gge(1)=std(sfpeak0(acells)) ./ sqrt(length(acells));
      gge(2)=std(sfpeak0(bcells)) ./ sqrt(length(bcells));
      gge(3)=std(sfpeak0(ccells)) ./ sqrt(length(ccells));
      svar='sfpeak0';
   end
   
   errorbar(1:length(gg),gg,gge,'k+');
   hold on
   bar(1:length(gg),gg)
   hold off
   title(['median ',svar,' per nc class']);
   xticks(1:length(gg),catstr2)
  
   subplot(4,3,2);
   errorbar(1:max(bestcat),g,ge,'k+');
   hold on
   bar(1:max(bestcat),g)
   hold off
   title('median circ std per nc class');
   xticks(1:max(bestcat),catstr)
   
   subplot(4,3,5);
   errorbar(1:max(bestcat),g2,g2e,'k+');
   hold on
   bar(1:max(bestcat),g2)
   hold off
   title('median 2-peakiness per class');
   xticks(1:max(bestcat),catstr)
   
   subplot(4,3,8);
   errorbar(1:max(bestcat),g3,g3e,'k+');
   hold on
   bar(1:max(bestcat),g3)
   hold off
   title('median sf-or sep per class');
   xticks(1:max(bestcat),catstr)
   
   subplot(4,3,11);
   %errorbar(1:max(bestcat),g4,g4e,'k+');
   %hold on
   %bar(1:max(bestcat),g4)
   %hold off
   %title('median sf-or pos sep per class');
   %xticks(1:max(bestcat),catstr)
   plot(sort(natfracbetter),'r');
   hold on
   plot(sort(natovercart));
   hold off
   title('fraction nat preferred');
   axis([0 length(natfracbetter)+1 0 0.2]);
   xlabel('cell');
   
   subplot(4,3,3);
   
   for bestidx=[3 1 2],
      hp=plot(mcstd(cc0>=ccrand & bestcat4==bestidx,1),...
           twopeakscore(cc0>=ccrand & bestcat4==bestidx,1),catfmt{bestidx});
      ccol=get(hp,'Color');
      set(hp,'MarkerFaceColor',ccol);
      hold on;
   end
   hold off
   xlabel('circ std'); ylabel('2-peakiness');
   legend(catstr2{[3 1 2]});
   
   if 0
   subplot(4,3,6);
   for bestidx=1:max(bestcat),
      hp=plot(mseprat(cc0>=ccrand & bestcat==bestidx,1),...
           twopeakscore(cc0>=ccrand & bestcat==bestidx,1),catfmt{bestidx});
      ccol=get(hp,'Color');
      set(hp,'MarkerFaceColor',ccol);
      hold on;
   end
   hold off
   xlabel('separability'); ylabel('2-peakiness');
   legend(catstr);
   end
   
   subplot(4,3,6);
   %bar(catcount);
   bar(catcount2e','stacked');
   title('cells per class');
   xticks(1:length(catstr2),catstr2);
   
   subplot(4,3,9);
   bar(catcount3e','stacked');
   title('cells per class');
   xticks(1:length(catstr3),catstr3);
   
   subplot(4,3,12);
   bar(catcount4e','stacked');
   title('cells per class');
   xticks(1:length(catstr2),catstr2);
   
   fullpage portrait
end

% dump stimulus class selectivity info

%[xxx,sortidx]=sortrows([-cc0 bestcat]);
sortidx=(1:cellcount)';
fprintf('CELL     PXC  ORP ORW   2PK SFP SFW   CART-SEL       NCART-SEL      NAT-SEL\n');
for ii=1:length(sortidx),
   cidx=sortidx(ii);
   
   fprintf('%-7s %5.2f',cellids{cidx},cc0(cidx));
   fprintf(' %3.0f',orpeak(cidx));
   fprintf(' %5.1f %4.2f',mcstd(cidx),twopeakscore(cidx));
   fprintf(' %3.1f %3.1f',sfpeak0(cidx),sfbw0(cidx));
   fprintf('%6.1f +/-%4.1f',res(cidx).cartvscc(1),res(cidx).cartvsccerr(1));
   if (res(cidx).cartvscc(1)-res(cidx).cartvsccerr(1))- ...
         (res(cidx).cartvscc(3)+res(cidx).cartvsccerr(3)) > 0,
      fprintf('*');
   else
      fprintf(' ');
   end
   
   fprintf('%6.1f +/-%4.1f',res(cidx).cartvscc(2), ...
           res(cidx).cartvsccerr(2));
   if (res(cidx).cartvscc(2)-res(cidx).cartvsccerr(2))- ...
         (res(cidx).cartvscc(1)+res(cidx).cartvsccerr(1)) > 0,
      fprintf('*');
   else
      fprintf(' ');
   end
   
   fprintf('%6.1f +/-%4.1f',res(cidx).cartvscc(3),res(cidx).cartvsccerr(3));
   if (res(cidx).cartvscc(3)-res(cidx).cartvsccerr(3))- ...
         (res(cidx).cartvscc(1)+res(cidx).cartvsccerr(1)) > 0,
      fprintf('*');
   else
      fprintf(' ');
   end
   
   fprintf('\n');
end

%disp('paused before displaying individual srfs');
m1=mean(twopeakscore0,2);  
s1=std(twopeakscore0,0,2).*sqrt(size(twopeakscore0,2)-1);
ccgood=find(cc0>ccrand);
[sum(m1>s1) sum(m1(ccgood)>s1(ccgood)) length(ccgood)]
[mean(m1(ccgood)) median(m1(ccgood))]

keyboard

disp('returned before displaying individual srfs');
return

%[xxx,sortidx]=sortrows([-cc0 bestcat]);
sortidx=(1:cellcount)';
rowcount=10;
colcount=6;
[cfilt,cfiltconj]=gencfilt(Xmax,Xmax);
for ii=1:length(sortidx),
   figure(ceil(ii/10)+1);
   subcell=mod(ii-1,10)+1;
   if subcell==1,
      clf
      colormap(redblue);
      fullpage portrait
   end
   
   subplot(rowcount,colcount,subcell*colcount-5);
   
   %tsf=pfft2sf(mean(res(sortidx(ii)).H(:,:,1),2),'pfft');
   tsf=sf2gr(mean(res(sortidx(ii)).H(:,:,1),2),15,8,0,0,'pfft')';
   imagesc(tsf,[-1 1].*max(abs(tsf(:))));
   axis image
   axis xy
   axis off
   title(sprintf('%s %.2f',res(sortidx(ii)).cellid,cc0(sortidx(ii))));
   
   if 1,
      subplot(rowcount,colcount,subcell*6-4);
      tfit=reshape(gaussfpN(beta2(sortidx(ii),1:10)',Xmax),Xmax,Xmax);
      tfitsf=sf2gr(tfit(cfilt),obincount,sfbincount,0,0,batchdata.kernfmt)';
      imagesc(obins,sfbins,tfitsf,[-1 1].*max(abs(tfitsf(:))));
      
      %axis image
      axis xy
      axis off
      title(sprintf('%s a1=%.2f',catstr{bestcat(sortidx(ii))},...
                    a(sortidx(ii),1)));
      
      subplot(rowcount,colcount,subcell*6-3);
      plot(obins,squeeze(ortuning(:,sortidx(ii),:)));
      hold on
      plot([obins(1) obins(end)],[0 0],'k--');
      hold off
      if subcell<rowcount
         set(gca,'XTickLabel',[]);
      end
      set(gca,'YTickLabel',[]);
      axis([obins(1) obins(end) [-1 1] .* ...
            max(max(abs(squeeze(ortuning(:,sortidx(ii),:)))))+eps]);
      title(sprintf('cstd=%.2f 2pk=%.2f',...
                    mean(cstd(sortidx(ii),:)),twopeakscore(sortidx(ii))));
      
   elseif 0,
      
      [fe,fo]=gaborcurve(beta2(sortidx(ii),:,1));
      f1=[fe fo];
      f1=f1-fo(1);
      if max(abs(f1(:)))>0,
         f1=f1./max(abs(f1(:)))./2 + 0.5;
      else
         f1=f1+0.5;
      end
      [fe,fo]=gaborcurve(beta2(sortidx(ii),:,2));
      f2=[fe fo];
      f2=f2-fo(1);
      if max(abs(f2(:)))>0,
         f2=f2./max(abs(f2(:)))./2 + 0.5;
      else
         f2=f2+0.5;
      end
      
      % show gabor pair fit
      subplot(rowcount,colcount,subcell*6-4);
      imagesc(repmat(f1,[1 1 3]));
      axis image; axis off;
      title(sprintf('%s a1=%.2f',catstr{bestcat(sortidx(ii))},...
                    a(sortidx(ii),1)));
      
      subplot(rowcount,colcount,subcell*6-3);
      imagesc(repmat(f2,[1 1 3]));
      axis image; axis off;
      title(sprintf('a1=%.2f',a(sortidx(ii),2)));
   else
      
      subplot(rowcount,colcount,subcell*6-4);
      plot(obins,squeeze(ortuning(:,sortidx(ii),:)));
      hold on
      plot([obins(1) obins(end)],[0 0],'k--');
      hold off
      if subcell<rowcount
         set(gca,'XTickLabel',[]);
      end
      set(gca,'YTickLabel',[]);
      axis([obins(1) obins(end) [-1 1] .* ...
            max(max(abs(squeeze(ortuning(:,sortidx(ii),:)))))+eps]);
      title(sprintf('cstd=%.2f 2pk=%.2f',...
                    mean(cstd(sortidx(ii),:)),twopeakscore(sortidx(ii))));
      
      subplot(rowcount,colcount,subcell*6-3);
      plot(sfbins,squeeze(sftuning(:,sortidx(ii),:)));
      hold on
      plot([sfbins(1) sfbins(end)],[0 0],'k--');
      hold off
      if subcell<rowcount
         set(gca,'XTickLabel',[]);
      end
      set(gca,'YTickLabel',[]);
      axis([sfbins(1) sfbins(end) [-1 1] .* ...
            max(max(abs(squeeze(sftuning(:,sortidx(ii),:)))))+eps]);
      title(sprintf('%s a1=%.2f',catstr{bestcat(sortidx(ii))},...
                    a(sortidx(ii),1)));
   end
   
   % show optimal cart/ncart stimuli
   [t,copt]=ncart_pfft(res(sortidx(ii)).cartparms,Xmax,'cart',...
                batchdata.stimfiltercmd,batchdata.stimfilterparms,[.5 .3]);
   copt(find(copt>1))=1;
   copt(find(copt<0))=0;
   subplot(rowcount,colcount,subcell*6-2);
   imagesc(repmat(copt,[1 1 3]));
   axis image; axis off;
   title(sprintf('cart=%.2f',res(sortidx(ii)).cartmax));
   
   [t,popt]=ncart_pfft([res(sortidx(ii)).polparms],Xmax,'pol',...
              batchdata.stimfiltercmd,batchdata.stimfilterparms,[.5 .3]);
   popt(find(popt>1))=1;
   popt(find(popt<0))=0;
   subplot(rowcount,colcount,subcell*6-1);
   imagesc(repmat(popt,[1 1 3]));
   axis image; axis off;
   title(sprintf('pol=%.2f',res(sortidx(ii)).polmax));
   
   [t,hopt]=ncart_pfft(res(sortidx(ii)).hypparms,Xmax,'hyp',...
              batchdata.stimfiltercmd,batchdata.stimfilterparms,[.5 .3]);
   hopt(find(hopt>1))=1;
   hopt(find(hopt<0))=0;
   subplot(rowcount,colcount,subcell*6);
   imagesc(repmat(hopt,[1 1 3]));
   axis image; axis off;
   title(sprintf('hyp=%.2f',res(sortidx(ii)).hypmax));
end

keyboard

for figidx=2:ceil(length(sortidx)/10+1),
   figure(figidx);
   %colormap(gray);
   %print -Plj2200
   print -Pgcolor
end


addpath /auto/k1/david/code/toolbox/kmeans

H=zeros(size(res(1).H,1),cellcount);
for ii=1:cellcount
   H(:,ii)=res(ii).H(:,1);
end
[z, c] = kmeans(H,clustercount);
n=hist(c,clustercount)

D=[sfpeak sf2(:,:) sfnpeak owid2(:,:) odiff sfdiff sort(a,2)];
clustercount=5;
[z, c] = kmeans(D',clustercount);
n=hist(c,clustercount)

% print out some stats.
fprintf('cell  cat  cc0 :   a1   a2    owid  ow1  ow2  odif sdif\n');
for ii=1:length(sortidx),
   fprintf(...
      '%5s %4s %.2f: %5.2f %5.2f  %4.1f %4.1f %4.1f  %4.1f %4.1f %4.1f\n',...
      res(sortidx(ii)).cellid,catstr{bestcat(sortidx(ii))},...
      cc0(sortidx(ii)),-sort(-a(sortidx(ii),:)),...
      owid(sortidx(ii)),owid2(sortidx(ii),:),...
      odiff(sortidx(ii)),sfpeak(sortidx(ii)),sfdiff(sortidx(ii)));
end

%odiff(catidx==2)=-1;
%sfdiff(catidx==2)=-1;

figure(1);
clf
%bar(g);
%legend(catstr);
%xlabel('opt number of gabors');
%ylabel('number of cells');

for bestidx=1:max(bestcat),
   plot(sf2(cc0>=ccrand & bestcat==bestidx,1),...
        owid(cc0>=ccrand & bestcat==bestidx),catfmt{bestidx})
   hold on;
end
hold off
legend(catstr);


sortrows([bestcat owid odiff a])

disp('need to load err variables out of res structure!');
keyboard


figure(3);
clf


%n1=histc(beta(:,2,1),0:10);
%n2=histc(beta(:,2,2),0:10);
%bar(0:10,[n1 n2]);






