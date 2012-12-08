% xcsepcore.m: wrapper for xccore to do space-time separable kernels
%
% assumes everything has been set up as if to call xccore
% directly.  returns the best fit strf
%
% params.nlidx determines which output nl gets saved
% strf(1) : full strf (strict regularization to underfit)
%
%

disp('xcsepcore.m: STARTING');

%
% 1. get rough estimate of the full kernel
%
parambak=params;

% force it to always re-calculate rough strf
if 1 | ~exist('bootidx','var') | bootidx==1 ,
   
   fprintf('1. Finding rough estimate of the full kernel:\n');
   
   % strict regularization to underfit
   %params.sfsstep=max([params.sfsstep/2 2.5]);
   params.sfscount=10;
   params.sffiltsigma=3;
   if params.resampcount>12,
      params.resampcount=12;
   end
   
   firstseg=1;
   
   % use xccorr and xcfit always since it's faster and we're not too
   % concerned about getting the perfect fit here
   
   %if params.repexclude,
   %   % do repeated exclusions
   %   xcitercore;
   %else
      xccore;
   %end
   
   %if isfield(params,'altfit'),
   %   eval(params.altfit);
   %else
      xcfit;
   %end
   
   nlidx=params.nlidxsave;
   % pull out spatial resposne function by summing over kernel...
   hfull=strf(nlidx).h;
   hsum=sum(hfull,2);
   htime=sum(hfull,1);
   hsump=sum(hfull.*(hfull>0),2);
   
   % spatial kernel is first pc of the full space-time kernel
   [u,s,v]=svd(hfull);
   
   % make sure the sign is correct... for display, really. 
   % (well, also for some of the adaptation index measurements too)      
   if (htime*v(:,1) > 0 & hsum'*u(:,1) > 0) | hsump'*u(:,1) > 0,
      hspace=u(:,1);
      tempresp0=v(:,1)'.*s(1);
   else
      hspace=-u(:,1);
      tempresp0=-v(:,1)'.*s(1);
   end
   
   eigrat=s(1).^2./sum(diag(s).^2);
   fprintf('eigrat=%.2f\n',eigrat);
   
   % save marginals to original kernel
   strf(nlidx).tempresp=tempresp0;
   strf(nlidx).tempresp0=tempresp0;
   strf(nlidx).hspace=hspace;
   sampidx=round(linspace(1,params.sfscount,6));
   strf(nlidx).mH=mH(:,:,sampidx);
   strf(nlidx).eH=eH(:,:,sampidx);
   strf(nlidx).sSA2=mSA2;
   strf(nlidx).hspacebiased=hspace;
   strf(nlidx).pctvar=pctvar(sampidx);
   
   savestrf=strf(nlidx);
   bexpxc(1)=expxc(nlidx);
   bxc{1}=xc(:,:,nlidx);
else
   fprintf('1. Using rough STRF from first jackknife.\n');
   
   bexpxc(1)=vexpxc(1,1,segidx);
   %savestrf=vstrf(1,1,segidx);
   savestrf=vstrf(1,1,end);
   mSA2=savestrf.sSA2;
end

clear sH H ttSA tsSA2 tSR tmR tn tSA1
clear tSA sSA1 sSA2 sSA0 SR
clear strf

%
% save un-projected stimulus
%
stimbak=stim;
movlen=size(stim,1);
mSA2save=mSA2;

if params.fitfrac>0,
   fdata0=fdata;
end

%
% 2. fit temporal response with projection onto optimal spatial response
%
fprintf('2. Finding optimal temporal kernel, given rough spatial kernel:\n');

params=parambak;
params.sfscount=8;

if params.meansub,
   mstimbak=savestrf(1).mS;
else
   mstimbak=zeros(size(stimbak,2),1);
end

hspace=savestrf(1).hspace;
if var(hspace)>0,
   
   % project stim onto spatial response function
   stim=(stimbak-repmat(mstimbak',size(stimbak,1),1))*hspace;
   
   if params.fitfrac>0,
      fdata.stim=(fdata0.stim-...
                  repmat(mstimbak',size(fdata.stim,1),1))*hspace;
   else
      clear fdata
   end
   
   % do rc
   firstseg=1;
   spacecount=1;
   
   if params.repexclude,
      xcitercore;
   else
      xccore;
   end
else
   mH=zeros(1,diff(params.maxlag)+1,params.sfscount);
   eH=zeros(1,diff(params.maxlag)+1,params.sfscount);
   mSall=0;
end

% set ends to zeros - hacky kludge, don't do it
%mHbak=mH;
%mH(:,:,1:end-3)=0;

% choose regularization and shrinkage for temporal kernel

% if fit set is external, need to pre-process fit data and feed in
% to xcfit

if isfield(params,'altfit'),
   eval(params.altfit);
else
   xcfit;
end

%mH=mHbak;

% record result appropriately in strf structure
nlidx=params.nlidxsave;
strf(nlidx).tempresp=strf(nlidx).h;
strf(nlidx).tempresp0=tempresp0;
strf(nlidx).hspace=hspace;
strf(nlidx).h=hspace*strf(nlidx).h;
strf(nlidx).mS=mstimbak;
strf(nlidx).powunbiased=savestrf(1).powunbiased;
sampidx=round(linspace(1,params.sfscount,6));
strf(nlidx).mH=mH(:,:,sampidx);
strf(nlidx).eH=eH(:,:,sampidx);
strf(nlidx).sSA2=mSA2save;

if 1,
   % re-measure threshold -- this may help if the last couple
   % bins of the temporal response were non-zero -- a la the
   % zeroing out that happened after xcfit
   
   % scale linear kernel to match stim space -- this isn't
   % really necessary (the threshold is) but why not be tidy?
   linpred=kernpredict(strf(nlidx).h,...
                       (stimbak'-repmat(mstimbak,1,size(stimbak,1))),1,0);
   tgoodidx=find(~isnan(resp));
   r0=resp(tgoodidx)-mean(resp(tgoodidx));
   r1=linpred(tgoodidx)-mean(linpred(tgoodidx));
   d1=sum(r1.^2);
   if d1>0,
      scf=sum(r0.*r1)./d1;
   else
      scf=1;
   end
   
   % adjust kernel by scaling factor
   strf(nlidx).h=strf(nlidx).h .* scf;
   strf(nlidx).hspace=strf(nlidx).hspace .* scf;
   
   % find optimal threshold, given the new scaling
   linpred=kernpredict(strf(nlidx).h,...
                       (stimbak'-repmat(mstimbak,1,size(stimbak,1))),1,0);
   strf(nlidx).nltype='post-sum thresh';
   strf(nlidx).nlparms=findthresh(linpred(tgoodidx),resp(tgoodidx),0);
end
strf(nlidx).hspacebiased=hspace;
strf(nlidx).pctvar=pctvar(sampidx);

savestrf(2)=strf(nlidx);
bexpxc(2)=expxc(nlidx);
bxc{2}=xc(:,:,nlidx);



%
% 3. now use optimal temporal response to get optimal spatial kernel!
%
fprintf('3. Finding optimal spatial kernel, given optimal timecourse:\n');

htime=savestrf(2).tempresp;
spacecount=size(stimbak,2);
th=ones(spacecount,1)*htime;

if params.meansub,
   mstimbak=mean(stimbak,1)';
else
   mstimbak=zeros(size(stimbak,2),1);
end

stim=kernpredict(th,(stimbak-repmat(mstimbak',size(stimbak,1),1))',...
                 spacecount,0,1);
if params.fitfrac>0,
   fdata.stim=kernpredict(th,(fdata0.stim-repmat(mstimbak',...
                    size(fdata0.stim,1),1))',spacecount,0,1);
else
   clear fdata
end

% do ic
firstseg=1;
params=parambak;

params.maxlag=[0 0];
if length(params.sffiltsigma)==1,
   fprintf(' sfs: %d/%.1f sffilt: %d\n',...
           params.sfscount,params.sfsstep,params.sffiltsigma);
else
   fprintf(' sfs: %d/%.1f sffilt: %.1f-%.1f (%d)\n',...
           params.sfscount,params.sfsstep,params.sffiltsigma(1),...
           params.sffiltsigma(end),length(params.sffiltsigma));
end

if params.repexclude,
   % do repeated exclusions
   xcitercore;
else
   % simply do the XC. standard
   xccore;
end

if 0,
   for resampidx=1:resampcount,
      fprintf('scg refinement: resampidx=%d\n',resampidx);
      s0=stim([1:rstartidx(resampidx)-1 rendidx(resampidx)+1:end],:);
      r0=resp([1:rstartidx(resampidx)-1 rendidx(resampidx)+1:end]);
      gidx=find(~isnan(r0));
      r0=r0(gidx);
      s0=s0(gidx,:)-repmat(mS(:,resampidx)',[length(gidx) 1]);
      
      sfsidx=round(params.sfscount./2);
      h0=H(:,:,sfsidx,1,resampidx)';
      h1=h0;
      %h1=h0+randn(size(h0)).*std(h0)./2;
      h2=xcscg(s0,r0,h1);
      H(:,1,params.sfscount,1,resampidx)=h2';
   end
end

% choose regularization and shrinkage for spatial kernel
clear tstrf strf
if isfield(params,'altfit'),
   eval(params.altfit);
else
   xcfit;
end

% skip re-fitting temporal kernel ... seems to mess things up
REDOTIME=0;
if REDOTIME,
   % now go back and optimally fit time to new spatial kernel!
   nlidx=params.nlidxsave;
   hspace=strf(nlidx).h;
   sampidx=round(linspace(1,params.sfscount,6));
   ttmH=mH(:,:,sampidx);
   tteH=eH(:,:,sampidx);
   
   % project stim onto spatial response function
   %stim=stimbak*hspace;
   stim=(stimbak-repmat(mstimbak',size(stimbak,1),1))*hspace;
   
   % do rc
   firstseg=1;
   spacecount=1;
   params.maxlag=parambak.maxlag;
   params.sfscount=8;
   xccore;
   
   % choose regularization and shrinkage for temporal kernel
   clear strf tstrf fdata
   if isfield(params,'altfit'),
      eval(params.altfit);
   else
      xcfit;
   end
   
   nlidx=params.nlidxsave;
   strf(nlidx).tempresp=strf(nlidx).h;
   strf(nlidx).tempresp0=htime;
   strf(nlidx).hspace=hspace;
   strf(nlidx).h=hspace*strf(nlidx).h;
   strf(nlidx).mS=mstimbak;
   strf(nlidx).mH=ttmH;
   strf(nlidx).eH=tteH;
   strf(nlidx).sSA2=mSA2;
else
   % skip re-fitting temporal kernel ... seems to mess things up
   nlidx=params.nlidxsave;
   strf(nlidx).tempresp=htime;
   strf(nlidx).tempresp0=htime;
   strf(nlidx).hspace=strf(nlidx).h;
   strf(nlidx).h=strf(nlidx).h*htime;
   strf(nlidx).mS=mstimbak;
   sampidx=round(linspace(1,params.sfscount,6));
   strf(nlidx).mH=mH(:,:,sampidx);
   strf(nlidx).eH=eH(:,:,sampidx);
   strf(nlidx).sSA2=mSA2;
end

if 1,
   % scale linear kernel to match stim space is this necessary?
   linpred=kernpredict(strf(nlidx).h,...
                       (stimbak'-repmat(mstimbak,1,size(stimbak,1))),1,0);
   linpred=thresh(strf(nlidx).nlparms,linpred);
   tgoodidx=find(~isnan(resp));
   r0=resp(tgoodidx)-mean(resp(tgoodidx));
   r1=linpred(tgoodidx)-mean(linpred(tgoodidx));
   d1=sum(r1.^2);
   if d1>0,
      scf=sum(r0.*r1)./d1;
   else
      scf=1;
   end
   
   % adjust kernel by scaling factor
   strf(nlidx).h=strf(nlidx).h .* scf;
   strf(nlidx).hspace=strf(nlidx).hspace .* scf;
   
   % find optimal threshold, given the new scaling
   linpred=kernpredict(strf(nlidx).h,...
                       (stimbak'-repmat(mstimbak,1,size(stimbak,1))),1,0);
   strf(nlidx).nltype='post-sum thresh';
   strf(nlidx).nlparms=findthresh(linpred(tgoodidx),resp(tgoodidx),0);
end
strf(nlidx).hspacebiased=hspace;
strf(nlidx).pctvar=pctvar(sampidx);

savestrf(3)=strf(nlidx);
bexpxc(3)=expxc(nlidx);
bxc{3}=xc(:,:,nlidx);
expxc=bexpxc;

% nlidx=3 is the interesting one
xc=bxc;

% convert to format for validation
strfcount=length(savestrf);

strf=savestrf(:);
stim=stimbak;
params=parambak;

% want third kernel
clear stimbak parambak seplinpred mSA2 mSA2save fdata

