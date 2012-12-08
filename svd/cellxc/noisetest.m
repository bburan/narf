

%cellid='R211A';
%batchid=31;
%latency=4;
cellid='R221A';
batchid=31;
latency=3;

if ~exist('truestim','var') | ~exist('trueresp','var'),
   params=xcsetparams(cellid,batchid);
   
   params.stimloadparms={0,0,16,0,1};
   params.stimfilterparms={0,0,1,1,0};
   
   times=params.times;
   tfileidx=times(1).fileidx;
   tstart=times(1).start;
   tstop=times(1).stop;
   
   [truestim,trueresp]=xcloadstimresp(tfileidx,tstart,tstop,params);
   
   spacecount=size(truestim,2);
   
   % white noise
   modelstim=randn(size(truestim));
   modelCss=eye(spacecount);
   modelCss(55,45)=0.5;
   modelCss(45,55)=0.5;
   modelstim=modelstim*modelCss;
   
   % just use natural stim
   modelstim=truestim;
   
   
else
   disp('assuming currently loaded truestim and trueresp are ok!');
end

cutoff=60;  %spacecount;  % MAX=spacecount


MODEL=0;
if MODEL,
   noiserat=0.5;
   hmodel=zeros(spacecount,1);
   hmodel(55)=1;
   modelresp=modelstim*hmodel;
   resp=modelresp./std(modelresp).*sqrt((1-noiserat)) + ...
        randn(size(resp)).*sqrt(noiserat);
   stim=modelstim;
   hmodel=hmodel./std(modelresp).*sqrt((1-noiserat));
   r0=resp;
else
   resp=trueresp;
   stim=truestim;
   r0=shift(resp,[-latency 1]);
end

goodidx=find(~isnan(r0));
T=length(goodidx);
r0=r0(goodidx);
s0=stim(goodidx,:);

r0=r0-mean(r0);
s0=s0-repmat(mean(s0,1),[T 1]);

sigresp=std(r0);

%
% sta is "true" spike-triggered average
%
sta=s0'*r0./T;
stv=sqrt(mean((s0.*repmat(r0,[1 spacecount])-...
               repmat(sta',[length(r0) 1])).^2,1)') ./ sqrt(T);
Css=s0'*s0./T;
Perr=0.05;

[u,s,v]=svd(Css);
ss=diag(s);

%
% h is "true" decorrelated kernel
%
esta=u'*sta;
Cssi=v * diag([1./ss(1:cutoff); zeros(spacecount-cutoff,1)]) * u';
h=Cssi * sta;
herr=sqrt(mean(((s0*Cssi').*repmat(r0,[1 spacecount])-...
               repmat(h',[length(r0) 1])).^2,1)') ./ sqrt(T);

bootcount=15;
bootsta=zeros(spacecount,bootcount);
booth=zeros(spacecount,bootcount);
bCss=zeros(spacecount,spacecount,bootcount);
totalstd=zeros(bootcount,1);
bootnoisestd=zeros(bootcount,1);
prednoisestd=zeros(bootcount,1);
predxc=zeros(bootcount,1);
fprintf('booting (%d):\n',bootcount);
for rlen=round([0.6 0.7 0.8 0.9 1.0].*length(r0)),
   bootstep=rlen./bootcount;
   fprintf('rlen=%d: ',rlen);
   for ii=1:bootcount,
      rridx=[1:round((ii-1).*bootstep) round(ii.*bootstep+1):rlen];
      predidx=round((ii-1).*bootstep+1):round(ii.*bootstep);
      rr=r0(rridx);
      
      bootsta(:,ii)=s0(rridx,:)'*rr./length(rr);
      bCss(:,:,ii)=s0(rridx,:)'*s0(rridx,:)./length(rr);
      [ui,si,vi]=svd(bCss(:,:,ii));
      
      ssi=diag(si);
      booth(:,ii)=vi(:,1:cutoff) * diag(1./ssi(1:cutoff)) * ...
          ui(:,1:cutoff)' * bootsta(:,ii);
      
      tr=s0(rridx,:)*booth(:,ii);
      d1=sum(tr.^2);
      if d1>0,
         scf=sum(rr.*tr)./d1;
      else
         scf=1;
      end
      booth(:,ii)=booth(:,ii) .* scf;
      
      rpred=r0(predidx);
      pp=s0(predidx,:)*booth(:,ii);
      
      totalstd(ii)=std(rr);
      bootnoisestd(ii)=std(tr.*scf-rr);
      prednoisestd(ii)=std(rpred-pp);
      predxc(ii)=xcov(rpred,pp,0,'coeff');
      
      fprintf('.');
   end
   fprintf(' expmse=%.3f mse=%.3f xc=%.3f / %.3f\n',...
           mean(bootnoisestd),mean(prednoisestd),mean(predxc),mean(totalstd));
end
ebootsta=u'*bootsta;


% model of noise:
% E(randsta)=E(erandsta)=0
% E(randsta^2)  = E( (s0'*rr).^2 )    = E(s0.^2) * E(rr.^2)
% E(erandsta^2) = E( (u'*s0'*rr).^2 ) = E[(u'*s0).^2] * E(rr.^2)
%
% ie, noise kernel does NOT appear to have correlation bias!!!

% estimate from bootstrapped data:
bstaerr=std(bootsta,0,2) .* sqrt(bootcount-1);
bestaerr=std(ebootsta,0,2) .* sqrt(bootcount-1);
bherr=std(booth,0,2) .* sqrt(bootcount-1);


%bherr=u(:,1:cutoff)*(bestaerr(1:cutoff)./ss(1:cutoff));


% estimate from the noise model:
predrstaerr=sqrt(var(stim)'.*sigresp.^2./T);
predrestaerr=sqrt(var(stim*u)'*sigresp.^2./T);

% estimate from the programmed noise model:
% note: diag(Css)==var(s0)'
% use sigresp or mean(bootnoisestd) for noise power?

predbstaerr=sqrt(diag(Css).*sigresp.^2./T);
predbestaerr=sqrt(ss.*sigresp.^2./T);
predbherr=sqrt(var(s0*Cssi')'.*sigresp.^2./T);

% this is confusing??? or is it?  it seems like it works out
% there's something very confusing about the T squared term.
%predbherr=sqrt(diag(diag(Cssi * Css)) * ones(spacecount,1).*...
%               mean(bootnoisestd).^2./T .^2);


figure(1);
clf
subplot(3,1,1);
semilogy(bstaerr,'g');
hold on
semilogy(predrstaerr);
semilogy(stv);
semilogy(predbstaerr,'r');
hold off
title(sprintf('Cell %s (%d) - STA: noise energy vs. spatial channel',...
              cellid,batchid));

subplot(3,1,2);
semilogy(bestaerr,'g');
hold on
%semilogy(restaerr);
%semilogy(predrestaerr);
semilogy(predbestaerr,'r');
hold off
title(sprintf('Cell %s (%d) - STA: noise energy vs. eigenvector',cellid,batchid));
xlabel('eigenvector');
ylabel('stderr');

subplot(3,1,3);
plot((bherr),'g');
hold on
plot((herr));
plot((predbherr),'r');
hold off
title(sprintf('Cell %s (%d) - h: noise energy vs. spatial channel',cellid,batchid));
xlabel('spatial channel');
ylabel('stderr');

legend('boot measured','sig h','boot pred');

set(gcf,'PaperOrientation','landscape','PaperPosition',[0.25 0.25 10.5 8]);

figure(2);
clf
subplot(3,1,1);
plot([1 spacecount],[0 0],'k--');
hold on
plot([1 spacecount],[3 3],'k--');
plot([1 spacecount],[-3 -3],'k--');
plot((esta./bestaerr));
hold off
title('stderrs from zero for each eigenvector');

subplot(3,1,2);
plot([1 spacecount],[0 0],'k--');
hold on
if MODEL,
   plot(h-hmodel,'g-');
   plot(bherr,'r:');
   plot(-bherr,'r:');
else
   plot(h);
   plot(h+bherr.*3,'r:');
   plot(h-bherr.*3,'r:');
end
hold off
title('kernel with 99% error bars');

subplot(3,1,3);
plot(diag(Cssi*Css));
title('fraction of each spatial channel preserved after regularization');


if MODEL,
   sqrt(ss'*(u'*(h-hmodel)).^2)
end
sqrt(ss(1:cutoff)'*(bestaerr(1:cutoff)./ss(1:cutoff)).^2)
sqrt(ss'*(u'*(bherr)).^2)./sqrt(128)


% show eigs
%showstim(u(:,1:64),'pfft',[16 16],8,8);


set(gcf,'PaperOrientation','landscape','PaperPosition',[0.25 0.25 10.5 8]);


return


%
% 
%
randcount=100;
randsta=zeros([spacecount randcount]);
for ii=1:randcount,
   if 1
      rr=shuffle(r0);   % randomly shuffle resp -- Poisson distribution
   else
      rr=randn(size(rr)).*sigresp;   % gaussian noise with same variance
   end
   
   randsta(:,ii)=s0'*rr./T;
end

erandsta=u'*randsta;

% estimate from the (simulated noise) data:
rstaerr=std(randsta,0,2);
restaerr=std(erandsta,0,2);

erandlo=zeros(size(esta));
erandhi=zeros(size(esta));
for xx=1:spacecount,
   tt=sort(erandsta(xx,:));
   
   erandlo(xx)=tt(round(randcount.*Perr));
   erandhi(xx)=tt(round(randcount.*(1-Perr)));
end
