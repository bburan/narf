% 5/17/03:
% mean one trial=0.240 (32.95%)
% mean maxxc=0.389 (86.83%)
% mean asympxc=0.418 (100.00%)
% expidx=1 nlinidx=1 cnf=1: 0.163-->0.278
% expidx=1 nlinidx=1 cnf=2: 0.029-->0.049
% expidx=1 nlinidx=1 cnf=3: 0.098-->0.167
% expidx=2 nlinidx=1 cnf=1: 0.060-->0.102
% expidx=2 nlinidx=1 cnf=2: 0.062-->0.105
% expidx=2 nlinidx=1 cnf=3: 0.033-->0.056
% expidx=3 nlinidx=1 cnf=1: 0.067-->0.113
% expidx=3 nlinidx=1 cnf=2: 0.043-->0.074
% expidx=3 nlinidx=1 cnf=3: 0.182-->0.310
% expidx=1 nlinidx=2 cnf=1: 0.234-->0.408
% expidx=1 nlinidx=2 cnf=2: 0.042-->0.074
% expidx=1 nlinidx=2 cnf=3: 0.098-->0.171
% expidx=2 nlinidx=2 cnf=1: 0.079-->0.138
% expidx=2 nlinidx=2 cnf=2: 0.313-->0.546
% expidx=2 nlinidx=2 cnf=3: 0.095-->0.165
% expidx=3 nlinidx=2 cnf=1: 0.115-->0.201
% expidx=3 nlinidx=2 cnf=2: 0.179-->0.313
% expidx=3 nlinidx=2 cnf=3: 0.311-->0.541

dbopen;
global BATQUEUEID
BATQUEUEID=[];

RECALC=0;
if ~RECALC,
   load /auto/k5/david/tmp/v1asympval.mat
   
else
   
   batches=[51; ...   % rev
            55; ...   % gr
            59];      % nr
   
   % load pred results for each batch
   res=loadbigdata(batches);
   predxc=permute(res.predxc,[1 2 6 4 5 3]);
   prederr=permute(res.prederr,[1 2 6 4 5 3]);
   predp=permute(res.predp,[1 2 6 4 5 3]);
   celllist=res.celllist;
   totlen=res.totlen;
   
   badcellidx=find(strcmp(celllist,'R150B') | strcmp(celllist,'93G83A'));
   predxc(badcellidx,1,:,:,:)=nan;
   
   s=size(predxc);
   respcount=s(3);
   batchcount=s(4);
   nlidx=1;
   
   asympxc=zeros(s(1:4));
   maxxc=zeros(s(1:4));
   onetrialxc=zeros(s(1:4));
   asympmse=zeros(s(1:4));
   maxmse=zeros(s(1:4));
   onetrialmse=zeros(s(1:4));
   asymprr=zeros(s(1:4));
   onetrialrr=zeros(s(1:4));
   trialcount=zeros(s(1:2));
   
   for cidx=1:length(celllist),
      for expidx=1:size(predxc,2),
         if ~isnan(predxc(cidx,expidx,1,expidx,nlidx)),
            
            fprintf('%s exp=%d\n',celllist{cidx},expidx);
            
            %[aa,xc,onexc,batch]=...
            %     xctestpreds(celllist{cidx},batches(expidx,nlinidx));
            
            res=xctestpreds(celllist{cidx},batches(expidx),nlidx);
            
            asympxc(cidx,expidx,:,:)=...
                reshape(res.asympxc,1,1,respcount,batchcount);
            maxxc(cidx,expidx,:,:)=...
                reshape(res.xc,1,1,respcount,batchcount);
            onetrialxc(cidx,expidx,:,:)=...
                reshape(res.onetrialxc,1,1,respcount,batchcount);
            asympmse(cidx,expidx,:,:)=...
                reshape(res.asympmse,1,1,respcount,batchcount);
            maxmse(cidx,expidx,:,:)=...
                reshape(res.mse,1,1,respcount,batchcount);
            onetrialmse(cidx,expidx,:,:)=...
                reshape(res.onetrialmse,1,1,respcount,batchcount);
            trialcount(cidx,expidx)=res.trialcount;
            
            asymprr(cidx,expidx,:,:)=...
                reshape(res.asymprr,1,1,respcount,batchcount);
            onetrialrr(cidx,expidx,:,:)=...
                reshape(res.onetrialrr,1,1,respcount,batchcount);
         end
      end
   end
   
   % so we can skip the long slow prediction step next time around.
   % ie, if RECALC=0, this gets loaded next time.
   save /auto/k5/david/tmp/v1asympval.mat
end

%
% Ok, we've calculated onetrialxc, xc, asympxc for each cell. we
% want to use this information to get contributions of infinite
% validation data and then incorporate v1asympsum analysis to see
% effect of infinite estimation data.
%
%

disp(' check asympxc matrix!!!! ');
goodidx=find(maxxc(:,1,nlinidx,1)~=0);
onetoasymppercent=squeeze(mean(asympxc(goodidx,1,:,1).^2 - ...
                               onetrialxc(goodidx,1,:,1).^2));

tlen=totlen(goodidx,1);
frac=[1.0 0.95 0.9 0.85 0.8];
datainv=1;
PMIN=0.05;

figure(1);
clf

prset=squeeze(onetrialxc(goodidx,1,:,1)).^2;
subplot(2,2,1);
[axc1,meanxc,pprset,bfitgood]=expnoise(prset,tlen,frac,datainv,PMIN);
fprintf('1: meanxc=%.2f axc=%.2f Ngood=%d\n',...
        (meanxc),(axc1),sum(bfitgood));

% reduce to only the cells that were good fits for onetrialxc
gidx=goodidx(find(bfitgood));

prset=squeeze(maxxc(gidx,1,:,1)).^2;
subplot(2,2,2);
[axcN,meanxc,pprset,bfitgood]=expnoise(prset,tlen,frac,datainv,PMIN);
fprintf('N: meanxc=%.2f axc=%.2f Ngood=%d\n',...
        (meanxc),(axcN),sum(bfitgood));

prset=squeeze(asympxc(gidx,1,:,1)).^2;
subplot(2,2,3);
[axcInf,meanxc,pprset,bfitgood]=expnoise(prset,tlen,frac,datainv,PMIN);
fprintf('inf: meanxc=%.2f axc=%.2f Ngood=%d\n',...
        (meanxc),(axcInf),sum(bfitgood));

subplot(2,2,4);
f=1./(frac.*mean(tlen));
plot(f,squeeze(1./mean(onetrialxc(gidx,1,:,1))),'k--');
hold on
a=axis;
axis([0 a(2) 1./.7 a(4)]);
plot(0,1./axc1,'o');
plot(f,squeeze(1./mean(maxxc(gidx,1,:,1))),'k-');
plot(0,1./axcN,'o');
plot(f,squeeze(1./mean(asympxc(gidx,1,:,1))),'k:');
plot(0,1./axcInf,'o');
hold off

set(gca,'YTickLabel',round(100./str2num(get(gca,'YTickLabel')))./100);
set(gca,'XTickLabel',round(100./str2num(get(gca,'XTickLabel')))./100);

figure(2);
f=frac.*mean(tlen);
plot(f,squeeze(mean(onetrialxc(gidx,1,:,1))),'k--');
hold on
a=axis;
axis([f(end) f(1).*1.5 0 .7]);
plot(f(1).*1.5,axc1,'o');
plot(f,squeeze(mean(maxxc(gidx,1,:,1))),'k-');
plot(f(1).*1.5,axcN,'o');
plot(f,squeeze(mean(asympxc(gidx,1,:,1))),'k:');
plot(f(1).*1.5,axcInf,'o');
hold off




disp('rev pred rev:');
fprintf('mean one trial=%.3f (%.2f%%)\n',...
        mean(onetrialxc(goodidx,1,2,1)),...
        mean(onetrialxc(goodidx,1,2,1)).^2./mean(asympxc(goodidx,1,2,1)).^2.*100);
fprintf('mean maxxc=%.3f (%.2f%%)\n',...
        mean(maxxc(goodidx,1,2,1)),...
        mean(maxxc(goodidx,1,2,1)).^2./mean(asympxc(goodidx,1,2,1)).^2.*100);
fprintf('mean asympxc=%.3f (%.2f%%)\n',...
        mean(asympxc(goodidx,1,2,1)),...
        mean(asympxc(goodidx,1,2,1)).^2./mean(asympxc(goodidx,1,2,1)).^2.*100);

for nlinidx=1:size(predxc,3),
   if onetoasymppercent(nlinidx)>0,
      for expidx=1:size(predxc,2),
         for cnfidx=1:size(predxc,5),
            goodidx=find(~isnan(predxc(:,expidx,nlinidx,cnfidx,nlidx)));
            
            txc=onetrialxc(goodidx,expidx,nlinidx,cnfidx);
            
            replaceidx=find(txc==0);
            
            txc(replaceidx)=predxc(goodidx(replaceidx),expidx,...
                                   nlinidx,cnfidx,nlidx);
            
            mxc=mean(txc);
            axc=mxc./sqrt(onetoasymppercent(nlinidx));
            
            fprintf('expidx=%d nlinidx=%d cnf=%d: %.3f-->%.3f\n',...
                    expidx,nlinidx,cnfidx,mxc,axc);
         end
      end
   end
end



figure
goodidx=find(maxxc(:,1,1,1)~=0);
subplot(2,2,1);
scatter(maxxc(goodidx,1,1,1),asympxc(goodidx,1,1,1));
hold on
plot([-0.2 1],[-0.2 1],'k--');
hold off
axis equal
title(sprintf('n trial 100%% (%.3f) v asymp 100%% (%.3f)',...
              mean(maxxc(goodidx,1,1,1)),...
              mean(asympxc(goodidx,1,1,1))));

subplot(2,2,2);
scatter(onetrialxc(goodidx,1,1,1),asympxc(goodidx,1,1,1));
hold on
plot([-0.2 1],[-0.2 1],'k--');
hold off
axis equal
title(sprintf('1 trial 100%% (%.3f) v asymp 100%% (%.3f)',...
              mean(onetrialxc(goodidx,1,1,1)),...
              mean(asympxc(goodidx,1,1,1))));

subplot(2,2,3);
scatter(onetrialxc(goodidx,1,4,1),onetrialxc(goodidx,1,1,1));
hold on
plot([-0.2 1],[-0.2 1],'k--');
hold off
axis equal
title(sprintf('1 trial 80%% (%.3f) v 100%% (%.3f)',...
              mean(onetrialxc(goodidx,1,4,1)),...
              mean(onetrialxc(goodidx,1,1,1))));

subplot(2,2,4);
scatter(asympxc(goodidx,1,4,1),asympxc(goodidx,1,1,1));
hold on
plot([-0.2 1],[-0.2 1],'k--');
hold off
axis equal
title(sprintf('asymp 80%% (%.3f) v 100%% (%.3f)',...
              mean(asympxc(goodidx,1,4,1)),...
              mean(asympxc(goodidx,1,1,1))));




%
%
% variable exploratory data analysis. ripped off of v1asympsum.m
%
%

badcellidx=find(strcmp(celllist,'R150B') | strcmp(celllist,'93G83A'));
predxc(badcellidx,1,:,:,:)=nan;


MINLEN=500;
expstr={'Nat','DGrat','DNat'};
instr={'pix','pow','psf'};
cnfstr=expstr;
outstr={'lin','rec','full'};

% if 1, take all cells
%PMIN=1.0;
PMIN=0.05;

DOMEDIAN=0;

exprange=[1 1 1 1;  3 3 3 3; 2 2 2 2 ];
%cnfrange=[1 1 1 1; 3 3 3 3; 2 2 2 2];
cnfrange=[1 1 1 1 ; 1 1 1 1; 1 1 1 1 ];
nlrange= [2 2 2 2; 2 2 2 2 ; 2 2 2 2 ];
trialcountrange=[1 2 3 4; 1 2 3 4; 1 2 3 4];
frac=[1.0 0.95 0.9 0.85];

%
% calculate ands display asymptotic confirmatory preds
%
figure(1);
clf
rowcount=size(exprange,2);
colcount=size(exprange,1);

mxc=zeros(size(exprange));
axc=zeros(size(exprange));
bfitgood=zeros([size(predxc,1) size(exprange)]);

clear r
r(colcount,rowcount).prset=[];

for ii=1:length(exprange(:)),
   if trialcountrange(ii)==1,
      r(ii).prset=squeeze(onetrialmse(:,exprange(ii),:,cnfrange(ii))).^2;
   elseif trialcountrange(ii)==2,
      r(ii).prset=squeeze(asympmse(:,exprange(ii),:,cnfrange(ii))).^2;
   elseif trialcountrange(ii)==3,
      r(ii).prset=squeeze(onetrialxc(:,exprange(ii),:,cnfrange(ii)));
      r(ii).prset=1./(r(ii).prset+(r(ii).prset==0)).*(r(ii).prset~=0);
   elseif trialcountrange(ii)==4,
      r(ii).prset=squeeze(asympxc(:,exprange(ii),:,cnfrange(ii)));
      r(ii).prset=1./(r(ii).prset+(r(ii).prset==0)).*(r(ii).prset~=0);
   end
   
   goodidx=find(r(ii).prset(:,1)~=0 & ~isnan(r(ii).prset(:,1)));
   
   if length(goodidx)>0,
      
      r(ii).prset=r(ii).prset(goodidx,:);
      r(ii).celllist={celllist{goodidx}};
      lens=totlen(goodidx,exprange(ii));
      ll=lens*frac;
      
      % find asymptote for each cell.
      r(ii).pp=zeros(length(goodidx),2);
      r(ii).ppe=zeros(length(goodidx),1);
      r(ii).pprset=zeros(size(r(ii).prset));
      r(ii).ppasymp=zeros(length(goodidx),1);
      for cc=1:length(goodidx),
         f=1./(frac.*lens(cc));  % ie, x-axis is 1/T
         r(ii).pp(cc,:)=polyfit(f,r(ii).prset(cc,:),1);   % y axis is predxc
         r(ii).pprset(cc,:)=(r(ii).pp(cc,2)+r(ii).pp(cc,1).*f);
         
         % error of linear fit, to determine if this cell qualifies
         if mean(r(ii).prset(cc,:).^2)==0,
            r(ii).ppe(cc)=1;
         else
            r(ii).ppe(cc)=sqrt(mean((r(ii).prset(cc,:)-...
                                     r(ii).pprset(cc,:)).^2) ./ ...
                         mean(r(ii).prset(cc,:).^2));
         end
         
         if trialcountrange(ii)<=2,
            r(ii).ppasymp(cc)=r(ii).pp(cc,2);
         else
            r(ii).ppasymp(cc)=1./r(ii).pp(cc,2);
         end
      end
      
      bfitgood(goodidx,ii)=r(ii).ppe<PMIN;
      fitidx=find(bfitgood(goodidx,ii));
      
      % find asymptotic pred for mean/median across all cells.
      
      f=1./(frac.*mean(lens(fitidx)));
      pp=polyfit(f,mean(r(ii).prset(fitidx,:)),1);
      pm=polyfit(f,median(r(ii).prset(fitidx,:)),1);
      if trialcountrange(ii)<=2,
         ppmean=pp;
         ppmed=pm;
      else
         ppmean=1./pp;
         ppmed=1./pm;
      end
      
      if DOMEDIAN,
         axc(ii)=ppmed(2);
         if trialcountrange(ii)<=2,
            mxc(ii)=median(r(ii).prset(fitidx,end));
         else
            mxc(ii)=median(1./r(ii).prset(fitidx,end));
         end
      else
         axc(ii)=ppmean(2);
         if trialcountrange(ii)<=2,
            mxc(ii)=mean(r(ii).prset(fitidx,end));
         else
            mxc(ii)=mean(1./r(ii).prset(fitidx,end));
         end
      end
      
      subplot(rowcount,colcount,ii);
      plot([f 0],pp(2)+pp(1).*[f 0],'r--');
      hold on
      plot(f,mean(r(ii).prset(fitidx,:)),'bx');
      errorbar(f,mean(r(ii).prset(fitidx,:)),...
               std(r(ii).prset(fitidx,:))./sqrt(length(fitidx)),'bx');
      plot(0,pp(2),'rx');
      hold off
      
      title(sprintf('exp=%s/cnf=%s, (%d) asymp %.3f = %.3f + %.3f ',...
                    expstr{exprange(ii)},expstr{cnfrange(ii)},...
                    trialcountrange(ii),axc(ii),mxc(ii),...
                    (axc(ii)-mxc(ii))));
      xt=get(gca,'XTick');
      xt(xt~=0)=1./xt(xt~=0);
      %set(gca,'XTickLabel',num2str(xt'));
      
   end
end

msedata=cat(2,r(1).ppasymp,r(1).prset(:,1),r(4).ppasymp,r(4).prset(:,1));
bfg=find(r(1).ppe<PMIN);
msedata=msedata(bfg,:);
sdata=sortrows(msedata,4);

figure(2)
clf
bar(sdata);

figure(3);

goodidx=find(onetrialmse(:,1,1,1)>0 & ~isnan( onetrialmse(:,1,1,1)));
goodidx=goodidx(bfg);

plot(squeeze(mean(onetrialmse(goodidx,1,:,1))),'k:');
hold on
plot(squeeze(mean(maxmse(goodidx,1,:,1))),'r--');
plot(squeeze(mean(asympmse(goodidx,1,:,1))));
hold off

% r.^2 one trial
otrr=repmat(onetrialrr(goodidx,1,1,1).^2,[1 2]);
r2one=(otrr-msedata(:,[2 1])) ./ otrr;

% r.^2 asymp
itrr=repmat(asymprr(goodidx,1,1,1).^2,[1 2]);
r2asymp=(itrr-msedata(:,[4 3])) ./ itrr;




