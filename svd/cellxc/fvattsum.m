% function fvattsum(batchid,reload[=0],recalc[=0])
%
% created SVD 6/21/04 - hacked from kvaparms.m
%
function fvattsum(batchid,reload,recalc);

FVBATCHES=[95 96];
RESPATH='/auto/k5/david/tmp/kvaparms/';

resfile=sprintf('%sdmsattsum.batch%d.mat',RESPATH,batchid);
if ~exist('reload','var'),
   reload=0;
end
if ~exist('recalc','var'),
   recalc=0;
end
if ~reload & ~exist(resfile,'file'),
   reload=1;
end

if reload | recalc
   dmsattsum(batchid,reload,recalc);
end

fprintf('loading %s\n',resfile);
load(resfile);

keyboard

keepidx=[];
for ii=1:length(res),
   if ~isempty(res(ii).cellid),
      keepidx=[keepidx; ii];
   end
end
res=res(keepidx);
cellcount=length(res);

PTHRESH=0.05;

predxc=zeros([size(res(1).predxc) cellcount]);
pxc=zeros([size(res(1).pxc) cellcount]);
for ii=1:cellcount,
   predxc(:,:,ii)=res(ii).predxc(1:size(predxc,1),1:size(predxc,2));
   pxc(:,:,ii)=res(ii).pxc(1:size(pxc,1),1:size(pxc,2));
   if isempty(res(ii).randxc),
      ii
      res(ii).randxc=1;
   end
end
pxc=squeeze(pxc)';
localsig=pxc(:,5)<PTHRESH;
randxc=cat(1,res.randxc);
goodpred=squeeze(predxc(1,1,:))>randxc;

%
% looks for trends in tuning
% 


% peak spatial frequency difference in feature attention
aidx=1;
shidx=find(localsig & goodpred)';
shidx0=find(~localsig & goodpred)';
shnorm=sqrt(length(shidx));
sh0norm=sqrt(length(shidx0));

aa=squeeze(min(norbw(2:end,:,:)));
bb=squeeze(max(norbw(2:end,:,:)));
aa=squeeze(min(sfpeak(2:end,:,:)));
bb=squeeze(max(sfpeak(2:end,:,:)));
aa=squeeze(min(psfpeak(2:end,:,:)));
bb=squeeze(max(psfpeak(2:end,:,:)));
aa=squeeze(min(nsfpeak(2:end,:,:)));
bb=squeeze(max(nsfpeak(2:end,:,:)));


figure(3);
clf
plot(aa,bb,'.');
hold on
plot(aa(shidx),bb(shidx),'r.');
plot([0 4],[0 4],'k--');
hold off
axis equal
axis square
[nanmean([abs(aa(shidx)-bb(shidx))]) ...
 nanstd(abs(aa(shidx)-bb(shidx)))./shnorm
 nanmean([abs(aa(shidx0)-bb(shidx0))]) ...
 nanstd(abs(aa(shidx0)-bb(shidx0)))./sh0norm
]

disp('variance across negative sf peak???');
aa=squeeze(std(((nsfbw(2:end,:,:)))));
aa=squeeze(std(log2((nsfpeak(2:end,:,:)))));
[nanmean(aa(shidx))  nanstd(aa(shidx))./shnorm
 nanmean(aa(shidx0)) nanstd(aa(shidx0))./sh0norm]

disp('need to look at trends across a bunch of tuning properties here!');
keyboard






DISPLAYON=0;
starstr={'','*'};
disp('  CELL    PRED XC     PEAK SF     PEAK PSF     PEAK NSF');
for ii=shidx,
   fprintf('%3d. %-7s %5.2f%5.2f  %4.2f %4.2f  %4.2f %4.2f  %4.2f %4.2f\n',...
           ii,res(ii).cellid,predxc(1,[1 aidx+1],ii),...
           sfpeak(aidx,2:3,ii),psfpeak(aidx,2:3,ii),nsfpeak(aidx,2:3,ii));
   if DISPLAYON,
      
      figure(3);
      clf
      
      ormin=min(min(min(nortuning(:,:,:,ii))))*1.1;
      ormax=max(max(max(portuning(:,:,:,ii))))*1.1;
      sfmin=min(min(min(nsftuning(:,:,:,ii))))*1.1;
      sfmax=max(max(max(psftuning(:,:,:,ii))))*1.1;
      
      showkern(permute(cat(3,res(ii).Hset,res(ii).Hset(:,:,1:2)),...
                       [1 3 2]),'pfftgr',[],{},1);
         
      subplot(2,5,1);
      title(sprintf('%s baseline',res(ii).cellid));
      subplot(2,5,6);
      cla; axis off
      kcount=size(res(ii).globaldc,1);
      bootcount=size(res(ii).globaldc,2);
      mdc=mean(res(ii).globaldc(2:end,:),2);
      edc=std(res(ii).globaldc(2:end,:),1,2) .* sqrt(bootcount);
      plot(linspace(0,kcount+1,bootcount)',res(ii).globaldc(2:end,:)','--');
      hold on
      errorbar(mdc,edc.*2,'ko');
      hold off
      title('global dc');
      axis square
      
      subplot(2,5,2);
      title(sprintf('Aa: %5.2f/%5.2f %s',...
                    res(ii).predxc(1:2,2),starstr{(res(ii).pxc(2)<PTHRESH)+1}));
      subplot(2,5,7);
      title(sprintf('AB: %5.2f/%5.2f %s',...
                    res(ii).predxc(1:2,3),starstr{(res(ii).pxc(3)<PTHRESH)+1}));
      subplot(2,5,3);
      title(sprintf('Bb'));
      subplot(2,5,8);
      title(sprintf('ab'));
      
      subplot(2,5,4);
      plot(obins,squeeze(portuning(:,1,2:3,ii)),'LineWidth',2);
      hold on
      plot(obins,squeeze(nortuning(:,1,2:3,ii)),'--','LineWidth',2);
      hold off
      axis([min(obins) max(obins) ormin ormax]);
      title(sprintf('Aa/Bb p0: %.1f/%.1f\nn0: %.1f/%.1f\norbw: %.1f/%.1f',...
                    porpeak(1,2:end,ii),norpeak(1,2:end,ii),orbw(1,2:end,ii)));
      axis square
      
      subplot(2,5,5);
      plot(sfbins,squeeze(psftuning(:,1,2:3,ii)),'LineWidth',2);
      hold on
      plot(sfbins,squeeze(nsftuning(:,1,2:3,ii)),'--','LineWidth',2);
      hold off
      axis([min(sfbins) max(sfbins) sfmin sfmax]);
      title(sprintf('Aa/Bb psf: %.1f/%.1f\nnsf: %.1f/%.1f\nsfbw: %.1f/%.1f',...
                    psfpeak(1,2:end,ii),nsfpeak(1,2:end,ii),sfbw(1,2:end,ii)));
      xlabel('spatial freq');
      axis square
      legend('Aa','Bb');
      
      subplot(2,5,9);
      plot(obins,squeeze(portuning(:,2,2:3,ii)),'LineWidth',2);
      hold on
      plot(obins,squeeze(nortuning(:,2,2:3,ii)),'--','LineWidth',2);
      hold off
      axis([min(res(ii).obins) max(res(ii).obins) ormin ormax]);
      title(sprintf('AB/ab p0: %.1f/%.1f\nn0: %.1f/%.1f\norbw: %.1f/%.1f',...
                    porpeak(2,2:end,ii),norpeak(2,2:end,ii),orbw(2,2:end,ii)));
      axis square
      
      subplot(2,5,10);
      plot(sfbins,squeeze(psftuning(:,2,2:3,ii)),'LineWidth',2);
      hold on
      plot(sfbins,squeeze(nsftuning(:,2,2:3,ii)),'--','LineWidth',2);
      hold off
      axis([min(res(ii).sfbins) max(res(ii).sfbins) sfmin sfmax]);
      title(sprintf('AB/ab psf: %.1f/%.1f\nnsf: %.1f/%.1f\nsfbw: %.1f/%.1f',...
                    psfpeak(2,2:end,ii),nsfpeak(2,2:end,ii),sfbw(2,2:end,ii)));
      xlabel('spatial freq');
      legend('AB','ab');
      axis square
      
      fullpage landscape
      pause
   end
end
fprintf('---. %-7s %5.2f%5.2f  %4.2f+%4.2f  %4.2f+%4.2f  %4.2f+%4.2f\n',...
        'AVG',mean(predxc(1,[1 aidx+1],shidx),3),...
        median(abs(diff(log2(sfpeak(aidx,2:3,shidx)),[],2)),3),...
        std(abs(diff(log2(sfpeak(aidx,2:3,shidx)),[],2)),3)./shnorm,...
        median(abs(diff(log2(psfpeak(aidx,2:3,shidx)),[],2)),3),...
        std(abs(diff(log2(psfpeak(aidx,2:3,shidx)),[],2)),3)./shnorm,...
        median(abs(diff(log2(nsfpeak(aidx,2:3,shidx)),[],2)),3),...
        std(abs(diff(log2(nsfpeak(aidx,2:3,shidx)),[],2)),3)./shnorm)
fprintf('---. %-7s %5.2f%5.2f  %4.2f+%4.2f  %4.2f+%4.2f  %4.2f+%4.2f\n',...
        'NAVG',mean(predxc(1,[1 aidx+1],shidx0),3),...
        median(abs(diff(log2(sfpeak(aidx,2:3,shidx0)),[],2)),3),...
        std(abs(diff(log2(sfpeak(aidx,2:3,shidx0)),[],2)),3)./sh0norm,...
        median(abs(diff(log2(psfpeak(aidx,2:3,shidx0)),[],2)),3),...
        std(abs(diff(log2(psfpeak(aidx,2:3,shidx0)),[],2)),3)./sh0norm,...
        median(abs(diff(log2(nsfpeak(aidx,2:3,shidx0)),[],2)),3),...
        std(abs(diff(log2(nsfpeak(aidx,2:3,shidx0)),[],2)),3)./sh0norm)

disp('need to look at trends across a bunch of tuning properties here!');
keyboard


aidx=1;

aa=squeeze(pamp(aidx,2,:));
bb=squeeze(pamp(aidx,3,:));

aa=squeeze(namp(aidx,2,:));
bb=squeeze(namp(aidx,3,:));

aa=squeeze(porpeak(aidx,2,:));
bb=squeeze(porpeak(aidx,3,:));
bb(find(bb-aa>90))=bb(find(bb-aa>90))-180;
bb(find(aa-bb>90))=bb(find(aa-bb>90))+180;

aa=squeeze(porbw(aidx,2,:));
bb=squeeze(porbw(aidx,3,:));
aa=squeeze(norbw(aidx,2,:));
bb=squeeze(norbw(aidx,3,:));
aa=log2(squeeze(sfpeak(aidx,2,:)));
bb=log2(squeeze(sfpeak(aidx,3,:)));
aa=log2(squeeze(psfpeak(aidx,2,:)));
bb=log2(squeeze(psfpeak(aidx,3,:)));
aa=log2(squeeze(nsfpeak(aidx,2,:)));
bb=log2(squeeze(nsfpeak(aidx,3,:)));
aa=(squeeze(sfbw(aidx,2,:)));
bb=(squeeze(sfbw(aidx,3,:)));
aa=(squeeze(psfbw(aidx,2,:)));
bb=(squeeze(psfbw(aidx,3,:)));
aa=(squeeze(nsfbw(aidx,2,:)));
bb=(squeeze(nsfbw(aidx,3,:)));

figure(3);
clf
plot(aa,bb,'.');
hold on
plot(aa(find(origlocalsig(:,aidx))),bb(find(origlocalsig(:,aidx))),'r.');
%plot([0 180],[0 180],'k--');
plot([0 4],[0 4],'k--');
hold off
axis equal
axis square
[nanmedian([abs(aa(find(origlocalsig(:,aidx)))-bb(find(origlocalsig(:,aidx))))])
 nanmedian([abs(aa(find(~origlocalsig(:,aidx)))-bb(find(~origlocalsig(:,aidx))))])
 nanstd(abs(aa(find(origlocalsig(:,aidx)))-...
            bb(find(origlocalsig(:,aidx)))))./...
 sqrt(sum(origlocalsig(:,aidx)))
 nanstd([abs(aa(find(~origlocalsig(:,aidx)))-...
             bb(find(~origlocalsig(:,aidx))))])./...
 sqrt(sum(~origlocalsig(:,aidx)))
]




ftall=[];
spall=[];
%for ii=find(ftsig(:,2) + spsig(:,2))',
for ii=1:length(res),
      
   if 0,
      sigparms=res(ii).sigparms;
      p0=linspace(res(ii).sigrange(1,1,1),res(ii).sigrange(2,1,1))';
      rft=[sigmoid(sigparms(:,1,1),p0) sigmoid(sigparms(:,2,1),p0)];
      rsp=[sigmoid(sigparms(:,1,2),p0) sigmoid(sigparms(:,2,2),p0)];
   elseif 1,
      p0=res(ii).rset(:,1);
      rft=[p0 hinge(hingeparms(:,1,ii),p0)];
      rsp=[p0 hinge(hingeparms(:,2,ii),p0)];
      
   elseif 1,
      
      p0=res(ii).rset(:,1);
      rft=[p0 p0+res(ii).rset(:,2) p0-res(ii).rset(:,2)];
      rsp=[p0 p0+res(ii).rset(:,3) p0-res(ii).rset(:,3)];
      
   else
      p0=linspace(res(ii).baserange(1),res(ii).baserange(2),100)';
      
      ftgg=mean(res(ii).globalgain(2,:));
      ftdc=mean(res(ii).globaldc(2,:))./60;
      rft=[p0 p0 + p0*ftgg + ftdc  p0 - p0*ftgg - ftdc];
      
      spgg=mean(res(ii).globalgain(3,:));
      spdc=mean(res(ii).globaldc(3,:))./60;
      rsp=[p0 p0 + p0*spgg + spdc  p0 - p0*spgg - spdc];
   end
   
   if 1,
      ftall=cat(3,ftall,rft);
      spall=cat(3,spall,rsp);
      
      
   else
      [ftdc.*60 ftgg spdc.*60 spgg]
      
      figure(1);
      clf
      
      subplot(1,2,1);
      plot(p0,rft);
      title(sprintf('%s ft Aa vs Bb %.2f-%.2f p<%.2f',res(ii).cellid,...
                    res(ii).fullxc([1 3],2),res(ii).fullp(3,2)));
      
      subplot(1,2,2);
      plot(p0,rsp);
      title(sprintf('%s sp AB vs ab %.2f-%.2f p<%.2f',res(ii).cellid,...
                    res(ii).fullxc([1 3],3),res(ii).fullp(3,3)));
      
      pause
   end
end
plot([mean(abs(ftall(:,2,:)),3) mean(abs(spall(:,2,:)),3)])



