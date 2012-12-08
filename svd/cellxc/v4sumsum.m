
% load from kvatunesum(77) and (99) then generate plots for final
% figures

%batch1=99   % 1.0 crf window
%batch2=82   % 1.1 crf window

batch1=77    % currently 1.0 crf window
batch2=99

tmppath='/auto/k5/david/tmp/kvaparms';

resfile1=sprintf('%s/kvatunesum.batch%d.mat',tmppath,batch1);
z1=load(resfile1);
resfile1=sprintf('%s/kvatunesum.batch%d.mat',tmppath,batch2);
z2=load(resfile1);

%
% FIGURE 5/6 STUFF: orientation and spatial frequency tuning summary
%

EXCLUDEBADCELLS=1

% don't discriminate good from bad predictors... easy to explain
%cc1good=ones(size(z1.cc0));
%cc2good=ones(size(z2.cc0));

% good predictors set to 1, bad predictors set to 2
cc1good=(z1.cc0<=z1.ccrand)+1;
cc2good=(z2.cc0<=z2.ccrand)+1;

% just exclude bad predictors by setting them to zero
if EXCLUDEBADCELLS,
   cc1good(cc1good==2)=0;
   cc2good(cc2good==2)=0;
   cc1good(42)=0; % v2 cell
else
   cc1good(cc1good==2)=1;
   cc2good(cc2good==2)=1;
end

%
% orientation tuning summary (FIGURE 5)
%
%scrat=sum(cc2good==1)./sum(cc1good==1)
scrat=1;
bwht=25;
tpkht=40;

orbw1=z1.mcstd;
orbw2=z2.mcstd;

figure(1);
clf
subplot(2,3,1);
data=orbw1;
data(data>=180)=220;
aa=find(cc1good==1 & data<180);
bb=find(cc1good==1 & data>=180);
histfracs(data(aa),data(bb),...
          sprintf('bat%d orbw',batch1),'bandwidth(deg)',[20 220 0 bwht],12);

subplot(2,3,2);
data=orbw2;
data(data>=180)=220;
aa=find(cc2good==1 & data<180);
bb=find(cc2good==1 & data>=180);
histfracs(data(aa),data(bb),...
          sprintf('bat%d orbw',batch2),'bandwidth(deg)',[20 220 0 bwht],12);


m1=mean(z1.twopeakscore0,2);
m1=mean(z1.twopeakscore,2);
s1=std(z1.twopeakscore0,0,2).*sqrt(size(z1.twopeakscore0,2)-5);
m2=mean(z2.twopeakscore0,2);
m2=mean(z2.twopeakscore,2);
s2=std(z2.twopeakscore0,0,2).*sqrt(size(z2.twopeakscore0,2)-5);

[sum(m1(cc1good==1)>s1(cc1good==1)) sum(m2(cc2good==1)>s2(cc2good==1))]

jjset=find(m1>s1 & cc1good==1);
for ii=1:length(jjset),
   fprintf('%s (%d): %.2f  ',z1.cellids{jjset(ii)},...
           jjset(ii),z1.twopeakscore(jjset(ii)));
end
disp('');
jjset=find(m2>s2 & cc2good==1);
for ii=1:length(jjset),
   fprintf('%s (%d): %.2f  ',z2.cellids{jjset(ii)},...
           z2.twopeakscore(jjset(ii)));
end
disp('');

subplot(2,3,4);
data=m1;
%data( data<=0) = -0.25;
histfracs(data(cc1good==1 & m1-s1>0),...
          data(cc1good==1 & m1-s1<=0),...
         sprintf('bat%d bimidx',batch1),'bimodal index',[0 0.75 0 tpkht],12);

subplot(2,3,5);
data=m2;
%data(data<=0) = -0.25;
histfracs(data(cc2good==1 & m2-s2>0),...
          data(cc2good==1 & m2-s2<=0),...
          sprintf('bat%d bimidx',batch2),'bimodal index',[0 0.75 0 tpkht],12);

mcstd1=median(orbw1(cc1good==1));
ecstd1=std(orbw1(cc1good==1))./sqrt(sum(cc1good==1));
mcstd2=median(orbw2(cc2good==1));
ecstd2=std(orbw2(cc2good==1))./sqrt(sum(cc2good==1));

SHOWMED2PK=1; % vs show count of significant cells
if ~SHOWMED2PK,
   m2pk1=nanmean(s1(cc1good==1)<m1(cc1good==1));
   e2pk1=nanstd(s1(cc1good==1)<m1(cc1good==1))./...
         sqrt(sum(~isnan(s1(cc1good==1))));
   m2pk2=nanmean(s2(cc2good==1)<m2(cc2good==1));
   e2pk2=nanstd(s2(cc2good==1)<m2(cc2good==1))./...
         sqrt(sum(~isnan(s2(cc2good==1))));
else
   
   m2pk1=nanmedian(m1(cc1good==1));
   e2pk1=nanstd(m1(cc1good==1))./...
         sqrt(sum(~isnan(m1(cc1good==1))));
   m2pk2=nanmedian(m2(cc2good==1));
   e2pk2=nanstd(m2(cc2good==1))./...
         sqrt(sum(~isnan(m2(cc2good==1))));
end

[porbw]=randttest(orbw1(cc1good==1),orbw2(cc2good==1),500,1)
[p2pk]=randttest(double(s1(cc1good==1)<m1(cc1good==1)),...
                 double(s2(cc2good==1)<m2(cc2good==1)),500,1)

subplot(2,3,3);
ht=errorbar([2 4],[mcstd1 mcstd2],[ecstd1 ecstd2],'k+');
set(ht,'LineWidth',1);
hold on
bar([2 4],[mcstd1 mcstd2]);
hold off
title('median cstd v4 vs v1');
axis square
axis([0 6 0 100]);

subplot(2,3,6);
ht=errorbar([2 4],[m2pk1 m2pk2],[e2pk1 e2pk2],'k+');
set(ht,'LineWidth',1);
hold on
bar([2 4],[m2pk1 m2pk2]);
hold off
title('median 2pk v4 vs v1');
axis square
if SHOWMED2PK,
   axis([0 6 0 0.15]);
else
   axis([0 6 0 0.35]);
end

colormap(gray);
set(gcf,'PaperOrientation','portrait','PaperPosition',[0.75 3 7 5]);

%
% SF tuning summary
%


windowdeg1=cat(1,z1.res.windowdeg);
windowdeg2=cat(1,z2.res.windowdeg);


sfcat1=cc1good;
sfcat2=cc2good;

if 1,
   disp('using gaussian sf stats');
   sfpeak1=z1.sfpeak0;
   sfpeak2=z2.sfpeak0;
   
   sfbw1=z1.sfbw0;
   sfbw2=z2.sfbw0;
   
else
   disp('using center of mass sf stats');
   sfpeak1=z1.sfpeak1;
   sfpeak2=z2.sfpeak1;
   
   sfbw1=z1.sfbw1;
   sfbw2=z2.sfbw1;
end

sfcat1(z1.sfbw1good<1 & cc1good==1)=2;
sfcat2(z2.sfbw1good<1 & cc2good==1)=2;

sfcat1(isnan(sfbw1))=0;
sfcat2(isnan(sfbw2))=0;

figure(2);
clf
if 1,
   disp('using cyc/RF')
   sfpeakrange=[1.0 7.0 0 25];
   
   subplot(2,3,1);
   histfracs(sfpeak1(sfcat1==1),sfpeak1(sfcat1==2),...
            sprintf('bat %d sfpeak0',batch1),'sf (cyc/RF)',sfpeakrange,12);
   
   subplot(2,3,2);
   histfracs(sfpeak2(sfcat2==1),sfpeak2(sfcat2==2),...
             sprintf('bat %d sfpeak0',batch2),'sf (cyc/RF)',sfpeakrange,12);
else
   disp('using cyc/deg')
   sfpeakrange1=[1 6 0 35];
   sfpeakrange2=[1 6 0 35];
   
   subplot(2,3,1);
   histfracs(sfpeak1(sfcat1==1)./windowdeg1(sfcat1==1),...
             sfpeak1(sfcat1==2)./windowdeg1(sfcat1==2),...
             sprintf('bat %d sfpeak0',batch1),'cyc/deg',sfpeakrange1,12);
   
   subplot(2,3,2);
   histfracs(sfpeak2(sfcat2==1)./windowdeg2(sfcat2==1),...
            sfpeak2(sfcat2==2)./windowdeg2(sfcat2==2),...
            sprintf('bat %d sfpeak0',batch2),'sf (cyc/deg)',sfpeakrange2,12);
end

sfbwrange=[0.5 2.5 0 25];

subplot(2,3,4);
histfracs(sfbw1(sfcat1==1),sfbw1(sfcat1==2),...
         sprintf('bat %d sfbw0',batch1),'sfbw (oct)',sfbwrange,12);

subplot(2,3,5);
histfracs(sfbw2(sfcat2==1),sfbw2(sfcat2==2),...
         sprintf('bat %d sfbw0',batch2),'sfbw (oct)',sfbwrange,12);

mpeak1=median(sfpeak1(sfcat1==1)./windowdeg1(sfcat1==1));
epeak1=std(sfpeak1(sfcat1==1)./windowdeg1(sfcat1==1)) ./ ...
       sqrt(sum(sfcat1==1));
mpeak2=median(sfpeak2(sfcat2==1)./windowdeg2(sfcat2==1));
epeak2=std(sfpeak2(sfcat2==1)./windowdeg2(sfcat2==1)) ./ ...
       sqrt(sum(sfcat2==1));
mpeak1rf=median(sfpeak1(sfcat1==1));
epeak1rf=std(sfpeak1(sfcat1==1)) ./ sqrt(sum(sfcat1==1));
mpeak2rf=median(sfpeak2(sfcat2==1));
epeak2rf=std(sfpeak2(sfcat2==1)) ./ sqrt(sum(sfcat2==1));

mbw1=median(sfbw1(sfcat1==1));
ebw1=std(sfbw1(sfcat1==1))./sqrt(z1.cellcount);
mbw2=median(sfbw2(sfcat2==1));
ebw2=std(sfbw2(sfcat2==1))./sqrt(z2.cellcount);
mbw1all=median(sfbw1(sfcat1 >0));
ebw1all=std(sfbw1(sfcat1>0))./sqrt(z1.cellcount);
mbw2all=median(sfbw2(sfcat2 >0));
ebw2all=std(sfbw2(sfcat2>0))./sqrt(z2.cellcount);

[psfpk]=randttest(sfpeak1(cc1good==1),sfpeak2(cc2good==1),500,1)
[psfbw]=randttest(sfbw1(cc1good==1),sfbw2(cc2good==1),500,1)

subplot(2,3,3);
%errorbar([1 2 4 5],[mpeak1rf mpeak2rf mpeak1 mpeak2 ],...
%         [ epeak1rf epeak2rf epeak1 epeak2],'k+');
ht=errorbar([2 4],[mpeak1rf mpeak2rf ],[epeak1rf epeak2rf],'k+');
set(ht,'LineWidth',1);
hold on
bar([2 4],[mpeak1rf mpeak2rf]);
hold off
title(sprintf('median sfpeak bat %d v %d',batch1,batch2));
axis square
axis([0 6 0 3.5]);

subplot(2,3,6);
%errorbar([1 2 4 5],[mbw1 mbw2 mbw1all mbw2all],...
%         [ebw1 ebw2 ebw1all ebw2all],'k+');
ht=errorbar([2 4],[mbw1 mbw2],[ebw1 ebw2],'k+');
set(ht,'LineWidth',1);
hold on
%bar([1 2 4 5],[mbw1 mbw2 mbw1all mbw2all]);
bar([2 4],[mbw1 mbw2]);
hold off
title(sprintf('median sfbw bat %d v %d',batch1,batch2));
axis square
axis([0 6 0 1.5]);

colormap(gray);
set(gcf,'PaperOrientation','portrait','PaperPosition',[0.75 3 7 5]);
drawnow


%
% FIGURE 7 : cart/ncart, synth/natural tuning summary (v4 and v1)
%

% sparseness measured using vinje selectivity metric from J Neurosci
% V4 - [cart pol hyp nat]
% V1 - [cart pol hyp nat]
%[nanmean(z1.sparseness) nanmean(z2.sparseness)]

% V4 - [cart pol+hyp cart+pol+hyp nat]
% V1 - [cart pol+hyp cart+pol+hyp nat] ...
%[nanmean(z1.sparseplus) nanmean(z2.sparseplus)]

% entropy meant to measure mutual information between stimulus and
% response. ie, I(R;S) = H(R) - H(R|S)
% currently, H(R|S)=0 because there's no noise in the response
% is this reasonable? 
% V4 - [cart pol+hyp cart+pol+hyp nat]
% V1 - [cart pol+hyp cart+pol+hyp nat] ...
%[nanmean(z1.entropy) nanmean(z2.entropy)]      


% set to -100 so that all cells are "good"
if ~EXCLUDEBADCELLS,
   z1.ccrand(:)=-100;
   z2.ccrand(:)=-100;
end

dg1=z1.cc0>z1.ccrand;
dg2=z2.cc0>z2.ccrand;
goodcount1=sum(dg1);
goodcount2=sum(dg2);

if 0,
   disp('special short best class analysis for cosyne poster');
   
   % special display for cosyne poster
   z1dnorm=cat(1,z1.res.cartvscc);
   z2dnorm=cat(1,z2.res.cartvscc);
   
   z1better=z1dnorm(dg1,[1 2 4])>repmat(z1dnorm(dg1,3),[1 3]);
   z2better=z2dnorm(dg2,[1 2 4])>repmat(z2dnorm(dg2,3),[1 3]);
   
   figure(3);
   clf
   
   subplot(1,2,1);
   bset=[sum(z1better==1); sum(z1better==0)]';
   bar(bset);
   title('v1');
   
   subplot(1,2,2);
   bset=[sum(z2better==1); sum(z2better==0)]';
   bar(bset);
   title('v4');
   
   set(gcf,'PaperPosition',[1 4 6 2.5]);
   
elseif 1,
   
   z1dnorm=[(cat(1,z1.res.cartvscc)) cat(1,z1.res.synthvs)];
   z1dnormerr=[(cat(1,z1.res.cartvsccerr)) cat(1,z1.res.synthvserr)];
   z2dnorm=[(cat(1,z2.res.cartvscc)) cat(1,z2.res.synthvs)];
   z2dnormerr=[(cat(1,z2.res.cartvsccerr)) cat(1,z2.res.synthvserr)];
   
   z1dmax=z1dnorm(:,1:3);
   z1dmaxerr=z1dnormerr(:,1:3);
   z2dmax=z2dnorm(:,1:3);
   z2dmaxerr=z2dnormerr(:,1:3);
   
   zcount=size(z1dmax,2);
   nf=sqrt(2);
   z1bestcat=zeros(size(z1dmax,1),1);
   z2bestcat=zeros(size(z2dmax,1),1);
   for ii=1:size(z1dmax,2),
      z1bestcat(z1dmax(:,ii)==max(z1dmax,[],2))=-ii;
      z1bestcat(z1dmax(:,ii)-z1dmaxerr(:,ii)>...
                max(z1dmax(:,[1:ii-1 ii+1:zcount]) + ...
                    z1dmaxerr(:,[1:ii-1 ii+1:zcount]),[],2))=ii;
      z2bestcat(z2dmax(:,ii)==max(z2dmax,[],2))=-ii;
      z2bestcat(z2dmax(:,ii)-z2dmaxerr(:,ii)>...
                max(z2dmax(:,[1:ii-1 ii+1:zcount]) + ...
                    z2dmaxerr(:,[1:ii-1 ii+1:zcount]),[],2))=ii;
   end
   
   mz1cart=nanmedian(z1dnorm(dg1,1:zcount));
   ez1cart=nanstd(z1dnorm(dg1,1:zcount))./sqrt(length(dg1));
   
   mz2cart=nanmedian(z2dnorm(dg1,1:zcount));
   ez2cart=nanstd(z2dnorm(dg1,1:zcount))./sqrt(length(dg1));
   
   figure(3);
   clf
   
   histht=25;
   
   for ii=1:zcount,
      subplot(zcount,3,ii*3-2);
      histfracs(z1dnorm(abs(z1bestcat)==ii & dg1,ii),...
                z1dnorm(abs(z1bestcat)~=ii & dg1,ii),...
                sprintf('bat%d cat %d sel',batch1,ii),...
                'norm pred',[0.05 0.95 0 histht]);
      
      subplot(zcount,3,ii*3-1);
      histfracs(z2dnorm(abs(z2bestcat)==ii & dg2,ii),...
                z2dnorm(abs(z2bestcat)~=ii & dg2,ii),...
                sprintf('bat%d cat %d sel',batch2,ii),...
                'norm pred',[0.05 0.95 0 histht]);
   end
   
   v1set1=find(abs(z1bestcat)==1 & dg1);
   v1set2=find(abs(z1bestcat)==2 & dg1);
   v1set3=find(abs(z1bestcat)==3 & dg1);
   set1=find(abs(z2bestcat)==1 & dg2);
   set2=find(abs(z2bestcat)==2 & dg2);
   set3=find(abs(z2bestcat)==3 & dg2);
   orwset=[median(orbw2(set1)) ...
           median(orbw2(set2)) ...
           median(orbw2(set3))];
   tpkset=[median(z2.twopeakscore(set1)) ...
           median(z2.twopeakscore(set2)) ...
           median(z2.twopeakscore(set3))];
   sfwset=[median(sfbw2(set1)) ...
           median(sfbw2(set2)) ...
           median(sfbw2(set3))];
   sfpset=[median(sfpeak2(set1)) ...
           median(sfpeak2(set2)) ...
           median(sfpeak2(set3))];
   orwsete=[std(orbw2(set1))./sqrt(length(set1)) ...
           std(orbw2(set2))./sqrt(length(set2)) ...
           std(orbw2(set3))./sqrt(length(set3))];
   tpksete=[std(z2.twopeakscore(set1))./sqrt(length(set1)) ...
           std(z2.twopeakscore(set2))./sqrt(length(set2)) ...
           std(z2.twopeakscore(set3))./sqrt(length(set3))];
   sfwsete=[std(sfbw2(set1))./sqrt(length(set1)) ...
           std(sfbw2(set2))./sqrt(length(set2)) ...
           std(sfbw2(set3))./sqrt(length(set3))];
   sfpsete=[std(sfpeak2(set1))./sqrt(length(set1)) ...
           std(sfpeak2(set2))./sqrt(length(set2)) ...
           std(sfpeak2(set3))./sqrt(length(set3))];
   
   
   v1cat=z1bestcat(find(dg1));
   v4cat=z2bestcat(find(dg2));
   v1frac=zeros(20,3);
   v4frac=zeros(20,3);
   stepsize1=goodcount1./20;
   stepsize2=goodcount2./20;
   for jnidx=1:20,
      useidx=[1:round((jnidx-1)*stepsize1) ...
              round(jnidx*stepsize1+1):goodcount1];
      
      set1=find(abs(v1cat(useidx))==1);
      set2=find(abs(v1cat(useidx))==2);
      set3=find(abs(v1cat(useidx))==3);
      v1frac(jnidx,:)=[length(set1) length(set2) length(set3)]./length(useidx);
      useidx=[1:round((jnidx-1)*stepsize2) ...
              round(jnidx*stepsize2+1):goodcount2];
      
      set1=find(abs(v4cat(useidx))==1);
      set2=find(abs(v4cat(useidx))==2);
      set3=find(abs(v4cat(useidx))==3);
      v4frac(jnidx,:)=[length(set1) length(set2) length(set3)]./length(useidx);
   end
   v1mm=mean(v1frac);
   v1ss=std(v1frac).*sqrt(20);
   v4mm=mean(v4frac);
   v4ss=std(v4frac).*sqrt(20);
   
   subplot(3,2,1);
   errorbar(1:3,v4mm,v4ss,'k+');
   hold on
   bar(v4mm);
   hold off
   title('v4 pref class');
   axis square
   
   subplot(3,2,2);
   errorbar(1:3,v1mm,v1ss,'k+');
   hold on
   bar(v1mm);
   hold off
   title('v1 pref class');
   axis square
   
   subplot(3,2,3);
   errorbar(1:3,orwset,orwsete,'k+');
   hold on
   bar(orwset);
   hold off
   title('median orw: car/nc/nat');
   axis square
   
   subplot(3,2,4);
   errorbar(1:3,tpkset,tpksete,'k+');
   hold on
   bar(tpkset);
   hold off
   title('median tpk: car/nc/nat');
   axis square
   
   subplot(3,2,5);
   errorbar(1:3,sfpset,sfpsete,'k+');
   hold on
   bar(sfpset);
   hold off
   title('median sfp: car/nc/nat');
   axis square
   
   subplot(3,2,6);
   errorbar(1:3,sfwset,sfwsete,'k+');
   hold on
   bar(sfwset);
   hold off
   title('median sfw: car/nc/nat');
   axis square
   
   colormap(gray);
   set(gcf,'PaperOrientation','portrait','PaperPosition',[2 2 3 6]);
   
elseif 1,
   
   %old: plot out histograms for each area/stimulus class
   
   z1dnorm=[(cat(1,z1.res.cartvscc)) cat(1,z1.res.synthvs)];
   z1dnormerr=[(cat(1,z1.res.cartvsccerr)) cat(1,z1.res.synthvserr)];
   z2dnorm=[(cat(1,z2.res.cartvscc)) cat(1,z2.res.synthvs)];
   z2dnormerr=[(cat(1,z2.res.cartvsccerr)) cat(1,z2.res.synthvserr)];
   
   if 0,
      z1dmax=[cat(1,z1.res.cartvs) cat(1,z1.res.synthvs)];
      z1dmaxerr=[cat(1,z1.res.cartvserr) cat(1,z1.res.synthvserr)];
      z2dmax=[cat(1,z2.res.cartvs) cat(1,z2.res.synthvs)];
      z2dmaxerr=[cat(1,z2.res.cartvserr) cat(1,z2.res.synthvserr)];
   elseif 0
      % compare cart/ncart/edstim
      z1dmax=z1dnorm(:,1:4);
      z1dmaxerr=z1dnormerr(:,1:4);
      z2dmax=z2dnorm(:,1:4);
      z2dmaxerr=z2dnormerr(:,1:4);
   else
      % stick to old comparison: cart/ncart/nat
      
      z1dmax=z1dnorm(:,1:3);
      z1dmaxerr=z1dnormerr(:,1:3);
      z2dmax=z2dnorm(:,1:3);
      z2dmaxerr=z2dnormerr(:,1:3);
   end
   
   zcount=size(z1dmax,2);
   nf=sqrt(2);
   z1bestcat=zeros(size(z1dmax,1),1);
   z2bestcat=zeros(size(z2dmax,1),1);
   for ii=1:size(z1dmax,2),
      z1bestcat(z1dmax(:,ii)==max(z1dmax,[],2))=-ii;
      z1bestcat(z1dmax(:,ii)-z1dmaxerr(:,ii)>...
                max(z1dmax(:,[1:ii-1 ii+1:zcount]) + ...
                    z1dmaxerr(:,[1:ii-1 ii+1:zcount]),[],2))=ii;
      z2bestcat(z2dmax(:,ii)==max(z2dmax,[],2))=-ii;
      z2bestcat(z2dmax(:,ii)-z2dmaxerr(:,ii)>...
                max(z2dmax(:,[1:ii-1 ii+1:zcount]) + ...
                    z2dmaxerr(:,[1:ii-1 ii+1:zcount]),[],2))=ii;
   end
   
   mz1cart=nanmedian(z1dnorm(dg1,1:zcount));
   ez1cart=nanstd(z1dnorm(dg1,1:zcount))./sqrt(length(dg1));
   
   mz2cart=nanmedian(z2dnorm(dg1,1:zcount));
   ez2cart=nanstd(z2dnorm(dg1,1:zcount))./sqrt(length(dg1));
   
   figure(3);
   clf
   
   histht=25;
   
   for ii=1:zcount,
      subplot(zcount,3,ii*3-2);
      histfracs(z1dnorm(abs(z1bestcat)==ii & dg1,ii),...
                z1dnorm(abs(z1bestcat)~=ii & dg1,ii),...
                sprintf('bat%d cat %d sel',batch1,ii),...
                'norm pred',[0.05 0.95 0 histht]);
      
      subplot(zcount,3,ii*3-1);
      histfracs(z2dnorm(abs(z2bestcat)==ii & dg2,ii),...
                z2dnorm(abs(z2bestcat)~=ii & dg2,ii),...
                sprintf('bat%d cat %d sel',batch2,ii),...
                'norm pred',[0.05 0.95 0 histht]);
   end
   
   subplot(3,3,3);
   nccat=ones(size(orbw2)).*3;
   nccat(find(z2dnorm(:,5)>z2dnorm(:,6)))=1;
   nccat(find(nccat==1 & z2.gmax(:,4)>z2.gmax(:,3) & ...
              z2.gmax(:,4)>z2.gmax(:,2)))=2;
   
   set1=find(abs(z2bestcat)==1 & dg2);
   set2=find(abs(z2bestcat)==2 & dg2);
   set3=find(abs(z2bestcat)==3 & dg2);
   orwset=[median(orbw2(set1)) ...
           median(orbw2(set2)) ...
           median(orbw2(set3))];
   tpkset=[median(z2.twopeakscore(set1)) ...
           median(z2.twopeakscore(set2)) ...
           median(z2.twopeakscore(set3))];
   sfwset=[median(sfbw2(set1)) ...
           median(sfbw2(set2)) ...
           median(sfbw2(set3))];
   sfpset=[median(sfpeak2(set1)) ...
           median(sfpeak2(set2)) ...
           median(sfpeak2(set3))];
   orwsete=[std(orbw2(set1))./sqrt(length(set1)) ...
           std(orbw2(set2))./sqrt(length(set2)) ...
           std(orbw2(set3))./sqrt(length(set3))];
   tpksete=[std(z2.twopeakscore(set1))./sqrt(length(set1)) ...
           std(z2.twopeakscore(set2))./sqrt(length(set2)) ...
           std(z2.twopeakscore(set3))./sqrt(length(set3))];
   sfwsete=[std(sfbw2(set1))./sqrt(length(set1)) ...
           std(sfbw2(set2))./sqrt(length(set2)) ...
           std(sfbw2(set3))./sqrt(length(set3))];
   sfpsete=[std(sfpeak2(set1))./sqrt(length(set1)) ...
           std(sfpeak2(set2))./sqrt(length(set2)) ...
           std(sfpeak2(set3))./sqrt(length(set3))];
   
   subplot(4,3,3);
   errorbar(1:3,orwset,orwsete,'k+');
   hold on
   bar(orwset);
   hold off
   title('median orw: car/nc/nat');
   axis square
   
   subplot(4,3,6);
   errorbar(1:3,tpkset,tpksete,'k+');
   hold on
   bar(tpkset);
   hold off
   title('median tpk: car/nc/nat');
   axis square
   
   subplot(4,3,9);
   errorbar(1:3,sfpset,sfpsete,'k+');
   hold on
   bar(sfpset);
   hold off
   title('median sfp: car/nc/nat');
   axis square
   
   subplot(4,3,12);
   errorbar(1:3,sfwset,sfwsete,'k+');
   hold on
   bar(sfwset);
   hold off
   title('median sfw: car/nc/nat');
   axis square
   
   if 0
      subplot(3,3,6);
      bar([z1catcount2 z2catcount2]','stacked');
      title(sprintf('pref cart/nat bat %d v %d',batch1,batch2));
      axis square
   end
   
   colormap(gray);
   set(gcf,'PaperOrientation','portrait','PaperPosition',[0.75 2 7 7]);
   
elseif 1,
   
   % differences 
   
   if 1,
      z1dmax=[cat(1,z1.res.cartvscc) cat(1,z1.res.synthvs)];
      z1dmaxerr=[cat(1,z1.res.cartvsccerr) cat(1,z1.res.synthvserr)];
      z2dmax=[cat(1,z2.res.cartvscc) cat(1,z2.res.synthvs)];
      z2dmaxerr=[cat(1,z2.res.cartvsccerr) cat(1,z2.res.synthvserr)];
   else
      z1dmax=[cat(1,z1.res.cartvs) cat(1,z1.res.synthvs)];
      z1dmaxerr=[cat(1,z1.res.cartvserr) cat(1,z1.res.synthvserr)];
      z2dmax=[cat(1,z2.res.cartvs) cat(1,z2.res.synthvs)];
      z2dmaxerr=[cat(1,z2.res.cartvserr) cat(1,z2.res.synthvserr)];
   end
   z1dmin=[cat(1,z1.res.cartmin) cat(1,z1.res.polmin) ...
           cat(1,z1.res.hypmin) cat(1,z1.res.natmin) ];
   
   z2dmin=[cat(1,z2.res.cartmin) cat(1,z2.res.polmin) ...
           cat(1,z2.res.hypmin) cat(1,z2.res.natmin) ];
   
   z1bestcat1=-ones(size(z1dmax,1),1);
   z1bestcat1(z1dmax(:,2)>z1dmax(:,1))=-2;
   z1bestcat1(z1dmax(:,2)+z1dmaxerr(:,2)<z1dmax(:,1)-z1dmaxerr(:,1))=1;
   z1bestcat1(z1dmax(:,2)-z1dmaxerr(:,2)>z1dmax(:,1)+z1dmaxerr(:,1))=2;
   z1bestcat1(z1.cc0<=z1.ccrand)=0;
   z1catcount1=hist(z1bestcat1,-2:2);
   z1catcount1=z1catcount1([4 5; 2 1]);
   z1catcount1=z1catcount1./sum(z1catcount1(:));
   
   z1bestcat2=-ones(size(z1dmax,1),1);
   z1bestcat2(z1dmax(:,3)>z1dmax(:,1))=-2;
   z1bestcat2(z1dmax(:,3)+z1dmaxerr(:,3)<z1dmax(:,1)-z1dmaxerr(:,1))=1;
   z1bestcat2(z1dmax(:,3)-z1dmaxerr(:,3)>z1dmax(:,1)+z1dmaxerr(:,1))=2;
   z1bestcat2(z1.cc0<=z1.ccrand)=0;
   z1catcount2=hist(z1bestcat2,-2:2);
   z1catcount2=z1catcount2([4 5; 2 1]);
   z1catcount2=z1catcount2./sum(z1catcount2(:));

   z1bestcat3=-ones(size(z1dmax,1),1);
   z1bestcat3(z1dmax(:,3)>z1dmax(:,2))=-2;
   z1bestcat3(z1dmax(:,3)+z1dmaxerr(:,3)<z1dmax(:,2)-z1dmaxerr(:,2))=1;
   z1bestcat3(z1dmax(:,3)-z1dmaxerr(:,3)>z1dmax(:,2)+z1dmaxerr(:,2))=2;
   z1bestcat3(z1.cc0<=z1.ccrand)=0;
   z1catcount3=hist(z1bestcat3,-2:2);
   z1catcount3=z1catcount3([4 5; 2 1]);
   z1catcount3=z1catcount3./sum(z1catcount3(:));
   
   z1bestcat4=-ones(size(z1dmax,1),1);
   z1bestcat4(z1dmax(:,5)>z1dmax(:,4))=-2;
   z1bestcat4(z1dmax(:,5)+z1dmaxerr(:,5)<z1dmax(:,4)-z1dmaxerr(:,4))=1;
   z1bestcat4(z1dmax(:,5)-z1dmaxerr(:,5)>z1dmax(:,4)+z1dmaxerr(:,4))=2;
   z1bestcat4(z1.cc0<=z1.ccrand)=0;
   z1catcount4=hist(z1bestcat4,-2:2);
   z1catcount4=z1catcount4([4 5; 2 1]);
   z1catcount4=z1catcount4./sum(z1catcount4(:));
   
   z2bestcat1=-ones(size(z2dmax,1),1);
   z2bestcat1(z2dmax(:,2)>z2dmax(:,1))=-2;
   z2bestcat1(z2dmax(:,2)+z2dmaxerr(:,2)<z2dmax(:,1)-z2dmaxerr(:,1))=1;
   z2bestcat1(z2dmax(:,2)-z2dmaxerr(:,2)>z2dmax(:,1)+z2dmaxerr(:,1))=2;
   z2bestcat1(z2.cc0<=z2.ccrand)=0;
   z2catcount1=hist(z2bestcat1,-2:2);
   z2catcount1=z2catcount1([4 5; 2 1]);
   z2catcount1=z2catcount1./sum(z2catcount1(:));
   
   z2bestcat2=-ones(size(z2dmax,1),1);
   z2bestcat2(z2dmax(:,3)>z2dmax(:,1))=-2;
   z2bestcat2(z2dmax(:,3)+z2dmaxerr(:,3)<z2dmax(:,1)-z2dmaxerr(:,1))=1;
   z2bestcat2(z2dmax(:,3)-z2dmaxerr(:,3)>z2dmax(:,1)+z2dmaxerr(:,1))=2;
   z2bestcat2(z2.cc0<=z2.ccrand)=0;
   z2catcount2=hist(z2bestcat2,-2:2);
   z2catcount2=z2catcount2([4 5; 2 1]);
   z2catcount2=z2catcount2./sum(z2catcount2(:));
   
   z2bestcat3=-ones(size(z2dmax,1),1);
   z2bestcat3(z2dmax(:,3)>z2dmax(:,2))=-2;
   z2bestcat3(z2dmax(:,3)+z2dmaxerr(:,3)<z2dmax(:,2)-z2dmaxerr(:,2))=1;
   z2bestcat3(z2dmax(:,3)-z2dmaxerr(:,3)>z2dmax(:,2)+z2dmaxerr(:,2))=2;
   z2bestcat3(z2.cc0<=z2.ccrand)=0;
   z2catcount3=hist(z2bestcat3,-2:2);
   z2catcount3=z2catcount3([4 5; 2 1]);
   z2catcount3=z2catcount3./sum(z2catcount3(:));
   
   z2bestcat4=-ones(size(z2dmax,1),1);
   z2bestcat4(z2dmax(:,5)>z2dmax(:,4))=-2;
   z2bestcat4(z2dmax(:,5)+z2dmaxerr(:,5)<z2dmax(:,4)-z2dmaxerr(:,4))=1;
   z2bestcat4(z2dmax(:,5)-z2dmaxerr(:,5)>z2dmax(:,4)+z2dmaxerr(:,4))=2;
   z2bestcat4(z1.cc0<=z1.ccrand)=0;
   z2catcount4=hist(z2bestcat4,-2:2);
   z2catcount4=z2catcount4([4 5; 2 1]);
   z2catcount4=z2catcount4./sum(z2catcount4(:));
   
   z1meanresp=cat(1,z1.res.meanactresp);
   z2meanresp=cat(1,z2.res.meanactresp);
   
   z1dmm=repmat(z1meanresp,[1 size(z1dmax,2)]);
   z2dmm=repmat(z2meanresp,[1 size(z2dmax,2)]);
   
   %z1dmod=(z1dmax-z1dmin)./(z1dmm);
   %z2dmod=(z2dmax-z2dmin)./(z2dmm);
   
   %z1dmod=(z1dmax)./abs(z1dmm);
   %z2dmod=(z2dmax)./abs(z2dmm);
   
   z1dmod=z1dmax;
   z2dmod=z2dmax;
   
   %z1dmod(z1dmod>3 | z1dmod<0)=nan;
   %z2dmod(z2dmod>3 | z2dmod<0)=nan;
   
   %z1dmod(z1dmod<0)=nan;
   %z2dmod(z2dmod<0)=nan;
   
   %z1dmod(z1.cc0<z1.ccrand,:)=nan;
   %z2dmod(z2.cc0<z2.ccrand,:)=nan;
   
   mz1cart=nanmedian(z1dmod(find(cc1good==1),1));
   ez1cart=nanstd(z1dmod((cc1good==1),1))./sqrt(sum(cc1good==1));
   mz1ncart=nanmedian(z1dmod((cc1good==1),2));
   ez1ncart=nanstd(z1dmod((cc1good==1),2))./sqrt(sum(cc1good==1));
   mz1nat=nanmedian(z1dmod((cc1good==1),3));
   ez1nat=nanstd(z1dmod((cc1good==1),3))./sqrt(sum(cc1good==1));
   
   mz2cart=nanmedian(z2dmod((cc2good==1),1));
   ez2cart=nanstd(z2dmod((cc2good==1),1))./sqrt(sum(cc1good==1));
   mz2ncart=nanmedian(z2dmod((cc2good==1),2));
   ez2ncart=nanstd(z2dmod((cc2good==1),2))./sqrt(sum(cc1good==1));
   mz2nat=nanmedian(z2dmod((cc2good==1),3));
   ez2nat=nanstd(z2dmod((cc2good==1),3))./sqrt(sum(cc1good==1));
   
   figure(3);
   clf
   
   histht=40;
   
   subplot(3,3,1);
   data1=(z1dmod(z1bestcat1>0 & dg1,2)-z1dmod(z1bestcat1>0 & dg1,1));
   data2=(z1dmod(z1bestcat1<=0 & dg1,2)-z1dmod(z1bestcat1<=0 & dg1,1));
   histfracs(data1,data2,sprintf('bat %d cart/ncart',batch1),...
             'frac diff',[-0.5 0.5 0 histht]);
   hold on; plot([0 0],[0 histht],'k--'); hold off;
   
   subplot(3,3,2);
   data1=(z2dmod(z2bestcat1>0 & dg2,2)-z2dmod(z2bestcat1>0 & dg2,1));
   data2=(z2dmod(z2bestcat1<=0 & dg2,2)-z2dmod(z2bestcat1<=0 & dg2,1));
   histfracs(data1,data2,sprintf('bat %d cart/ncart',batch2),...
             'frac diff',[-0.5 0.5 0 histht]);
   hold on; plot([0 0],[0 histht],'k--'); hold off;
   
   subplot(3,3,3);
   bar([z1catcount1 z2catcount1]','stacked');
   title(sprintf('pref cart/ncart bat %d v %d',batch1,batch2));
   axis square
   
   subplot(3,3,4);
   data1=(z1dmod(z1bestcat2>0 & dg1,3)-z1dmod(z1bestcat2>0 & dg1,1));
   data2=(z1dmod(z1bestcat2<=0 & dg1,3)-z1dmod(z1bestcat2<=0 & dg1,1));
   histfracs(data1,data2,sprintf('bat %d cart/nat',batch1),...
            'frac diff',[-0.5 0.5 0 histht]);
   hold on; plot([0 0],[0 histht],'k--'); hold off;
   
   subplot(3,3,5);
   data1=(z2dmod(z2bestcat2>0 & dg2,3)-z2dmod(z2bestcat2>0 & dg2,1));
   data2=(z2dmod(z2bestcat2<=0 & dg2,3)-z2dmod(z2bestcat2<=0 & dg2,1));
   histfracs(data1,data2,sprintf('bat %d cart/nat',batch2),...
             'frac diff',[-0.5 0.5 0 histht]);
   hold on; plot([0 0],[0 histht],'k--'); hold off;
   
   subplot(3,3,6);
   bar([z1catcount2 z2catcount2]','stacked');
   title(sprintf('pref cart/nat bat %d v %d',batch1,batch2));
   axis square
   
   if 1,
      subplot(3,3,7);
      data1=(z1dmod(z1bestcat3>0 & dg1,3)-z1dmod(z1bestcat3>0 & dg1,2));
      data2=(z1dmod(z1bestcat3<=0 & dg1,3)-z1dmod(z1bestcat3<=0 & dg1,2));
      histfracs(data1,data2,sprintf('bat %d ncart/nat',batch1),...
                'frac diff',[-0.5 0.5 0 histht]);
      hold on; plot([0 0],[0 histht],'k--'); hold off;
      
      subplot(3,3,8);
      data1=(z2dmod(z2bestcat3>0 & dg2,3)-z2dmod(z2bestcat3>0 & dg2,2));
      data2=(z2dmod(z2bestcat3<=0 & dg2,3)-z2dmod(z2bestcat3<=0 & dg2,2));
      histfracs(data1,data2,sprintf('bat %d ncart/nat',batch2),...
                'frac diff',[-0.5 0.5 0 histht]);
      hold on; plot([0 0],[0 histht],'k--'); hold off;
      
      subplot(3,3,9);
      bar([z1catcount3 z2catcount3]','stacked');
      title(sprintf('pref ncart/nat bat %d v %d',batch1,batch2));
      axis square
   else
      % do synth vs nat
      subplot(3,3,7);
      data1=(z1dmod(z1bestcat3>0 & dg1,5)-z1dmod(z1bestcat3>0 & dg1,4));
      data2=(z1dmod(z1bestcat3<=0 & dg1,5)-z1dmod(z1bestcat3<=0 & dg1,4));
      histcomp(data1,data2,'z2 ncart','nat','pct diff',[-0.5 0.5 0 histht]);
      hold on; plot([0 0],[0 histht],'k--'); hold off;
      
      subplot(3,3,8);
      data1=(z2dmod(z2bestcat3>0 & dg2,5)-z2dmod(z2bestcat3>0 & dg2,4));
      data2=(z2dmod(z2bestcat3<=0 & dg2,5)-z2dmod(z2bestcat3<=0 & dg2,4));
      histcomp(data1,data2,'z1 ncart','nat','pct diff',[-0.5 0.5 0 histht]);
      hold on; plot([0 0],[0 histht],'k--'); hold off;
      
      subplot(3,3,9);
      bar([z1catcount4 z2catcount4]','stacked');
      title(sprintf('pref synth/nat bat %d v %d',batch1,batch2));
      axis square
   end
   
   subplot(3,3,9);
   bar([z1catcount4 z2catcount4]','stacked');
   title(sprintf('pref synth/nat bat %d v %d',batch1,batch2));
   axis square
   
   
   colormap(gray);
   set(gcf,'PaperOrientation','portrait','PaperPosition',[0.75 2 7 7]);
   
   keyboard
else
   v1stat=z1.entropy;
   v4stat=z2.entropy;
   statname='entropy';
   %v1stat=z1.sparseplus;
   %v4stat=z2.sparseplus;
   %statname='sparseplus';
   
   v1stat(z1.cc0<z1.ccrand,:)=nan;
   v4stat(z2.cc0<z2.ccrand,:)=nan;
   
   mv1cart=nanmean(v1stat(:,1));
   ev1cart=nanstd(v1stat(:,1))./sqrt(z1.cellcount);
   mv1ncart=nanmean(v1stat(:,2));
   ev1ncart=nanstd(v1stat(:,2))./sqrt(z1.cellcount);
   mv1nat=nanmean(v1stat(:,4));
   ev1nat=nanstd(v1stat(:,4))./sqrt(z1.cellcount);
   mv4cart=nanmean(v4stat(:,1));
   ev4cart=nanstd(v4stat(:,1))./sqrt(z2.cellcount);
   mv4ncart=nanmean(v4stat(:,2));
   ev4ncart=nanstd(v4stat(:,2))./sqrt(z1.cellcount);
   mv4nat=nanmean(v4stat(:,4));
   ev4nat=nanstd(v4stat(:,4))./sqrt(z2.cellcount);
   
   figure(3);
   clf
   subplot(2,3,1);
   plotcomp(v1stat(:,1),v1stat(:,2),['v1 cart ',statname],'v1 Ncart');
   subplot(2,3,2);
   plotcomp(v4stat(:,1),v4stat(:,2),'v4 cart','v4 Ncart');
   subplot(2,3,3);
   errorbar([mv1cart mv1ncart mv4cart mv4ncart],...
            [ev1cart ev1ncart ev4cart ev4ncart],'k+');
   hold on
   bar([mv1cart mv1ncart mv4cart mv4ncart]);
   hold off
   title('mean cart/ncart v1 vs v4');
   
   subplot(2,3,4);
   plotcomp(v1stat(:,1),v1stat(:,4),['v1 cart ',statname],'v1 nat');
   subplot(2,3,5);
   plotcomp(v4stat(:,1),v4stat(:,4),'v4 cart','v4 nat');
   subplot(2,3,6);
   errorbar([mv1cart mv1nat mv4cart mv4nat],[ev1cart ev1nat ev4cart ev4nat],'k+');
   hold on
   bar([mv1cart mv1nat mv4cart mv4nat]);
   hold off
   title('mean cart/nat v1 vs v4');
end

%
% ED STIM STUFF
%

edparms1=cat(1,z1.res.edparms);
edparms2=cat(1,z2.res.edparms);

edccparms1=cat(1,z1.res.edparms);
edccparms2=cat(1,z2.res.edparms);

for ii=1:length(z1.res),
   edccparms1(ii,:)=z1.res(ii).edccparms(1,:);
end
for ii=1:length(z2.res),
   if length(z2.res(ii).edccparms)>0,
      edccparms2(ii,:)=z2.res(ii).edccparms(1,:);
   end
end
%edparms1=cat(1,z1.res.edccparms);
%edparms2=cat(1,z2.res.edccparms);

edmatch=cat(4,z2.res.edmatch);
es=size(edmatch);

angmax=zeros(es(4),1);
for ii=1:es(4),
   tset=edmatch(:,:,:,ii);
   %tset(:,4,:)=0;
   %tset(:,:,2:3:end)=0;
   %tset(:,:,3:3:end)=0;
   tt=min(find(tset(:)==max(tset(:))));
   [e1,e2,e3]=ind2sub(es,tt);
   
   angmax(ii)=e2;
end

m2p=zeros(max(edparms2(:,2)),2);
mbw=zeros(max(edparms2(:,2)),2);
mn=zeros(max(edparms2(:,2)),2);
e2p=zeros(max(edparms2(:,2)),2);
ebw=zeros(max(edparms2(:,2)),2);
for ii=1:size(m2p,1),
   
   % 1st column is v1, 2nd is v4
   
   if length(find(edparms1(:,2)==ii))>0
      m2p(ii,1)=median(z1.twopeakscore(find(edparms1(:,2)==ii & dg1)));
      mbw(ii,1)=median(orbw1(find(edparms1(:,2)==ii & dg1)));
      mn(ii,1)=sum(edparms1(:,2)==ii & dg1);
   end
   %uidx=find(edparms2(:,2)==ii & dg2);
   uidx=find(angmax==ii & dg2);
   if length(uidx)>0
      m2p(ii,2)=median(z2.twopeakscore(uidx));
      mbw(ii,2)=median(orbw2(uidx));
      mn(ii,2)=length(uidx);
      e2p(ii,2)=std(z2.twopeakscore(uidx))./sqrt(length(uidx));
      ebw(ii,2)=std(orbw2(uidx))./sqrt(length(uidx));
   end
end

uidx=find(angmax<4 & dg2);
sharpbmi=median(z2.twopeakscore(uidx));
sharpbmierr=std(z2.twopeakscore(uidx))./sqrt(length(uidx));
uidx=find(angmax>4 & dg2);
roundbmi=median(z2.twopeakscore(uidx));
roundbmierr=std(z2.twopeakscore(uidx))./sqrt(length(uidx));

uidx=find(ismember(angmax,[1 7]) & dg2);
narrowbw=median(orbw2(uidx));
narrowbwerr=std(orbw2(uidx))./sqrt(length(uidx));
uidx=find(ismember(angmax,[3 5]) & dg2);
widebw=median(orbw2(uidx));
widebwerr=std(orbw2(uidx))./sqrt(length(uidx));


figure(4);
clf

subplot(3,1,1);

ht=errorbar(mbw(:,2),ebw(:,2),'k-');
set(ht,'LineWidth',2);
ylabel('med or bw');
axis([0.5 length(m2p)+0.5 0 200]);

subplot(3,1,2);

ht=errorbar(m2p(:,2),e2p(:,2),'k-');
set(ht,'LineWidth',2);
hold on
%plot([1 length(m2p)],median(z2.twopeakscore(dg2)).*[1 1],'--');
%errorbar(mbw(:,2)./500,ebw(:,2)./500,'r');
%plot(mn(:,2)./sum(mn(:,2)),'g');
hold off
ylabel('med bimod');
axis([0.5 length(m2p)+0.5 0 0.2]);

ffb=loadimfile('/auto/k5/david/sample/stim/edstim.imsm');
ffs=loadimfile('/auto/k5/david/sample/stim/edstim.imsm',0,0,20,0,20,0,1);
epix=size(ffb,1);
ffs=(ffs-mean(ffs(:)))./std(ffs(:));
ffs=ffs*21+40;
ffb=(ffb-mean(ffb(:)))./std(ffb(:));
ffb=ffb*21+40;

userange=56*0+(1:8:56)+1;
userange2=56*1+(1:8:56)+1;
ffg1=reshape(ffb(:,:,userange),epix*epix,7);
ffg2=reshape(ffb(:,:,userange2),epix*epix,7);
pffs=movpower((ffs-mean(ffs(:))./std(ffs(:)))*21+40,0,0,1,0.5,0);
ppg1=pffs(:,userange)-repmat(mean(pffs(:,1:end),2),[1 length(userange)]);
ppg2=pffs(:,userange2)-repmat(mean(pffs(:,1:end),2),[1 length(userange)]);

showkern(cat(3,nan.*ones([size(ffg1) 6]),ffg1,nan.*ones([size(ffg1)]),...
             nan.*ones([size(ffg1)])),'space');

tppg=sf2gr(ppg1,15,10,0,0,'pfft');
mm=max(abs(tppg(:)));

for ii=1:size(ppg1,2),
   ta=flipud(tppg(:,:,ii)');
   ta=repmat(ta,[1 1 3])./mm;
   ta(:,:,1)=1+ta(:,:,1).*(ta(:,:,1)<0);
   ta(:,:,2)=1-abs(ta(:,:,2));
   ta(:,:,3)=1-ta(:,:,3).*(ta(:,:,3)>0);
   hs=subplot(9,size(ppg1,2),7*size(ppg1,2)+ii);
   imagesc(ta);
   axis image;
   set(hs,'XTickLabel',[],'YTickLabel',[]);
end

if 0
   tppg=sf2gr(ppg2,15,10,0,0,'pfft');
   mm=max(abs(tppg(:)));
   
   for ii=1:size(tppg,3),
      ta=flipud(tppg(:,:,ii)');
      ta=repmat(ta,[1 1 3])./mm;
      ta(:,:,1)=1+ta(:,:,1).*(ta(:,:,1)<0);
      ta(:,:,2)=1-abs(ta(:,:,2));
      ta(:,:,3)=1-ta(:,:,3).*(ta(:,:,3)>0);
      hs=subplot(9,size(ppg1,2),8*size(ppg1,2)+ii);
      imagesc(ta);
      axis image;
      set(hs,'XTickLabel',[],'YTickLabel',[]);
   end
end

set(gcf,'PaperOrientation','portrait','PaperPosition',[2.25 4 4 4]);




sharp=median(z2.twopeakscore(find(edparms2(:,2)<6)));
smooth=median(z2.twopeakscore(find(edparms2(:,2)>6)));
esharp=std(z2.twopeakscore(find(edparms2(:,2)<6)))./sqrt(sum(edparms2(:,2)<6));
esmooth=std(z2.twopeakscore(find(edparms2(:,2)>6)))./sqrt(sum(edparms2(:,2)>6));
return

gidx=find(cc2good);
mm=zeros(length(gidx),1);
ee=zeros(length(gidx),1);
for ii=1:length(gidx),
   mm(ii)=mean(z2.sftuning(1,gidx(ii),:));
   ee(ii)=std(z2.sftuning(1,gidx(ii),:)).*sqrt(size(z2.sftuning,3)- ...
                                               1);
   if mm(ii)+ee(ii).*2 <= 0,
      fprintf('%s\n',z2.cellids{gidx(ii)});
   end
end


keyboard

nac=zeros(max(edparms2(:,2)),2);
nac(:,1)=hist(edparms1(:,2),1:max(edparms1(:,2)))';
nac(:,2)=hist(edparms2(:,2),1:max(edparms2(:,2)))';

nac2=nac(1:4,:);
nac2(1:3,:)=nac2(1:3,:)+nac(7:-1:5,:);
nac2(:,1)=nac2(:,1)./sum(nac2(:,1));
nac2(:,2)=nac2(:,2)./sum(nac2(:,2));
plot(nac2);


figure(5);
clf;

% orbw vs ncart enhancment (over cart)
oridx=find(dg2);

subplot(2,2,1);
plot(z2.mcstd(oridx),z2dnorm(oridx,2)-z2dmax(oridx,1),'.')
title(sprintf('orbw v nc-cart r=%.2f',...
              xcov(z2.mcstd(oridx),z2dnorm(oridx,2)-z2dmax(oridx,1),0,'coeff')));

subplot(2,2,2);
plot(z2.mcstd(oridx),z2dnorm(oridx,3)-z2dmax(oridx,1),'.');
title(sprintf('orbw v nat-cart r=%.2f',...
              xcov(z2.mcstd(oridx),z2dnorm(oridx,3)-z2dmax(oridx,1),0,'coeff')));

sfidx=find(sfcat2==1);
sfidx=find(dg2);

subplot(2,2,3);
plot(sfbw2(sfidx),z2dnorm(sfidx,2)-z2dmax(sfidx,1),'.')
title(sprintf('sfbw v nc-cart r=%.2f',...
              xcov(sfbw2(sfidx),z2dnorm(sfidx,3)-z2dmax(sfidx,1),0,'coeff')));

subplot(2,2,4);
plot(sfbw2(sfidx),z2dnorm(sfidx,3)-z2dmax(sfidx,2),'.')
title(sprintf('sfbw v nat-nc r=%.2f',...
              xcov(sfbw2(sfidx),z2dnorm(sfidx,3)-z2dmax(sfidx,2),0,'coeff')));



% 2pk vs hyp enhancment (over pol)
pp=[cat(1,z2.res.cartmax) cat(1,z2.res.polmax) ...
    cat(1,z2.res.hypmax)  cat(1,z2.res.natmax)];
plot(z2.twopeakscore,pp(:,6)-pp(:,4),'.')




