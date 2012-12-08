% v1asympsum.m

% callmaster : exp X in-model X cell
% predxc: exp X in-model X cell X out-model X cnf

dbopen;

figure(1);
figure(2);

SPACE=0;
OVERLAPONLY=1;

if SPACE,
   batches=[52; ...   % rev: 1.0 0.95 0.9 0.85
            56; ...   % gr:  ""
            60];      % nr:  ""
elseif 1,
   batches=[51; ...   % rev: 1.0 0.95 0.9 0.85
            55; ...   % gr:  ""
            59];      % nr:  ""
end

yn=input('reload and re-calc fredfrac (y/[n])? ','s');
fname=sprintf('/auto/k5/david/tmp/v1asympsum.s%d.o%d.mat',...
              SPACE,OVERLAPONLY);

if length(yn)>0 & yn(1)=='y',
   
   % load pred results for each batch
   res=loadbigdata(batches);
   predxc=permute(res.predxc,[1 2 6 4 5 3]);
   prederr=permute(res.prederr,[1 2 6 4 5 3]);
   predp=permute(res.predp,[1 2 6 4 5 3]);
   predmse=permute(res.predmse,[1 2 6 4 5 3]);
   celllist=res.celllist;
   totlen=res.totlen;
   spikes=res.spikes;
   expxc=permute(res.expxc,[1 2 5 4  3]);
   
   badcellidx=find(strcmp(celllist,'R150B') | strcmp(celllist,'93G83A'));
   predxc(badcellidx,1,:,:,:)=nan;
   
   % skip cells that don't have data from all batches
   if OVERLAPONLY,
      nonoveridx=find(sum(~isnan(predxc(:,:,1,1,1)),2) < size(batches,1));
      predxc(nonoveridx,:)=nan;
      expxc(nonoveridx,:)=0;
   end
   
   cellrange=find(sum(~isnan(predxc(:,:)),2) > 0 )';
   
   ccr=ones(length(celllist),length(batches)).*nan;
   ccrM=ones(length(celllist),length(batches)).*nan;
   MtoInf=ones(length(celllist),length(batches)).*nan;
   OnetoInf=ones(length(celllist),length(batches)).*nan;
   Mval=zeros(length(celllist),length(batches));
   predr2err=zeros(size(predxc)).*nan;
   for cellidx=cellrange,
      batchidx=1;
      [tccr,tccrM,tM,tvpredxc,txcvstrials]=...
          xcfredfrac(celllist{cellidx},batches(batchidx));
      if length(tccr)>0,
         ccr(cellidx,:)=tccr;
         ccrM(cellidx,:)=tccrM;
         Mval(cellidx,:)=tM;
         tccr=tccr+(tccr==0);
         MtoInf(cellidx,:)=tccrM./tccr;
         OnetoInf(cellidx,:)=1./tccr;
         
         ppe=squeeze(std(tvpredxc.*abs(tvpredxc),0,3).*sqrt(size(tvpredxc,3)));
         predr2err(cellidx,batchidx,:,:,:)=permute(ppe,[3 1 2]);
         
         for predidx=1:length(txcvstrials),
            xcvstrials(cellidx,batchidx,predidx)=txcvstrials(predidx);
         end
      end
      
      for batchidx=2:length(batches),
         if ~isnan(predxc(cellidx,batchidx,1,batchidx,1)),
            sql=['SELECT * FROM sRunData WHERE cellid="',celllist{cellidx},'"',...
                 ' AND batch=',num2str(batches(batchidx))];
            rundata=mysql(sql);
            
            z=zload([rundata.respath,rundata.resfile,'.gz']);
            tvpredxc=z.vpredxc;
            ppe=squeeze(std(tvpredxc.*abs(tvpredxc),0,3).*sqrt(size(tvpredxc,3)));
            predr2err(cellidx,batchidx,:,:,:)=permute(ppe,[3 1 2]);
         end
      end
   end

   % alternative, non-fred way of estimating validation noise ceiling
   MtoInf2=zeros(size(MtoInf)).*nan;
   OnetoInf2=zeros(size(OnetoInf)).*nan;
   fitopts=optimset('Display','off');
   nlidx=2;

   for cellidx=1:size(xcvstrials,1),
      %fprintf('cellid %d (%s)\n',cellidx,celllist{cellidx});
      batchidx=1;
      %for batchidx=1:length(batches),
      if sum(isnan(xcvstrials(cellidx,1,batchidx).n))==0 & ...
             length(xcvstrials(cellidx,1,batchidx).n)>0,
         
         n=xcvstrials(cellidx,1,batchidx).n;
         cxy=xcvstrials(cellidx,1,batchidx).cxy;
         
         bootcount=size(cxy,3);
         vfit=zeros(bootcount,2);
         for bootidx=1:bootcount,
            m=mean(cxy(:,:,bootidx),2);
            s=std(cxy(:,:,bootidx),0,2);
            if sum(abs(m))>0,
               vfit(bootidx,:)=polyfit(1./n,1./m,1);
               vfit(bootidx,:)=fminsearch('invlineerr',vfit(bootidx,:),...
                                          fitopts,n,m,s);
            end
         end
         
         vfitgood=find(vfit(:,2)>0);
         if length(vfitgood)>0,
            onexc=squeeze(mean(cxy(1,:,vfitgood),2));
            Mxc=squeeze(mean(cxy(end,:,vfitgood),2));
            infxc=1./vfit(vfitgood,2);
            
            OnetoInf2(cellidx,batchidx)=mean(infxc./onexc);
            MtoInf2(cellidx,batchidx)=mean(infxc./Mxc);
         end
         
         xcobs=mean(cxy(end,:));
         fprintf('%d cellid %s: %.3f  %.3f --> %.3f v %.3f\n',...
                 cellidx,celllist{cellidx},...
                 predxc(cellidx,batchidx,1,batchidx,nlidx) .* ...
                 abs(predxc(cellidx,batchidx,1,batchidx,nlidx)),...
                 xcobs,xcobs.*MtoInf(cellidx,batchidx),...
                 xcobs.*MtoInf2(cellidx,batchidx));
      end
   end
   
   save(fname);
else
   load(fname);
end
   
expstr={'Nat','DGrat','DNat'};
instr={'pix','pow','psf'};
cnfstr=expstr;
outstr={'lin','rec','full'};
datastr={'cnfxc2','expxc','cnfmse','cnfinf2'};

% if 1, take all cells
%PMIN=2.0;
PMIN=0.2;
%PMIN=0.05;

frac=[1.0 0.85 0.7 0.55 0.4 0.25 0.1];
nfit=length(frac)-1;
DOMEDIAN=0;

fprintf('SPACE=%d PMIN=%.2f nfit=%d\n',SPACE,PMIN,nfit);

exprange= [1 1 1 1 1 1 ; 3 3 3 3 3 3 ; 2 2 2 2 2 2 ];
cnfrange= [1 1 3 3 2 2 ; 1 1 3 3 2 2 ; 1 1 3 3 2 2 ];
nlrange=  [2 2 2 2 2 2 ; 2 2 2 2 2 2 ; 2 2 2 2 2 2 ];
datarange=[1 4 1 4 1 4 ; 1 4 1 4 1 4 ; 1 4 1 4 1 4 ];
%datainv=  [0 0 0 0 0 0 ; 0 0 0 0 0 0 ; 0 0 0 0 0 0 ];
datainv=  [1 1 1 1 1 1 ; 1 1 1 1 1 1 ; 1 1 1 1 1 1 ];

%
% calculate asymptotic confirmatory preds
%

rowcount=size(exprange,2);
colcount=size(exprange,1);

meanxc=zeros(size(exprange));
asympxc=zeros(size(exprange));
bfitgood=zeros([size(predxc,1) size(exprange)]);

for ii=1:length(exprange(:)),
   %subplot(rowcount,colcount,ii);
   
   if datarange(ii)==1,
      % confirmatory xc^2
      prset=squeeze(predxc(:,exprange(ii),:,cnfrange(ii),nlrange(ii)));
      prset=prset.*abs(prset); % preserve sign in r^2
      
      prerr=squeeze(predr2err(:,exprange(ii),:,cnfrange(ii),nlrange(ii)));
      goodidx=find(~isnan(prset(:,1)));
   elseif datarange(ii)==2,
      % exploratory xc^2
      prset=squeeze(expxc(:,exprange(ii),:,nlrange(ii)).^2);
      goodidx=find(prset(:,1)~=0);
   elseif datarange(ii)==3,
      % mse frac of total variance
      prset=squeeze(prederr(:,exprange(ii),:,cnfrange(ii),nlrange(ii)));
      goodidx=find(~isnan(prset(:,1)));
   elseif datarange(ii)==4,
      % cnf xc scaled to compensate for finite cnf trials
      if cnfrange(ii)==1,
         mtoi=MtoInf2(:,1);
      else
         mtoi=OnetoInf2(:,1);
      end
      over10idx=find(Mval(:,1)>=10);
      mtoi(find(Mval(:,1)<10))=mean(OnetoInf2(over10idx,1));
      
      onetoi=OnetoInf2(:,1);
      onetoi(find(Mval(:,1)<10))=mean(OnetoInf2(over10idx,1));
      
      prset=squeeze(predxc(:,exprange(ii),:,cnfrange(ii),nlrange(ii)));
      prset=prset .* abs(prset) .* repmat(mtoi,1,size(predxc,3));
      prerr=squeeze(predr2err(:,exprange(ii),:,cnfrange(ii),nlrange(ii))) ...
            .* repmat(mtoi,1,size(predxc,3));
      goodidx=find(~isnan(prset(:,1)));
      
   end
   
   prset=prset(goodidx,:);
   prerr=prerr(goodidx,:);
   lens=totlen(goodidx,exprange(ii));
   ll=lens*frac(1:nfit);
   
   % find asymptote for each cell.
   pp=zeros(length(goodidx),2);
   ppe=zeros(length(goodidx),1);
   pprset=zeros(size(prset));
   ppasymp=zeros(length(goodidx),1);
   for cc=1:length(prset),
      f=1./(frac.*lens(cc));  % ie, x-axis is 1/T
      
      if datainv(ii) & ~sum(prset(cc,:)==0),
         pp(cc,:)=polyfit(f(1:nfit),1./prset(cc,1:nfit),1);
         estart=invlineerr(pp(cc,:),frac(1:nfit).*lens(cc),...
                            prset(cc,1:nfit),prerr(cc,1:nfit));
         
         %fitopts=optimset('Display','off','LargeScale','off');
         fitopts=optimset('Display','off');
         ppbetter=fminsearch('invlineerr',pp(cc,:),...
                             fitopts,frac(1:nfit).*lens(cc),...
                             prset(cc,1:nfit),prerr(cc,1:nfit));
         estop=invlineerr(ppbetter,frac(1:nfit).*lens(cc),...
                            prset(cc,1:nfit),prerr(cc,1:nfit));
         %fprintf('invlineer: estart=%.3f estop=%.3f\n',...
         %        estart,estop);
         
         pp(cc,:)=ppbetter;
         pprset(cc,:)=1./(pp(cc,2)+pp(cc,1).*f);
         ppasymp(cc)=1./pp(cc,2);
      else
         pp(cc,:)=polyfit(f(1:nfit),prset(cc,1:nfit),1);
         pprset(cc,:)=(pp(cc,2)+pp(cc,1).*f);
         ppasymp(cc)=pp(cc,2);
      end
      
      
      % error of linear fit, to determine if this cell qualifies
      if sum(prset(cc,:)==0),
         ppe(cc)=1;
      else
         ppe(cc)=invlineerr(pp(cc,:),frac(1:nfit).*lens(cc),...
                            prset(cc,1:nfit),prerr(cc,1:nfit));
      end
      
      
      
      %bfitgood(goodidx(cc),ii)=(ppe(cc)<PMIN);
      %bfitgood(goodidx(cc),ii)=(ppe(cc)<PMIN & ppasymp(cc)<2 & lens(cc)>600);
      bfitgood(goodidx(cc),ii)=(ppe(cc)<PMIN & ...
                                ppasymp(cc)<1.5 & ppasymp(cc)>-0.01);
      
      if 0 & ii==1,
         figure(2);
         subplot(ceil(length(goodidx)./10),10,cc);
         plot(1./f,prset(cc,:));
         hold on
         plot(1./f,pprset(cc,:),'r--');
         plot(1./f,ppasymp(cc).*ones(size(f)),'k--');
         hold off
         ht=title(sprintf('%s %.2f',celllist{goodidx(cc)},ppe(cc)));
         if bfitgood(goodidx(cc),ii),
            set(ht,'FontWeight','bold');
         end
      end
   end
   
   fitidx=find(bfitgood(goodidx,ii));
   
   if datainv(ii),
      % find asymptotic pred for 1./mean across all cells.
      
      f=1./(frac.*mean(lens(fitidx)));
      tpr=prset(fitidx,:);
      mtpr=mean(prset(fitidx,:),1);
      ppmean=mean(1./pp(fitidx,:));
      ppmedian=median(1./pp(fitidx,:));
      
      amax=ppmean(2);
            
   elseif datainv(ii),
      % find asymptotic pred for 1./mean across all cells.
      
      f=1./(frac.*mean(lens(fitidx)));
      tpr=prset(fitidx,:);
      mtpr=mean(prset(fitidx,:),1);
      ppavg=polyfit(f(1:nfit),1./mean(tpr(:,1:nfit),1),1);
      ppmean=1./ppavg;
      
      if ppmean(2)<max(mtpr),
         
         ppmean(2)=max(mtpr);
      end
      
      pm=polyfit(f(1:nfit),1./median(prset(fitidx,1:nfit),1),1);
      ppmed=1./pm;
      
      amax=ppmean(2);
            
   elseif ~datainv(ii),
      % find asymptotic pred for mean across all cells.
      
      tpr=prset(fitidx,:);
      mtpr=mean(prset(fitidx,:),1);         
      f=1./(frac.*mean(lens(fitidx)));
      pp=polyfit(f(1:nfit),mean(tpr(:,1:nfit),1),1);
      ppmean=pp;
      
      if ppmean(2)<max(mtpr),
         ppmean(2)=max(mtpr);
      end
      
      pm=polyfit(f(1:nfit),median(prset(fitidx,1:nfit),1),1);
      ppmed=pm;
      
      amax=ppmean(2);
      
      cla
      fill([0 1.1 1.1 0],[0 0 1 1],[0.6 0.6 0.6],'LineStyle','none');
      hold on
      xx=[0     1.1  1.1           frac       0               0    ];
      yy=[amax  amax (mtpr(1)) (mtpr) (mtpr(end)) amax ];
      fill(xx,yy,[0.75 0.75 0.75],'LineStyle','none');
      fill([0 1.1 1.1 0],[amax amax 1 1],[0.9 0.9 0.9],'LineStyle','none');
      
      plot(frac,(pp(2)+pp(1).* f),'k--');
      plot([0 1.1],[amax amax],'k--');      
      ht=errorbar(frac,(mtpr),nanstd(tpr)./sqrt(size(tpr,1)),'k-');
      set(ht,'LineWidth',2);
      hold off
      a=axis;
      axis([0 1.1 0 1]);
      
      if datarange(ii)==4,
         subplot(rowcount,colcount,ii-colcount);
         hold on
         plot([0 1.1],[amax amax],'k--');
         fill([0 1.1 1.1 0],[amax amax 1 1],[1 1 1],'LineStyle','none');
         hold off
         subplot(rowcount,colcount,ii);
      end
      
   elseif 0,
      % fit single line to all data. screwey
      
      ppmean=mean(pp(fitidx,:),1);
      ppmed=median(pp(fitidx,:),1);
      
      f0=1./(lens(fitidx)*frac);
      p0=prset(fitidx,:);
      
      scatter(f0(:),p0(:),'bx');
      hold on
      plot([max(f0(:)) 0],ppmean(2)+ppmean(1).*[max(f0(:)) 0],'r--');
      hold off
   elseif 1,
      
      f0=sqrt(1./(lens(fitidx)*frac));
      p0=prset(fitidx,1:nfit);
      
      pp=polyfit(f0,p0,1);
      ppmean=pp;
      
      pm=polyfit(f0,p0,1);
      ppmed=pm;
      
      scatter(f0(:),p0(:),'bx');
      hold on
      plot([max(f0(:)) 0],ppmean(2)+ppmean(1).*[max(f0(:)) 0],'r--');
      hold off
      
   end
   if DOMEDIAN,
      asympxc(ii)=(ppmed(2));
      meanxc(ii)=(median(prset(fitidx,1)));
   else
      asympxc(ii)=(ppmean(2));
      meanxc(ii)=(mean(prset(fitidx,1)));
   end
   percentinc=(asympxc(ii)-meanxc(ii))./meanxc(ii);
   
   fprintf(['ii=%2d %5s/%5s %6s (n=%2d/%2d) xc: %.2f-->%.2f',...
            ' extrap: %.2f-->%.2f\n'],...
           ii,expstr{exprange(ii)},expstr{cnfrange(ii)},...
           datastr{datarange(ii)},...
           length(fitidx),length(goodidx),meanxc(ii),asympxc(ii),...
           (mean(prset(:,1))),((1+percentinc).*mean(prset(:,1))));
   
   valres(ii).asympxc=asympxc(ii);
   valres(ii).meanxc=meanxc(ii);
   valres(ii).percentinc=(asympxc(ii)-meanxc(ii))./meanxc(ii);
   valres(ii).expstr=expstr{exprange(ii)};
   valres(ii).cnfstr=expstr{cnfrange(ii)};
   valres(ii).datastr=datastr{datarange(ii)};
   valres(ii).goodidx=goodidx;
   valres(ii).fitidx=fitidx;
   valres(ii).prset=prset;
   valres(ii).lens=lens;
   valres(ii).pp=pp;
   valres(ii).ppmean=ppmean;
   valres(ii).ppmedian=ppmedian;
   
end

% population summary noise for within class (figure 6)

figure(1);
clf
if OVERLAPONLY
   showdat=[1 2 3 7 8 9 13 14 15];
else
   showdat=[1 8 15];
end

for showdatidx=1:length(showdat),
   ii=showdat(showdatidx);
   
   subplot(ceil(length(showdat)/3),3,showdatidx);
   
   % average over all cells (good and bad exp noise model fits)
   nall=size(valres(ii).prset,1);
   meanall=mean(valres(ii).prset(:,1));
   ppmeanall=mean(1./valres(ii).pp(:,2));
   ppmeanall2=mean(1./valres(ii+3).pp(:,2));
   
   % average over all good cells (good exp noise model fit)
   fitidx=(1:nall)';
   %fitidx=valres(ii).fitidx;
   n=length(fitidx);
   prset=valres(ii).prset(fitidx,:);
   prset2=valres(ii+3).prset(fitidx,:);
   obsxc=mean(prset);
   obsxc2=mean(prset2);
   
   if showdatidx==1,
      
      mtoi=OnetoInf2(:,1);
      over10idx=find(Mval(:,1)>=10 & ~isnan(mtoi));
      under10idx=find(Mval(:,1)<10 | isnan(mtoi));
      mtoi(under10idx)=mean(mtoi(over10idx,1));
      
      a=repmat(mtoi(valres(ii).goodidx),1,size(prset2,2));
      a(find(isnan(a)))=1;
      prset1=prset2./a;
   else
      prset1=prset;
   end
   
   obsxcboot=zeros(size(prset,1),length(obsxc));
   pfitboot=zeros(size(prset,1),2);
   obsxc1boot=zeros(size(prset,1),length(obsxc));
   pfit1boot=zeros(size(prset,1),2);
   obsxc2boot=zeros(size(prset,1),length(obsxc));
   pfit2boot=zeros(size(prset,1),2);
   bootcount=size(prset,1);
   for jj=1:bootcount,
      useidx=[1:jj-1 jj+1:bootcount];
      obsxcboot(jj,:)=mean(prset(useidx,:));
      
      % fit exp noise model to pred data
      a=1./frac;
      b=1./obsxcboot(jj,:);
      pfitboot(jj,:)=polyfit(a(1:nfit),b(1:nfit),1);
      
      obsxc1boot(jj,:)=mean(prset1(useidx,:));
      b=1./obsxc1boot(jj,:);
      pfit1boot(jj,:)=polyfit(a(1:nfit),b(1:nfit),1);
      
      obsxc2boot(jj,:)=mean(prset2(useidx,:));
      b=1./obsxc2boot(jj,:);
      pfit2boot(jj,:)=polyfit(a(1:nfit),b(1:nfit),1);
   end
   obsxcerr=std(obsxcboot).*sqrt(bootcount); 
   obsxc1err=std(obsxc2boot).*sqrt(size(prset1,1)); 
   obsxc2err=std(obsxc2boot).*sqrt(size(prset2,1)); 
   
   pfitmean=mean(1./pfitboot(:,2));
   pfiterr=std(1./pfitboot(:,2)).*sqrt(bootcount);
   pfit1mean=mean(1./pfit1boot(:,2));
   pfit1err=std(1./pfit1boot(:,2)).*sqrt(bootcount);
   pfit2mean=mean(1./pfit2boot(:,2));
   pfit2err=std(1./pfit2boot(:,2)).*sqrt(bootcount);
   
   goodtoallpercent=ppmeanall/pfitmean;
   tovalcorrpercent=pfit2mean/pfitmean;
   
   bootxx=frac;
   bootyy=1./(pfitboot(:,1)*(1./bootxx) + ...
              repmat(pfitboot(:,2),1,length(bootxx)));
   yyerr=std(bootyy).*sqrt(bootcount);
   yymean=mean(bootyy);
   
   bootyy1=1./(pfit1boot(:,1)*(1./bootxx) + ...
              repmat(pfit1boot(:,2),1,length(bootxx)));
   yy1err=std(bootyy1).*sqrt(bootcount);
   yy1mean=mean(bootyy1);
   
   bootyy2=1./(pfit2boot(:,1)*(1./bootxx) + ...
              repmat(pfit2boot(:,2),1,length(bootxx)));
   yy2err=std(bootyy2).*sqrt(bootcount);
   yy2mean=mean(bootyy2);
   
   fill([0 1.5 1.5 0],[0 0 1 1],[0.6 0.6 0.6],'LineStyle','none');
   hold on
   xx=[0     1.5  1.5           frac       0               0    ];
   yy=[pfitmean  pfitmean (yymean(1)) (yymean) (yymean(end)) pfitmean ];
   fill(xx,yy,[0.75 0.75 0.75],'LineStyle','none');
   fill([0 1.5 1.5 0],[pfitmean pfitmean 1 1],[0.9 0.9 0.9],...
        'LineStyle','none');
   
   %plot([0 1.5],[pfit2mean pfit2mean],'k--');
   fill([0 1.5 1.5 0],[pfit2mean pfit2mean 1 1],[1 1 1],...
        'LineStyle','none');
   
   ht=errorbar(bootxx,yy1mean,yy1err,'k:');
   ht=errorbar(max(frac)*1.4,pfit1mean,pfit1err,'ko');
   
   ht=errorbar(bootxx,yymean,yyerr,'k-');
   ht=errorbar(max(frac)*1.4,pfitmean,pfiterr,'ko');
   
   ht=errorbar(bootxx,yy2mean,yy2err,'k--');
   ht=errorbar(max(frac)*1.4,pfit2mean,pfit2err,'ko');
   
   hold off
   
   a=axis;
   axis([0 1.5 0 1]);
   
   fprintf('%d %s/%s %d: %.2f-->%.2f-->%.2f %d: %.2f-->%.2f-->%.2f\n',...
           ii,valres(ii).expstr,valres(ii).cnfstr,...
           n,obsxc(1),pfitmean,pfit2mean,...
           nall,meanall,ppmeanall,ppmeanall2);
   title(sprintf('%s/%s %d: %.2f-->%.2f-->%.2f',...
                 valres(ii).expstr,valres(ii).cnfstr,...
                 n,obsxc(1),pfitmean,pfit2mean));
   
end

xlabel(sprintf('space=%d olap=%d',SPACE,OVERLAPONLY));
if OVERLAPONLY
   set(gcf,'PaperPosition',[0.25 2.5 8 6],'PaperOrientation','portrait')
else
   set(gcf,'PaperPosition',[0.25 4.5 8 2],'PaperOrientation','portrait')
end

disp('skipping demo cell');
return


% demo cells (figure 5)
dbopen;
figure(2);
clf
%crr=[ 66 66 66 67 67 67 75 75 75];
%brr=[  1  3  2  1  3  2  1  3  2 ];
crr=[  67 ];
brr=[   1 ];
nfitdemo=length(frac);
for demoidx=1:length(crr),
   
   cellidx=crr(demoidx);
   batchidx=brr(demoidx);
   sql=['SELECT * FROM sRunData WHERE cellid="',celllist{cellidx},'"',...
        ' AND batch=',num2str(batches(batchidx))];
   rundata=mysql(sql);
   z=zload([rundata.respath,rundata.resfile,'.gz']);
   
   nlidx=2;
   bootcount=z.params.fitboot;
   pfit=zeros(bootcount,2);
   for bootidx=1:bootcount,
      useidx=[1:bootidx-1 bootidx+1:bootcount];
      pp=squeeze(mean(z.vpredxc(batchidx,nlidx,useidx,:) .* ...
                      abs(z.vpredxc(batchidx,nlidx,useidx,:)),3));
      %pp=squeeze(mean(z.vpredxc(batchidx,nlidx,bootidx,:) .* ...
      %                abs(z.vpredxc(batchidx,nlidx,bootidx,:)),3));
      ppe=squeeze(std(z.vpredxc(batchidx,nlidx,useidx,:) .* ...
                      abs(z.vpredxc(batchidx,nlidx,useidx,:)),0,3)) .* ...
          sqrt((z.params.fitboot-1).*frac(:));
      if batchidx==1,
         ppInf=MtoInf2(cellidx,1).*pp;
      else
         ppInf=OnetoInf2(cellidx,1).*pp;
      end
      
      % fit exp noise model to pred data
      a=1./(frac.*res.totlen(cellidx,1));
      b=1./pp;
      
      pfit(bootidx,:)=polyfit(a(:),b(:),1);
      tpfit=fminsearch('invlineerr',pfit(bootidx,:),fitopts,...
                      frac(1:nfitdemo).*res.totlen(cellidx,1),...
                       pp(1:nfitdemo)',ppe(1:nfitdemo)');
      pfit(bootidx,:)=tpfit(:)';
      
      %[tpfit,r,j]=nlinfit(frac(1:nfitdemo).*res.totlen(cellidx,1),...
      %                       pp(1:nfitdemo)','invline',pfit(bootidx,:));
      %ci = nlparci(pfit(bootidx,:),r,j);
   end
   
   pfiterr=std(1./pfit(:,2)).*sqrt(bootcount);
   pfitmean=mean(pfit,1);
   
   ff=repmat(frac.*res.totlen(cellidx,1),[z.params.fitboot 1]);
   ppall=squeeze(z.vpredxc(batchidx,nlidx,:,:) .* ...
                 abs(z.vpredxc(batchidx,nlidx,:,:)));
   
   % fit exp noise model to val noise-corrected pred data
   pp=squeeze(mean(ppall));
   
   ppe=squeeze(std(ppall .* abs(ppall))) .* sqrt((z.params.fitboot));
   if batchidx==1,
      ppInf=MtoInf2(cellidx,1).*pp;
   else
      ppInf=OnetoInf2(cellidx,1).*pp;
   end
   b2=1./ppInf;
   pfit2=polyfit(a(:),b2(:),1);
   pfit2=fminsearch('invlineerr',pfit2,fitopts,...
                    frac(1:nfitdemo).*res.totlen(cellidx,1),...
                    ppInf(1:nfitdemo),ppe(1:nfitdemo));
   
   ff=repmat(frac.*res.totlen(cellidx,1),[z.params.fitboot 1]);
   ppall=squeeze(z.vpredxc(batchidx,nlidx,:,:) .* ...
                 abs(z.vpredxc(batchidx,nlidx,:,:)));
   
   if length(crr)<3,
      valres=xctestpreds(rundata.id,rundata.batch,nlidx,brr(demoidx));
      subplot(length(crr),2,demoidx*2);
      a=repmat(valres.n(:,1),1,size(valres.cxy,2));
      b=mean(valres.cxy,3);
      
      bootcount=size(b,2);
      vfit=zeros(2,bootcount);
      vpred=zeros(size(b));
      for bootidx=1:bootcount,
         n=valres.n(:,1);
         m=mean(b(1:end-1,[1:bootidx-1 bootidx+1:size(b,2)]),2);
         vfit(:,bootidx)=polyfit(1./n(1:end-1),1./m,1)';
         vpred(:,bootidx)=1./(vfit(1,bootidx)*(1./n) + vfit(2,bootidx));
      end
      
      scatter(a(:),b(:),2,'k','filled');
      hold on
      errorshade(n,mean(vpred,2),std(vpred,0,2).*sqrt(bootcount),...
                 'k',[0.8 0.8 0.8]);
      %plot(valres.n(:,1),mean(valres.m,2),'k-');
      
      xinf=max(valres.n(:,1))+4;
      
      errorbar(xinf,mean(1./vfit(2,:)),...
               std(1./vfit(2,:)).*sqrt(bootcount),'k+');
      
      hold off
      axis([0 max(valres.n(:,1))+5 0 0.65]);
      
      title(sprintf('val: 1: %.2f M: %.2f Inf %.2f',...
                    mean(valres.m(1,:)),mean(valres.m(end,:)),...
                    ppInf(1)));
      
      subplot(length(crr),2,demoidx*2-1);
      
   else
      subplot(round(length(crr)./3),3,demoidx);
   end   
   
   xx=[max(frac).*[2.0 1.8 1.5 1.2] frac].*res.totlen(cellidx,1);
   
   bootxx=xx; % frac.*res.totlen(cellidx,1);
   bootyy=1./(pfit(:,1)*(1./bootxx) + repmat(pfit(:,2),1,length(bootxx)));
   yyerr=std(bootyy).*sqrt(bootcount);
   yymean=mean(bootyy);
   
   %errorbar(bootxx./1000,yymean,yyerr,'k-');
   errorshade(bootxx./1000,yymean,yyerr,'k',[0.8 0.8 0.8]);
   hold on
   
   %plot(frac.*res.totlen(cellidx,1)./1000,pp,'k+');
   
   xinf=max(frac).*res.totlen(cellidx,1).*2.4./1000;
   %scatter(xinf,1./pfit(2),'ko','filled');
   errorbar(xinf,1./pfitmean(2),pfiterr,'k+');
   
   %scatter(xinf,1./pfit2(2),'ks','filled');
   scatter(ff(:)./1000,ppall(:),2,'k','filled');
   
   %errorbar(frac.*res.totlen(cellidx,1),pp,ppe,'k-');
   %hold on
   %plot(frac.*res.totlen(cellidx,1),ppInf,'k--');
   %hold off
   %plot(xx./1000,1./(pfitmean(1)./xx + pfitmean(2)),'k-');
   %plot(xx./1000,1./(pfit2(1)./xx + pfit2(2)),'k--');
   
   %errorbar(frac.*res.totlen(cellidx,1)./1000,pp,ppe,'k+');
   
   hold off
   
   title(sprintf('cell %s (%d), r2=%.3f +e %.3f +v %.3f',...
                 celllist{cellidx},batches(batchidx),...
                 mean(ppall(:,1)),1./pfitmean(2),1./pfit2(2)));
   xlabel('frames (1000s)');
   ylabel('r^2');
   axis([0 xinf+1 0 0.65]);
   
   drawnow;
   
end

set(gcf,'PaperPosition',[1.25 4.5 6 2],'PaperOrientation','portrait')



figure(3);
clf
Mgoodidx=find(Mval(:,1)>=10);
scatter(OnetoInf(Mgoodidx,1),OnetoInf2(Mgoodidx,1),Mval(Mgoodidx));
%plot(Mval(Mgoodidx),((OnetoInf(Mgoodidx,1)-OnetoInf2(Mgoodidx,1))),'o');
axis([0.9 6 0.9 6]);
axis square



return


% figure(2);

subplot(2,3,4);
m1=diag(expmeanxc);
m2=diag(expasympxc);

b2=bar(m2);
hold on
set(b2(1),'CDataMapping','direct','FaceVertexCData',zeros(16,1)+64)
set(b2(2),'CDataMapping','direct','FaceVertexCData',zeros(16,1)+64)
set(b2(3),'CDataMapping','direct','FaceVertexCData',zeros(16,1)+64)
b1=bar(m1);
set(b1(1),'CDataMapping','direct','FaceVertexCData',zeros(16,1)+0)
set(b1(2),'CDataMapping','direct','FaceVertexCData',zeros(16,1)+25)
set(b1(3),'CDataMapping','direct','FaceVertexCData',zeros(16,1)+50)
hold off

axis([0 4 0 1.0]);
title(sprintf('asymp exp preds (n=%d)',length(goodidx)));
set(gca,'XTickLabel',{'Nat','DNat','DGrat'});
xlabel('validation class');
%legend(b1,'Nat','DNat','DGrat');

subplot(2,3,5);
m1=diag(diag(meanxc));
m2=diag(diag(asympxc));

b2=bar(m2);
hold on
set(b2(1),'CDataMapping','direct','FaceVertexCData',zeros(16,1)+64)
set(b2(2),'CDataMapping','direct','FaceVertexCData',zeros(16,1)+64)
set(b2(3),'CDataMapping','direct','FaceVertexCData',zeros(16,1)+64)
b1=bar(m1);
set(b1(1),'CDataMapping','direct','FaceVertexCData',zeros(16,1)+0)
set(b1(2),'CDataMapping','direct','FaceVertexCData',zeros(16,1)+25)
set(b1(3),'CDataMapping','direct','FaceVertexCData',zeros(16,1)+50)
hold off

axis([0 4 0 1.0]);
title(sprintf('asymp within preds (n=%d)',length(goodidx)));
set(gca,'XTickLabel',{'Nat','DNat','DGrat'});
xlabel('validation class');
%legend(b1,'Nat','DNat','DGrat');

subplot(2,3,6);
m1=meanxc;
m2=asympxc;

b2=bar(m2);
hold on
set(b2(1),'CDataMapping','direct','FaceVertexCData',zeros(16,1)+64)
set(b2(2),'CDataMapping','direct','FaceVertexCData',zeros(16,1)+64)
set(b2(3),'CDataMapping','direct','FaceVertexCData',zeros(16,1)+64)
b1=bar(m1);
set(b1(1),'CDataMapping','direct','FaceVertexCData',zeros(16,1)+0)
set(b1(2),'CDataMapping','direct','FaceVertexCData',zeros(16,1)+25)
set(b1(3),'CDataMapping','direct','FaceVertexCData',zeros(16,1)+50)
hold off

axis([0 4 0 1.0]);
title(sprintf('asymp cross preds (n=%d)',length(goodidx)));
set(gca,'XTickLabel',{'Nat','DNat','DGrat'});
xlabel('validation class');
%legend(b1,'Nat','DNat','DGrat');

colormap(gray)
set(gcf,'PaperPosition',[0.25 2.5 8 6],'PaperOrientation','portrait')

