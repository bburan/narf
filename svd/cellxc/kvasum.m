% function res=kvasum(batchid,dump)
%
function res=kvasum(batchid,dump)

dbopen;

if ~exist('batchid','var'),
   batchid=79
end
if ~exist('dump','var'),
   dump=0;
end
if ~exist('latidx','var'),
   latidx=1;
end

batchdata=dbget('sBatch',batchid);

sql=['SELECT sResults.*',...
     ' FROM sResults INNER JOIN sRunData',...
     ' ON sResults.runid=sRunData.id',...
     ' WHERE sResults.batch=',num2str(batchid),...
     ' AND not(sResults.matstr like "%m0000%")',...
     ' ORDER BY cellid'];
resdata=mysql(sql);

for ii=1:length(resdata),
   eval(char(resdata(ii).matstr));
end

okidx=zeros(length(preddata),1);
tt=0;
for ii=1:length(preddata),
   if ~isempty(preddata(ii).pxc),
      okidx(ii)=1;
      tt=tt+1;
      kerndata(tt)=preddata(ii);
   end
end
ttt=find(okidx)';

cellids={kerndata.cellid};
[sc,scidx]=sort(cellids);
%scidx=1:length(cellids);

fprintf('batchid=%d cellcount=%d\n',batchid,length(cellids));

if 0,
   keyboard
end

PFILT=0.01;
PATT=0.008; % 0.008;
PXCT=0.05;
% line from kernvsatt.m:

kerncount=length(kerndata);
predxc=zeros([size(kerndata(1).predxc),kerncount]).*nan;
if isfield(kerndata(1),'predinf'),
   predinf=zeros([size(kerndata(1).predinf),kerncount]).*nan;
end
if isfield(kerndata(1),'predxcpair'),
   predxcpair=zeros([size(kerndata(1).predxcpair),kerncount]).*nan;
   pxcpair=ones([size(kerndata(1).pxcpair),kerncount]);
end
%predxccross=zeros([size(kerndata(1).predxccross),kerncount]).*nan;
pxc=ones([size(kerndata(1).pxc),kerncount]);
%pxccross=ones([size(kerndata(1).pxccross),kerncount]);
randxc=cat(1,kerndata.randxc);

for kidx=1:kerncount,
   %s=size(kerndata(kidx).predxccross);
   s=size(kerndata(kidx).predxc);
   predxc(1:s(1),1:s(2),kidx)=kerndata(kidx).predxc;
   %predxccross(1:s(1),1:s(2),1:s(3),kidx)=kerndata(kidx).predxccross;
   pxc(1:s(2),1,kidx)=kerndata(kidx).pxc;
   %pxccross(1:s(2),1:s(3),kidx)=kerndata(kidx).pxccross;
   if isfield(kerndata(kidx),'predinf'),
      s=size(kerndata(kidx).predinf);
      predinf(1:s(1),1:s(2),kidx)=kerndata(kidx).predinf;
   end
   if isfield(kerndata(kidx),'predxcpair'),
      s=[size(kerndata(kidx).predxcpair) 1];
      
      if prod(s)>0,
         predxcpair(1:s(1),1:s(2),1:s(3),kidx)=kerndata(kidx).predxcpair;
         pxcpair(1:s(2),1:s(3),kidx)=kerndata(kidx).pxcpair;
      end
   end
end

nlcount=size(predxc,2);
xccount=pxc<PXCT;

if nlcount>=5,
   
   OUTNLMODE=4;
   if OUTNLMODE==0,
      nlnames={'no att','dc att','gain att','dc/g att',...
               'hng att','all att'};
   elseif OUTNLMODE==1,
      nlnames={'no att','x10','x90','D','A','all att'};
   elseif OUTNLMODE==2,
      nlnames={'no att','dc att','gain att','thr att',...
               'hng att','all att'};
   elseif OUTNLMODE==3,
      nlnames={'glob none','glob dc','glob dcg',...
               'dc att','gain att','all att'};
   elseif OUTNLMODE==4,
      if nlcount>5,
         nlnames={'glob none','glob dc','glob g','glob dcg',...
                  'allodd1','allodd2','allodd3','allodd4'};
         if nlcount>length(nlnames),
            nlcount=length(nlnames);
            predxc=predxc(:,1:nlcount,:,:);
            pxc=pxc(1:nlcount,:,:);
            PCXT2=1-(1-PXCT)^(1./(nlcount-4));
            xccount=[pxc(1:4,:,:)<PXCT; pxc(5:end,:,:)<PCXT2];
         end
      else
         nlnames={'visual','DC','gain','DC/g','tuning'};
      end
   end
   
   goodcheck=~isnan(squeeze(predxc(1,:,:))) & ...
             squeeze(predxc(1,:,:)) > repmat(randxc',[size(predxc,2) 1]);
elseif nlcount==2,
   % pairwise in-out comparison for fvvs
   
   OUTNLMODE=6;
   nlnames={'mean','diff'};
   
   predxc=squeeze(predxcpair(:,2,:,:));
   if size(predxc,2)>6,
      predxc=predxc(:,1:6,:);
   end
   nlcount=size(predxc,2);
   pxc=reshape(pxcpair(2,1:nlcount,:),nlcount,1,size(pxcpair,3));
   xccount=pxc<PXCT;
   
   for ii=1:nlcount,
      nlnames{ii}=sprintf('pair%d',ii);
   end
   
   goodcheck=~isnan(squeeze(predxc(1,:,:))) & ...
             squeeze(predxc(1,:,:)) > repmat(randxc',[size(predxc,2) 1]);
else
   % xcdms
   goodcheck=~isnan(squeeze(predxc(1,:,:)));
   nlnames={'ABab','Aa-Bb','AB-ab','A-B'};
   %nlnames={'AB','A-B','AB-ab','Aa-Bb'};
   OUTNLMODE=5;
end

fprintf('note: assuming OUTNLMODE=%d\n',OUTNLMODE);
if OUTNLMODE==4,
   
   fidx=find(goodcheck(5,:));
   fidx=1:size(goodcheck,2);
   fracxc=squeeze(predxc(1,[1 2 4 5],fidx))';
   fracxc(find(fracxc<0))=0;
   for ii=1:length(fracxc),
      fracmax=max(fracxc(ii,:));
      if fracmax==0,
         fracxc(ii,1)=1;
      else
         if fracxc(ii,4)==fracmax & xccount(5,1,fidx(ii)),
            fracxc(ii,4)=(fracxc(ii,4).^2-fracxc(ii,3).^2)./fracmax.^2;
         else
            fracxc(ii,4)=0;
         end
         if fracxc(ii,3)>fracxc(ii,2) & xccount(3,1,fidx(ii)),
            fracxc(ii,3)=(fracxc(ii,3).^2-fracxc(ii,2).^2)./fracmax.^2;
         else
            fracxc(ii,3)=0;
         end
         if fracxc(ii,2)>fracxc(ii,1) & xccount(2,1,fidx(ii)),
            fracxc(ii,2)=(fracxc(ii,2).^2-fracxc(ii,1).^2)./fracmax.^2;
         else
            fracxc(ii,2)=0;
         end
         if sum(fracxc(ii,2:end))>1,
            ii
         end
         fracxc(ii,1)=1-sum(fracxc(ii,2:end));
      end
   end
   
end

% header
fprintf('CELL:    ');
for nlidx=1:nlcount,
   ss=nlnames{nlidx}(1:min([6 length(nlnames{nlidx})]));
   fprintf(' %-6s ',ss);
end
fprintf('(in/out)\n');

% list pred results, on cell per line
cellid={};
fracxc=zeros(length(kerndata),4);
for ii=1:length(kerndata),
   kidx=scidx(ii);
   if OUTNLMODE==4,
      % rescale local atts to have a baseline pred equal to dcg
      % pred... compensates for extra noise in local strf fits
      if predxc(2,1,kidx)>0,
         predxc(1:2,2,kidx)=predxc(1:2,2,kidx)./...
             predxc(2,2,kidx).*predxc(1,1,kidx);
         predxc(1:2,3,kidx)=predxc(1:2,3,kidx)./...
             predxc(2,3,kidx).*predxc(1,1,kidx);
         predxc(1:2,4,kidx)=predxc(1:2,4,kidx)./...
             predxc(2,4,kidx).*predxc(1,1,kidx);
      end
      if predxc(2,end,kidx)>0,
         predxc(1:2,end,kidx)=predxc(1:2,end,kidx)./...
             predxc(2,end,kidx).*predxc(1,end-1,kidx);
      end
      
      if 1| goodcheck(5,kidx),
         tfracxc=predxc(1,[1 2 4 5],kidx);
         tfracxc(find(tfracxc<0))=0;
         tpct=zeros(size(tfracxc));
         
         fracmax=max(tfracxc);
         if fracmax==0,
            tfracxc(1)=1;
         else
            if tfracxc(2)>tfracxc(1) & xccount(2,1,kidx),
               tpct(2)=(tfracxc(2).^2-tfracxc(1).^2)./fracmax.^2;
            else
               tfracxc(2)=tfracxc(1);
            end
            if tfracxc(3)>tfracxc(2) & xccount(3,1,kidx),
               tpct(3)=(tfracxc(3).^2-tfracxc(2).^2)./fracmax.^2;
            else
               tfracxc(3)=tfracxc(2);
            end
            if tfracxc(4)==fracmax & xccount(5,1,kidx),
               tpct(4)=(tfracxc(4).^2-tfracxc(3).^2)./fracmax.^2;
            end
            if sum(tpct(2:end))>1.00001,
               ii
               keyboard
            end
            tpct(1)=1-sum(tpct(2:end));
         end
         fracxc(kidx,:)=tpct;
      else
         fracxc(kidx,:)=nan;
      end
   end
   
   fprintf('%-7s:',kerndata(kidx).cellid);
   cellid{ii}=kerndata(kidx).cellid;
   
   tpredxc=predxc(:,:,kidx);
   tpredxc(find(isnan(tpredxc)))=0;
   
   for nlidx=1:nlcount,
      fprintf(' .%0.2d/.%0.2d',...
              round(tpredxc(1,nlidx)*100*(tpredxc(1,nlidx)>0)),...
              round(tpredxc(2,nlidx)*100*(tpredxc(2,nlidx)>0)));
   end
   
   fprintf(' p<');
   %fprintf(' .%0.2d',round(pxc(2:nlcount,kidx).*100));
   fprintf('%.2f ',pxc(2:nlcount,kidx));
   
   if OUTNLMODE==4,
      for nlidx=1:size(fracxc,2),
         if goodcheck(5,kidx),
            fprintf('%.2f ',fracxc(kidx,nlidx));
         else
            fprintf('     ');
         end
      end
   end
      
   fprintf('%d',xccount(1:nlcount,kidx));
   
   fprintf('\n');
end

if nargout>0,
   res.cellid=cellid;
   res.sigflag=squeeze(xccount)';
   res.p=squeeze(pxc)';
   res.predxc=permute(predxc(1:2,:,:),[3 2 1]);
   res.goodcheck=goodcheck';
   return
end


mpxcin=sqrt(nanmean(squeeze(predxc(1,:,:).^2)'));
mpxcout=sqrt(nanmean(squeeze(predxc(2,:,:).^2)'));

m2pxcin=squeeze(predxc(1,:,:));
m2pxcin(find(~goodcheck))=nan;
m2pxcin=sqrt(nanmean((m2pxcin.^2)'));
m2pxcout=squeeze(predxc(2,:,:));
m2pxcout(find(~goodcheck))=nan;
m2pxcout=sqrt(nanmean((m2pxcout.^2)'));

if OUTNLMODE~=6 & exist('predinf','var'),
   m3pxcin=squeeze(predinf(1,:,:));
   m3pxcin(find(~goodcheck))=nan;
   m3pxcin=sqrt(nanmean((m3pxcin.^2)'));
   m3pxcout=squeeze(predinf(2,:,:));
   m3pxcout(find(~goodcheck))=nan;
   m3pxcout=sqrt(nanmean((m3pxcout.^2)'));
end

fprintf('mean :');
fprintf(' .%0.2d/.%0.2d',round(cat(1,mpxcin,mpxcout).*100));
fprintf('\n');
fprintf('modn :');
fprintf(' %3d    ',sum(xccount,3));
fprintf('\n');
fprintf('tmod%%:');
fprintf('   %.1f ',sum(xccount,3)./size(predxc,3).*100);
fprintf('\n');

fprintf('gmean:');
fprintf(' .%0.2d/.%0.2d',round(cat(1,m2pxcin,m2pxcout).*100));
if OUTNLMODE==4,
   fprintf('                       ');
   fprintf('%.2f ',nanmean(fracxc));
end
fprintf('\n');
fprintf('goodn:');
fprintf(' %3d    ',sum(goodcheck,2));
fprintf(' /%3d\n',size(predxc,3));
fprintf('gmod%%:');
if OUTNLMODE==5,
   fprintf('   %.1f ',100.*sum(xccount,3)./...
           sum(goodcheck | ...
               repmat(goodcheck(1,:),[size(goodcheck,1) 1]),2));
else
   fprintf('   %.1f ',sum(xccount,3)./sum(goodcheck,2).*100);
end
fprintf('\n');
if OUTNLMODE~=6 & exist('predinf','var'),
   fprintf('ginf :');
   fprintf(' .%0.2d/.%0.2d',round(cat(1,m3pxcin,m3pxcout).*100));
   fprintf('\n');
end

if 1,
   goodcheck(:)=1;
   disp('including all cells (even bad ones!) in summary figures');
end

if OUTNLMODE==4,
   rac=[2 3 5];
elseif OUTNLMODE==5,
   rac=[2 3 4];
else
   rac=[2 3 6];
end

fprintf('olap: ');
%rac=[4 5 6];
if OUTNLMODE==5,
   vdata=zeros(2,2);
   cn=0;
   for ii=0:1,
      for jj=0:1,
         n=sum((xccount(rac(1),:)==ii) & ...
               (xccount(rac(2),:)==jj));
         fprintf('%d%d: %d ',ii,jj,n);
         cn=cn+n;
         vdata(ii+1,jj+1)=n;
      end
   end
   fprintf('tot: %3d\n',cn);
   
else
   vdata=zeros(2,2,2);
   cn=0;
   for ii=0:1,
      for jj=0:1,
         for kk=0:1,
            n=sum((xccount(rac(1),:)==ii) & ...
                  (xccount(rac(2),:)==jj) & ...
                  (xccount(rac(3),:)==kk) & goodcheck(rac(1),:));
            fprintf('%d%d%d: %d ',ii,jj,kk,n);
            cn=cn+n;
            vdata(ii+1,jj+1,kk+1)=n;
      end
      end
   end
   fprintf('tot: %3d\n',cn);
end

figure(1)
clf

for nlidx=2:length(nlnames),
   subplot(3,3,nlidx-1);
   
   if OUTNLMODE==4,
      
      bb=find(~xccount(nlidx,:) & goodcheck(nlidx,:) & ...
              ~isnan(fracxc(:,1)') & fracxc(:,1)'<1 & ...
              max(squeeze(predxc(:,nlidx,:) > 0),[],1));
      aa=find(xccount(nlidx,:) & goodcheck(nlidx,:) & ...
              ~isnan(fracxc(:,1)') & fracxc(:,1)'<1 & ...
              max(squeeze(predxc(:,nlidx,:) > 0),[],1));
      pa=(predxc(1,nlidx,aa).^2-predxc(2,nlidx,aa).^2)./...
         max(predxc(:,nlidx,aa),[],1).^2;
      pb=(predxc(1,nlidx,bb).^2-predxc(2,nlidx,bb).^2)./...
         max(predxc(:,nlidx,bb),[],1).^2;
      %pa=pa./sqrt(abs(pa+(pa==0)));
      %pb=pb./sqrt(abs(pb+(pb==0)));
      arange=[-1.2 1.2 0 30];
   else
      bb=find(~xccount(nlidx,:) & goodcheck(nlidx,:));
      aa=find(xccount(nlidx,:) & goodcheck(nlidx,:));
      pa=predxc(1,nlidx,aa);
      pb=predxc(1,nlidx,bb);
      arange=[-0.3 0.3 0 30];
   end
   histcomp(pa,pb,'no att',nlnames{nlidx},'xc',arange);
   hold on
   plot([0 0],[arange(3) arange(4)],'k--');
   hold off
end

hs=subplot(3,3,5);
tf=fracxc(find(~isnan(fracxc(:,1)) & fracxc(:,1)<1),:);

if 0,
   tf=nanmean(fliplr(tf));
   sleg={};
   rrange=fliplr([1 rac]);
   for ii=1:length(rrange),
      sleg{ii}=sprintf('%s %.2f',nlnames{rrange(ii)},tf(ii));
   end
   
   pie(tf+0.001,sleg);
   title('pct explainable variance');
else
   bar(nanmean(tf));
   xticks(1:length(rac)+1,{nlnames{[1 rac]}});
   ylabel('frac pred var');
   axis([0 length(rac)+2 0 1]);
   axis square
end

hs=subplot(3,3,6);
if OUTNLMODE==5,
   counts=sum([xccount(2,1,:) ...
               ~xccount(2,1,:)],3);
   pie(counts+0.001);
elseif OUTNLMODE==4,
   % only show local , gain, dc and not (not dcg)
   gidx=find(goodcheck(5,:));
   counts=sum([xccount(5,1,gidx) ...
               ~xccount(5,1,gidx) & xccount(3,1,gidx) ...
               ~xccount(5,1,gidx) & ~xccount(3,1,gidx) & xccount(2,1,gidx) ...
               ~xccount(5,1,gidx) & ~xccount(3,1,gidx) & ~xccount(2,1,gidx)],3);
   pie(counts+0.001,...
       {['local ',num2str(counts(1))],['gain ',num2str(counts(2))],...
        ['DC ',num2str(counts(3))],['None ',num2str(counts(4))]});
   %counts=sum([xccount(5,1,:) ...
   %            ~xccount(5,1,:) & xccount(3,1,:) & xccount(2,1,:) ...
   %            ~xccount(5,1,:) & ~xccount(2,1,:) & xccount(3,1,:) ...
   %            ~xccount(5,1,:) & ~xccount(3,1,:) & xccount(2,1,:) ...
   %            ~xccount(5,1,:) & ~xccount(3,1,:) & ~xccount(2,1,:)],3);
   %pie(counts+0.001,...
   %    {['local ',num2str(counts(1))],['g+dc ',num2str(counts(2))],...
   %     ['glob ',num2str(counts(3))],['DC ',num2str(counts(4))],...
   %     ['None ',num2str(counts(5))]});
else
   counts=sum([xccount(6,1,:) ...
               ~xccount(6,1,:) & xccount(3,1,:) & xccount(2,1,:) ...
               ~xccount(6,1,:) & ~xccount(2,1,:) & xccount(3,1,:) ...
               ~xccount(6,1,:) & ~xccount(3,1,:) & xccount(2,1,:) ...
               ~xccount(6,1,:) & ~xccount(4,1,:) & ~xccount(3,1,:) ...
               & ~xccount(2,1,:)],3);
   pie(counts+0.001,...
       {['local ',num2str(counts(1))],['g+dc ',num2str(counts(2))],...
        ['glob ',num2str(counts(3))],['DC ',num2str(counts(4))],...
        ['None ',num2str(counts(5))]});
end

title('sig cell counts');

hs=subplot(3,3,7);
hist(pxc(rac(1),1,:),20);
title([nlnames{rac(1)},' pred p']);
axis square

hs=subplot(3,3,8);
hist(pxc(rac(2),1,:),20);
title([nlnames{rac(2)},' pred p']);
axis square

hs=subplot(3,3,9);
hist(pxc(rac(3),1,:),20);
title([nlnames{rac(3)},' pred p']);
axis square

fullpage('portrait');
colormap(gray);

figure(2);
clf
venn(vdata,1);
legend(sprintf('%s (%d)',nlnames{rac(1)},sum(xccount(rac(1),:))),...
       sprintf('%s (%d)',nlnames{rac(2)},sum(xccount(rac(2),:))),...
       sprintf('%s (%d)',nlnames{rac(3)},sum(xccount(rac(3),:))));
title('overlap between effects');





if OUTNLMODE==4,
   figure(3);
   clf
   
   %totset=find(xccount(2,:) | xccount(3,:) | xccount(4,:) | xccount(5,:));
   totset=find(xccount(4,:) | xccount(5,:));
   txc0=predxc(1:2,2:5,totset);
   %txc0(2,1:3,:)=repmat(predxc(2,1,totset),[1 3]);
   txc0(find(txc0<0))=0;
   txc0=txc0.^2;
   
   txcmax=squeeze(max(txc0(1,:,:),[],2));
   txcmax(txcmax==0)=1;
   txc=squeeze(txc0(1,:,:)-txc0(2,:,:))' ./ repmat(txcmax,[1 4]);
   
   arange=[-0.85 0.85 0 25];
   
   for nlidx=1:4,
      subplot(2,2,nlidx);
      bb=find(~xccount(nlidx+1,totset));
      aa=find(xccount(nlidx+1,totset));
      pa=txc(aa,nlidx);
      pb=txc(bb,nlidx);
      tt=sprintf('med %.2f',median([pa;pb]));
      histcomp(pa,pb,nlnames{nlidx+1},tt,'xc',arange);
      hold on
      plot([0 0],[arange(3) arange(4)],'k--');
      hold off
   end
   colormap(gray);
   
   keyboard
end


if 0,
   figure(3);
   
   subplot(1,2,1);
   tf=fracxc(find(~isnan(fracxc(:,1)) & fracxc(:,1)<1),:);
   
   tf(tf==0)=-0.1;
   [x,n]=hist(tf(:,2:end));
   bar(n,x,1.2,'grouped');
   %tf=[tf(:,2:end) tf(:,1)];
   %bar(flipud(sortrows(tf)),'stacked');
   title('pct explainable variance per att model');
   
   subplot(1,2,2);
   
   %nlidx=[2 5];
   %aa=find(goodcheck(5,:));
   %pa=squeeze((predxc(1,nlidx,aa).^2-predxc(2,nlidx,aa).^2)./...
   %   predxc(1,nlidx,aa).^2)';
   %scatter(pa(:,1),pa(:,2),'.');
   tf=fracxc(find(~isnan(fracxc(:,1)) & fracxc(:,1)<1),:);
   scatter(tf(:,2),tf(:,4),'.');
   axis square
   axis([-0.1 1 -0.1 1]);
   title('dc vs local pct pred');
   
   %keyboard
   
end

if dump,
   fprintf('dumping summary to eps files');
   outpath=sprintf('/auto/k5/david/data/batch%d/output',batchid);
   if ~exist(outpath,'dir'),
      unix(['mkdir ',outpath]);
   end
   cd(outpath);
   drawnow
   print -f1 -depsc summary.eps
   print -f2 -depsc venn.eps
   
   for ii=1:length(cellids),
      
      if strcmp(batchdata.matcmd,'xcdms'),
         dmsres(cellids{ii},batchid);
      else
         kvares(cellids{ii},batchid);
      end
      
      drawnow
      print('-f2','-depsc',sprintf('%s.%d.kern.eps',cellids{ii},batchid));
   end
end



return

keyboard

cellid={preddata(ttt).cellid};
cellcount=length(cellid);
cellcount=cellcount+1;
cellid{cellcount}='mean';

runidx=find(okidx);
runidx(cellcount)=0;

predxc=cat(5,preddata(ttt).predxc);
tpredxc=predxc;
tpredxc(find(isnan(tpredxc)))=0;
predxc=cat(5,predxc,mean(tpredxc,5));
%attsimxc=cat(4,preddata(ttt).attsimxc);
%attsimxc=cat(4,attsimxc,mean(attsimxc,4));

TSON=0;
if TSON,
   targsim=cat(4,preddata(ttt).targsim);
   targsim=cat(4,targsim,mean(targsim,4));
end

% make the last entry in each summary statistic the mean across
% cells
if length(preddata(ttt(1)).p)>1,
   p=repmat(zeros(size(preddata(ttt(1)).p)),[1 1 1 length(ttt)]);
   for ii=1:length(ttt),
      pcount0=min([size(preddata(ttt(ii)).p,1) size(p,1)]);
      pcount=min([size(preddata(ttt(ii)).p,2) size(p,2)]);
      pcount2=min([size(preddata(ttt(ii)).p,3) size(p,3)]);
      p(1:pcount0,1:pcount,1:pcount2,ii)=...
          preddata(ttt(ii)).p(1:pcount0,1:pcount,1:pcount2);
   end
   p=cat(4,p,mean(p,4));
else
   p=zeros(size(predxc));
end

% display prediction scores for this latency bin (4 means 0-100 ms after
% fixation onset)
latcount=size(predxc,3);
latbase=min([latcount 1])-1;
if ~exist('latidx','var'),
   latidx=min([latcount 6]);
end
PCLOW=10;
if ~exist('pccount','var'),
   pccount=18
end
nlidx=1;

% PTHRESH=0.012; % p < 1-(1-0.05)^(1/4)
PTHRESH=0.008; % p < 1-(1-0.05)^(1/6
%PTHRESH=0.05
figure(1);
clf
minx=min(min(min(predxc(:,:,latidx,:))));
maxx=max(max(max(predxc(:,:,latidx,:))));
if minx<-maxx & maxx>0,
   minx=-maxx;
end
for ii=1:cellcount,
   subplot(ceil(cellcount/4),4,ii);
   imagesc(predxc(:,:,latidx,nlidx,ii),[minx maxx]);
   axis image
   axis off
   title(sprintf('%s (%d)',cellid{ii},runidx(ii)));
end
colorbar
colormap(hot);

attcount=size(predxc,1);
attpcount=size(p,1);

disp('cell  SMT (pred in/out) (att comp 1-2 1-3 1-4 2-3 2-4 3-4)');
pxc=zeros(cellcount,attcount);
pxcout=zeros(cellcount,attcount);
pxcatt=[];
sigcount=[0 0];
attpairs=[1 2; 1 3; 1 4; 2 3; 2 4; 3 4];

xcallin=[];
xcallout=[];
pall=[];

attbetter=zeros(cellcount-1,attcount-1);
for ii=1:cellcount,
   fprintf('%-5s',cellid{ii});
   for attidx=1:attcount,
      %tpredxc=predxc(:,:,9,nlidx,ii);
      tpredxc=predxc(:,:,latidx,nlidx,ii);
      %tpredxc(find(isnan(tpredxc)))=0;
      pxc(ii,attidx)=tpredxc(attidx,attidx);
      pxcout(ii,attidx)=mean(tpredxc(attidx,[2:attidx-1 (attidx+1):end]));
      
      if attidx>1 & ii<cellcount,
         xcallin=[xcallin tpredxc(attidx,attidx).*ones(1,attcount-2)];
         xcallout=[xcallout tpredxc(attidx,[2:attidx-1 (attidx+1):end])];
      end
      
      pxcout(ii,attidx)=mean(tpredxc(attidx,[2:attidx-1 (attidx+1):end]));
      fprintf(' %5.2f/%5.2f',pxc(ii,attidx),pxcout(ii,attidx));
      
      if ii<cellcount & attidx>1 & pxc(ii,attidx)>pxcout(ii,attidx),
         attbetter(ii,attidx-1)=1;
      end
   end
   
   % print out in vs out pred xcs. transfer over to pairwise
   % comparison if it looks like it's gonna help.
   
   if preddata(ttt(1)).p~=0,
      for attidx=1:attpcount,
         pall=[pall; p(attidx,pccount,latidx-latbase,ii)];
         if p(attidx,pccount,latidx-latbase,ii)<=PTHRESH,
            sextra='*';
            %pxcatt[pxcatt; pxc(ii,attidx) pxcout(ii,attidx)];
            pxcatt=[pxcatt; pxc(ii,attpairs(attidx,1)+1) ...
                    pxc(ii,attpairs(attidx,2)+1)];
         else
            sextra=' ';
         end
         fprintf('%s%5.3f',sextra,p(attidx,pccount,latidx-latbase,ii));
      end
      
      % if qualifies, mark as cell with significant attentional effect
      if min(p(1:attpcount,PCLOW,latidx-latbase,ii))<=PTHRESH,
         fprintf('*%d',PCLOW);
         sigcount(2)=sigcount(2)+1;
      end
      if min(p(1:attpcount,pccount,latidx-latbase,ii))<=PTHRESH,
         fprintf('*%d',pccount);
         sigcount(1)=sigcount(1)+1;
      end
   else
      pxcatt=[0 0];
   end
   fprintf('\n');
end
fprintf(['sig att: %d/%d %.2f%% (%dPC) %d/%d %.2f%% (%dPC)',...
         ' p<%.3f respidx=%d\n'],...
        sigcount(2),cellcount-1,sigcount(2)./(cellcount-1).*100,...
        PCLOW,...
        sigcount(1),cellcount-1,sigcount(1)./(cellcount-1).*100,...
        pccount,PTHRESH,latidx);
xcallin=xcallin(find(~isnan(xcallout)));
xcallout=xcallout(find(~isnan(xcallout)));
[pp,mm,ss]=randpairtest(xcallin,xcallout,2000);

fprintf('mean xc in/out: %.3f/%.3f (p<%.3f): %d/%d (%.1f%%) attstates\n',...
        mean(xcallin),mean(xcallout),pp,...
        sum(xcallin>xcallout),length(xcallout),...
        sum(xcallin>xcallout)/length(xcallout)*100);

figure(2);
clf
subplot(2,3,1);
scatter(xcallout(:),xcallin(:));
hold on
plot([0 1.0],[0 1.0],'k--');
plot(nanmean(xcallout),nanmean(xcallin),'r+','LineWidth',2);
hold off
axis equal
xlabel('pcx out');
ylabel('pcx in');
title('attention preds in vs. out');

subplot(2,3,2);
scatter(mean(pxcout(:,2:end),2),pxc(:,1));
hold on
plot([0 1.0],[0 1.0],'k--');
plot(median(mean(pxcout(:,2:end),2)),median(pxc(:,1)),...
     'r+','LineWidth',2);
hold off
axis equal
xlabel('pcx att in');
ylabel('pcx all att');
title('data lim: all vs. atten preds');

subplot(2,3,3);
mxcin=zeros(latcount,1);
mxcout=zeros(latcount,1);
for attidx=2:attcount,
   mxcin=mxcin+squeeze(predxc(attidx,attidx,:,nlidx,end)) ./ (attcount-1);
   mxcout=mxcout+squeeze(sum(predxc(attidx,[2:attidx-1 (attidx+1):end],...
                                     :,nlidx,end),2)) ./ ...
          ((attcount-2)*(attcount-1));
end
plot([mxcout mxcin]);
legend('out','in');
%scatter(pxcatt(:,2),pxcatt(:,1));
%hold on
%plot([0 1.0],[0 1.0],'k--');
%plot(median(pxcatt(:,2)),median(pxcatt(:,1)),...
%     'r+','LineWidth',2);
%hold off
%axis equal
%xlabel('xc att out');
%ylabel('xc att in');
title(sprintf('att mod preds (latidx=%d,pccount=%d)',latidx,pccount));

subplot(2,3,4);
%plot(squeeze(sum((min(p(2:end,:,:,:),[],1)<=PTHRESH),4)),'LineWidth',2);
plot(squeeze(sum((min(p(1:end,:,:,:),[],1)<=PTHRESH),4)),'LineWidth',2);
hold on
%plot(squeeze(sum((min(p(2:end,:,latidx-latbase,:),[],1)<=PTHRESH),4)),...
%     'kx','LineWidth',2);
plot(squeeze(sum((min(p(1:end,:,latidx-latbase,:),[],1)<=PTHRESH),4)),...
     'kx','LineWidth',2);
hold off
xlabel('pc count');
ylabel('N signif atten');
legend('lat1','lat2');
title('sig att cells vs. pc count');

subplot(2,3,5);
hist(pall,15);
title('pair p values')

subplot(2,3,6);
hist(xcallin(:)-xcallout(:),15);
title('pair pred diff (in - out)');

if TSON,
   subplot(2,3,5);
   ts1=targsim(:,pccount,1,:);
   ts2=targsim(:,pccount,2,:);
   tp=p(:,pccount,latidx-latbase,:);
   sigsim=[ts1(:) ts2(:) tp(:)];
   scatter(sigsim(:,1),sigsim(:,3));
   title('att diff vs targ diff (unnorm)');
   xlabel('targ1.*targ2');
   ylabel('att p');
   
   subplot(2,3,6);
   scatter(sigsim(find(tp<PTHRESH),2),sigsim(find(tp<PTHRESH),3));
   title('att diff vs targ diff (norm)');
   xlabel('angle(targ1,targ2)');
   ylabel('att p');
end

%keyboard
return

figure(3);
clf
atxc=squeeze(attsimxc(:,:,pccount,:));
minx=min(atxc(:));
maxx=max(atxc(:));
if minx<-maxx & maxx>0,
   minx=-maxx;
end
for ii=1:cellcount-1,
   subplot(ceil(cellcount/4),4,ii);
   imagesc(atxc(:,:,ii),[minx maxx]);
   axis image
   axis off
   title(sprintf('atxc %s (%d)',cellid{ii},runidx(ii)));
end
colorbar
colormap(hot);



disp('cell  T   (pred in/out)');
for ii=1:cellcount,
   fprintf('%-5s',cellid{ii});
   for attidx=1:attcount,
      pxc=predxc(attidx,attidx,2,ii);
      pxcout=mean(predxc(attidx,[2:attidx-1 (attidx+1):end],latidx,ii));
      fprintf(' %5.3f (%4.1f/%4.1f)',p(attidx,1,1,ii),pxc,pxcout);
   end
   
   % if qualifies, mark as cell with significant attentional effect
   if min(p(:,1,1,ii))<PTHRESH,
      fprintf('*');
   end
   fprintf('\n');
end


