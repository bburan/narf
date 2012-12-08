% v1predsum.m

% callmaster : exp X in-model X cell
% predxc: exp X in-model X cell X out-model X cnf

%function v1predsum;

dbopen;

%batches=[24 26 23; ...   % rev: pix, pfft, psfft
%         27 28 29; ...   % gr:  ""
%         30 31 32];      % nr:  ""
batches=[24 26; ...   % rev: pix, pfft
         27 28; ...   % gr:  ""
         30 31];      % nr:  ""
%batches=[24 23; ...   % rev: pix, pfft
%         27 29; ...   % gr:  ""
%         30 32];      % nr:  ""


%close all
for ii=1:2,
   figure(ii);
end
drawnow;

% load pred results for each batch
res=loadbigdata(batches);
predxc=res.predxc;
prederr=res.prederr;
predp=res.predp;
celllist=res.celllist;

badcellidx=find(strcmp(celllist,'R150B') | strcmp(celllist,'93G83A'));
predxc(badcellidx,1,:,:,:)=nan;


MINLEN=500;
PMIN=0.01;
expstr={'Rev','GR','NR'};
instr={'pix','pow','psf'};
cnfstr=expstr;
outstr={'lin','thr','ethr','full'};

% review-only model comparison
figure(1);
clf
expidx1=[1 3 2 ];
expidx2=[1 3 2 ];
nlinidx1=[1 1 1 ]; 
nlinidx2=[2 2 2 ];
cnfidx=[1 3 2 ];
nloutidx1=[2 2 2 ];
nloutidx2=[2 2 2 ];
%expidx1=[1 2 3; 2 1 3; 3 1 2];
%expidx2=[1 1 1; 2 2 2; 3 3 3];
%nlinidx1=[1 1 1; 1 1 1; 1 1 1]; 
%nlinidx2=[2 1 1; 2 1 1; 2 1 1];
%cnfidx=[1 1 1; 2 2 2; 3 3 3];
%nloutidx1=[5 5 5; 5 5 5; 5 5 5];   % alt sub 9 for 3, 5 for 2
%nloutidx2=[5 5 5; 5 5 5; 5 5 5];

rowcount=size(expidx1,1);
colcount=size(expidx2,2);
for ii=1:length(expidx1(:)),
   subplot(rowcount,colcount,ii);
   n1=[expstr{expidx1(ii)},'-',cnfstr{cnfidx(ii)},'-',instr{nlinidx1(ii)},...
       '-',outstr{nloutidx1(ii)}];
   n2=[expstr{expidx2(ii)},'-',cnfstr{cnfidx(ii)},'-',instr{nlinidx2(ii)},...
       '-',outstr{nloutidx2(ii)}];
   set1=predxc(:,expidx1(ii),nlinidx1(ii),cnfidx(ii),nloutidx1(ii));
   set2=predxc(:,expidx2(ii),nlinidx2(ii),cnfidx(ii),nloutidx2(ii));
   goodidx=find(~isnan(set1) & ~isnan(set2));
   %goodidx=find(~isnan(set1) & ~isnan(set2) & ...
   %        totlen(:,expidx1(ii),1)>MINLEN & totlen(:,expidx2(ii),1)>MINLEN);
   set1=set1(goodidx);
   set2=set2(goodidx);
   p1=predp(goodidx,expidx1(ii),nlinidx1(ii),cnfidx(ii),nloutidx1(ii));
   p2=predp(goodidx,expidx2(ii),nlinidx2(ii),cnfidx(ii),nloutidx2(ii));
   e1=prederr(goodidx,expidx1(ii),nlinidx1(ii),cnfidx(ii),nloutidx1(ii));
   e2=prederr(goodidx,expidx2(ii),nlinidx2(ii),cnfidx(ii),nloutidx2(ii));
   
   set1hi=find(set1-e1 > set2+e2);
   set2hi=find(set1+e1 < set2-e2);
   sigidx=ones(size(set1))*2;
   sigidx(set1hi)=3;
   sigidx(set2hi)=1;
   
   %sigidx=2-(set1(goodidx)+e1(goodidx) < set2(goodidx)-e2(goodidx));
   %sigidx=2-(p1<PMIN | p2<PMIN);
   
   [p,m]=plotcomp(set1,set2,n1,n2,[-0.4 1.0 -0.4 1.0],sigidx);

   fprintf('%s v %s: %.2f v %.2f (p<%.3f) %d > %d < %d / %d \n',...
           n1,n2,mean(set1),mean(set2),p(end),sum(sigidx==3),...
           sum(sigidx==2),sum(sigidx==1),length(sigidx))
   fprintf('                    %.2f v %.2f\n',...
           mean(set1.*abs(set1)),mean(set2.*abs(set2)))
end
set(gcf,'PaperPosition',[0.25 0.25 8 10.5],'PaperOrientation','portrait');
drawnow


% within class analysis
figure(2)
clf
%expidx1=[3 2];
%expidx2=[1 1];
expidx1=[3 2 2];
expidx2=[1 1 3];
nlinidx1=[2 2 2]; 
nlinidx2=[2 2 2];
cnfidx=[1 1 1];
nloutidx1=[2 2 2];
nloutidx2=[2 2 2];

rowcount=size(expidx1,1)+1;
colcount=size(expidx2,2);
for ii=1:length(expidx1(:)),
   subplot(rowcount,colcount,ii);
   n1=[expstr{expidx1(ii)},'-',cnfstr{cnfidx(ii)},'-',instr{nlinidx1(ii)},...
       '-',outstr{nloutidx1(ii)}];
   n2=[expstr{expidx2(ii)},'-',cnfstr{cnfidx(ii)},'-',instr{nlinidx2(ii)},...
       '-',outstr{nloutidx2(ii)}];
   set1=predxc(:,expidx1(ii),nlinidx1(ii),cnfidx(ii),nloutidx1(ii));
   set2=predxc(:,expidx2(ii),nlinidx2(ii),cnfidx(ii),nloutidx2(ii));
   goodidx=find(~isnan(set1) & ~isnan(set2));
   %goodidx=find(~isnan(set1) & ~isnan(set2) & ...
   %        totlen(:,expidx1(ii),1)>MINLEN & totlen(:,expidx2(ii),1)>MINLEN);
   set1=set1(goodidx);
   set2=set2(goodidx);
   p1=predp(goodidx,expidx1(ii),nlinidx1(ii),cnfidx(ii),nloutidx1(ii));
   p2=predp(goodidx,expidx2(ii),nlinidx2(ii),cnfidx(ii),nloutidx2(ii));
   e1=prederr(goodidx,expidx1(ii),nlinidx1(ii),cnfidx(ii),nloutidx1(ii));
   e2=prederr(goodidx,expidx2(ii),nlinidx2(ii),cnfidx(ii),nloutidx2(ii));
   
   set1hi=find(set1-e1 > set2+e2);
   set2hi=find(set1+e1 < set2-e2);
   sigidx=ones(size(set1))*2;
   sigidx(set1hi)=3;
   sigidx(set2hi)=1;
   
   %sigidx=2-(set1(goodidx)+e1(goodidx) < set2(goodidx)-e2(goodidx));
   %sigidx=2-(p1<PMIN | p2<PMIN);
   
   [p,m]=plotcomp(set1,set2,n1,n2,[-0.4 1.0 -0.4 1.0],sigidx);
   
   fprintf('%s v %s: %.2f v %.2f (p<%.3f) %d > %d < %d / %d \n',...
           n1,n2,mean(set1),mean(set2),p(end),sum(sigidx==3),...
           sum(sigidx==2),sum(sigidx==1),length(sigidx))
   fprintf('                    %.2f v %.2f\n',...
           mean(set1.*abs(set1)),mean(set2.*abs(set2)))
end
set(gcf,'PaperPosition',[0.25 0.25 8 10.5],'PaperOrientation','portrait');
drawnow

hs=subplot(rowcount,colcount,length(expidx1(:))+1);
nlinidx=2; nloutidx=2;
set1=squeeze(predxc(:,[1 3 2],nlinidx,[1 3 2],nloutidx));
goodidx=find(sum(isnan(set1(:,:,1)),2)==0);
m=squeeze(mean(set1(goodidx,:,:),1))';
b1=bar(m);
set(b1(1),'CDataMapping','direct','FaceVertexCData',zeros(16,1)+0)
set(b1(2),'CDataMapping','direct','FaceVertexCData',zeros(16,1)+25)
set(b1(3),'CDataMapping','direct','FaceVertexCData',zeros(16,1)+50)
axis([0 4 0 0.6]);
axis square
title(sprintf('power domain cross preds (n=%d)',length(goodidx)));
set(hs,'XTickLabel',{'Nat','DNat','DGrat'});
xlabel('validation class');
disp('legend turned off');
%legend('Nat','DNat','DGrat');

hs=subplot(rowcount,colcount,length(expidx1(:))+2);
nlinidx=1; nloutidx=2;
set1=squeeze(predxc(:,[1 3 2],nlinidx,[1 3 2],nloutidx));
goodidx=find(sum(isnan(set1(:,:,1)),2)==0);
m=squeeze(mean(set1(goodidx,:,:),1))';
b1=bar(m);
set(b1(1),'CDataMapping','direct','FaceVertexCData',zeros(16,1)+0)
set(b1(2),'CDataMapping','direct','FaceVertexCData',zeros(16,1)+25)
set(b1(3),'CDataMapping','direct','FaceVertexCData',zeros(16,1)+50)
axis([0 4 0 0.6]);
axis square
title(sprintf('pixel domain cross preds (n=%d)',length(goodidx)));
set(hs,'XTickLabel',{'Nat','DNat','DGrat'});
xlabel('validation class');

set(gcf,'PaperPosition',[0.25 2.5 8 6],'PaperOrientation','portrait')
colormap(gray);
drawnow

%
% find all matches that pred rev
set1=predxc(:,:,:,:,2);
all3goodidx=find(sum(isnan(set1(:,:,1,1)),2)==0);
set1=set1(all3goodidx,:,:,:);
fprintf('%d cells have all stim classes\n',length(all3goodidx));
fprintf('Exp Dom Val\n');
for nlinidx=1:2,
   for expidx=1:3,
      for cnfidx=1:3,
         fprintf('%3s %3s %3s  r=%.3f  r2=%.3f\n',...
                 expstr{expidx},instr{nlinidx},expstr{cnfidx},...
                 mean(set1(:,expidx,nlinidx,cnfidx)),...
                 mean(set1(:,expidx,nlinidx,cnfidx).* ...
                      abs(set1(:,expidx,nlinidx,cnfidx))));
      end
   end
end


if 1,
nloutidx=2;
nlinidx=2;
disp('              within class...            Nat pred            Pix-w/in pred ...');
disp('cellid:     Nat   DGrat    DNat     Nat   DGrat    DNat     Nat   DGrat    DNat');
for cellidx=1:length(celllist),
   fprintf('%6s:',celllist{cellidx});
   for expidx=1:3,
      fprintf(' %7.3f',predxc(cellidx,expidx,nlinidx,expidx,nloutidx));
   end
   for expidx=1:3,
      fprintf(' %7.3f',predxc(cellidx,expidx,nlinidx,1,nloutidx));
   end
   for expidx=1:3,
      fprintf(' %7.3f',predxc(cellidx,expidx,1,expidx,nloutidx));
   end
   fprintf('\n');
   
end
end

return

% within class analysis
figure(3)
clf
expidx1=[2 3 3; 1 3 3; 1 2 2];
expidx2=[1 1 2; 2 2 1; 3 3 1];
nlinidx1=[2 2 2; 2 2 2; 2 2 2]; 
nlinidx2=[2 2 2; 2 2 2; 2 2 2];
cnfidx=[1 1 1; 2 2 2; 3 3 3];
nloutidx1=[5 5 5; 5 5 5; 5 5 5];
nloutidx2=[5 5 5; 5 5 5; 5 5 5];

rowcount=size(expidx1,1);
colcount=size(expidx2,2);
for ii=1:length(expidx1(:)),
   subplot(rowcount,colcount,ii);
   n1=[expstr{expidx1(ii)},'-',cnfstr{cnfidx(ii)},'-',instr{nlinidx1(ii)},...
       '-',outstr{nloutidx1(ii)}];
   n2=[expstr{expidx2(ii)},'-',cnfstr{cnfidx(ii)},'-',instr{nlinidx2(ii)},...
       '-',outstr{nloutidx2(ii)}];
   set1=predxc(:,expidx1(ii),nlinidx1(ii),cnfidx(ii),nloutidx1(ii));
   set2=predxc(:,expidx2(ii),nlinidx2(ii),cnfidx(ii),nloutidx2(ii));
   p1=predp(:,expidx1(ii),nlinidx1(ii),cnfidx(ii),nloutidx1(ii));
   p2=predp(:,expidx2(ii),nlinidx2(ii),cnfidx(ii),nloutidx2(ii));
   
   goodidx=find(~isnan(set1) & ~isnan(set2));
   sigidx=2-(p1(goodidx)<PMIN | p2(goodidx)<PMIN);
   %goodidx=find(~isnan(set1) & ~isnan(set2) & ...
   %        totlen(:,expidx1(ii),1)>MINLEN & totlen(:,expidx2(ii),1)>MINLEN);
   
   plotcomp(set1(goodidx),set2(goodidx),n1,n2,[-0.4 1.0 -0.4 1.0],sigidx);
end
set(gcf,'PaperPosition',[0.25 0.25 10.5 8],'PaperOrientation','landscape');
drawnow

return

% review-only model comparison
figure(1);
clf
expidx1=[1 1 1; 1 1 1; 1 1 1];
expidx2=[1 1 1; 1 1 1; 1 1 1];
nlinidx1=[1 2 1; 3 3 3; 1 2 1]; 
nlinidx2=[2 3 3; 3 3 3; 2 3 3];
cnfidx=[1 1 1; 1 1 1; 1 1 1];
nloutidx1=[9 9 9; 1 5 1; 5 5 5];   % alt sub 9 for 3, 5 for 2
nloutidx2=[9 9 9; 5 9 9; 5 5 5];
%nloutidx1=[3 3 3; 1 2 1; 1 1 1];
%nloutidx2=[3 3 3; 2 3 3; 1 1 1];

rowcount=size(expidx1,1);
colcount=size(expidx2,2);
for ii=1:length(expidx1(:)),
   subplot(rowcount,colcount,ii);
   n1=[expstr{expidx1(ii)},'-',cnfstr{cnfidx(ii)},'-',instr{nlinidx1(ii)},...
       '-',outstr{nloutidx1(ii)}];
   n2=[expstr{expidx2(ii)},'-',cnfstr{cnfidx(ii)},'-',instr{nlinidx2(ii)},...
       '-',outstr{nloutidx2(ii)}];
   set1=predxc(:,expidx1(ii),nlinidx1(ii),cnfidx(ii),nloutidx1(ii));
   set2=predxc(:,expidx2(ii),nlinidx2(ii),cnfidx(ii),nloutidx2(ii));
   p1=predp(:,expidx1(ii),nlinidx1(ii),cnfidx(ii),nloutidx1(ii));
   p2=predp(:,expidx2(ii),nlinidx2(ii),cnfidx(ii),nloutidx2(ii));
   
   goodidx=find(~isnan(set1) & ~isnan(set2));
   sigidx=2-(p1(goodidx)<PMIN | p2(goodidx)<PMIN);
   %goodidx=find(~isnan(set1) & ~isnan(set2) & ...
   %        totlen(:,expidx1(ii),1)>MINLEN & totlen(:,expidx2(ii),1)>MINLEN);
   
   plotcomp(set1(goodidx),set2(goodidx),n1,n2,[-0.4 1.0 -0.4 1.0],sigidx);
end
%set(gcf,'PaperPosition',[0.25 0.25 8 10.5],'PaperOrientation','portrait');
set(gcf,'PaperPosition',[0.25 0.25 10.5 8],'PaperOrientation','landscape');
drawnow


% within class analysis
figure(2)
clf
expidx1=[1 2 3; 1 2 3; 1 2 3];
expidx2=[1 2 3; 1 2 3; 1 2 3];
nlinidx1=[1 1 1; 2 2 2; 3 3 3]; 
nlinidx2=[2 2 2; 3 3 3; 3 3 3];
cnfidx=[1 2 3; 1 2 3; 1 2 3];
nloutidx1=[9 9 9; 9 9 9; 5 5 5];   % alt sub 9 for 3, 5 for 2
nloutidx2=[9 9 9; 9 9 9; 9 9 9];
%nloutidx1=[3 3 3; 3 3 3; 2 2 2];   % alt sub 9 for 3, 5 for 2
%nloutidx2=[3 3 3; 3 3 3; 3 3 3];

rowcount=size(expidx1,1);
colcount=size(expidx2,2);
for ii=1:length(expidx1(:)),
   subplot(rowcount,colcount,ii);
   n1=[expstr{expidx1(ii)},'-',cnfstr{cnfidx(ii)},'-',instr{nlinidx1(ii)},...
       '-',outstr{nloutidx1(ii)}];
   n2=[expstr{expidx2(ii)},'-',cnfstr{cnfidx(ii)},'-',instr{nlinidx2(ii)},...
       '-',outstr{nloutidx2(ii)}];
   set1=predxc(:,expidx1(ii),nlinidx1(ii),cnfidx(ii),nloutidx1(ii));
   set2=predxc(:,expidx2(ii),nlinidx2(ii),cnfidx(ii),nloutidx2(ii));
   goodidx=find(~isnan(set1) & ~isnan(set2) & ...
             totlen(:,expidx1(ii),1)>MINLEN & totlen(:,expidx2(ii),1)>MINLEN);
   
   plotcomp(set1(goodidx),set2(goodidx),n1,n2);
end
%set(gcf,'PaperPosition',[0.25 0.25 8 10.5],'PaperOrientation','portrait');
set(gcf,'PaperPosition',[0.25 0.25 10.5 8],'PaperOrientation','landscape');
drawnow

% cross-class prediction
figure(3)
clf
expidx1=[1 1 1; 2 2 2; 1 1 1];
expidx2=[2 2 2; 3 3 3; 3 3 3];
nlinidx1=[3 3 3; 3 3 3; 3 3 3]; 
nlinidx2=[3 3 3; 3 3 3; 3 3 3];
cnfidx=[1 2 3; 1 2 3; 1 2 3];
nloutidx1=[9 9 9; 9 9 9; 9 9 9];
nloutidx2=[9 9 9; 9 9 9; 9 9 9];
%nloutidx1=[3 3 3; 3 3 3; 3 3 3];
%nloutidx2=[3 3 3; 3 3 3; 3 3 3];

rowcount=size(expidx1,1);
colcount=size(expidx2,2);
for ii=1:length(expidx1(:)),
   subplot(rowcount,colcount,ii);
   n1=[expstr{expidx1(ii)},'-',cnfstr{cnfidx(ii)},'-',instr{nlinidx1(ii)},...
       '-',outstr{nloutidx1(ii)}];
   n2=[expstr{expidx2(ii)},'-',cnfstr{cnfidx(ii)},'-',instr{nlinidx2(ii)},...
       '-',outstr{nloutidx2(ii)}];
   set1=predxc(:,expidx1(ii),nlinidx1(ii),cnfidx(ii),nloutidx1(ii));
   set2=predxc(:,expidx2(ii),nlinidx2(ii),cnfidx(ii),nloutidx2(ii));
   goodidx=find(~isnan(set1) & ~isnan(set2) & ...
             totlen(:,expidx1(ii),1)>MINLEN & totlen(:,expidx2(ii),1)>MINLEN);
   
   plotcomp(set1(goodidx),set2(goodidx),n1,n2);
end
%set(gcf,'PaperPosition',[0.25 0.25 8 10.5],'PaperOrientation','portrait');
set(gcf,'PaperPosition',[0.25 0.25 10.5 8],'PaperOrientation','landscape');


return


if 0,
   keyboard
end

cellid={preddata(ttt).cellid};
cellcount=length(cellid);
cellcount=cellcount+1;
cellid{cellcount}='mean';

runidx=find(okidx);
runidx(cellcount)=0;

s=size(preddata(end).predxc);
nlcount=s(4);
predxc=zeros(s(1),s(2)*(cellcount-1),s(3),s(4));
for ii=1:cellcount-1,
   nltemp1=min([size(preddata(ttt(ii)).predxc,1) s(1)]);
   nltemp4=min([size(preddata(ttt(ii)).predxc,4) s(4)]);
   
   predxc(1:nltemp1,ii,:,1:nltemp4)=...
       preddata(ttt(ii)).predxc(1:nltemp1,1,:,1:nltemp4);
end
predxc=cat(2,predxc,permute(nanmean(permute(predxc,[2 1 3 4])),[2 1 3 4]));

predcount=size(predxc,1);
nlcount=size(predxc,4);

for ii=1:cellcount,
   if 0,
      predrange=[1 2 3];
   elseif ismember(batchid,[30 31 32]),
      predrange=3;
   elseif ismember(batchid,[27 28 29]),
      predrange=2;
   else
      predrange=1;
   end
   for predidx=predrange,
      if predidx==predrange(1),
         fprintf('%-6s',cellid{ii});
      else
         fprintf('%-6s','');
      end
      for nlidx=1:nlcount,
         if ~isnan(predxc(predidx,ii,:,nlidx)),
            fprintf(' %7.3f',predxc(predidx,ii,:,nlidx));
         else
            fprintf(' %-7s','  x.xxx');
         end
      end
      fprintf('\n');
   end
end



