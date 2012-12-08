
% bad within class cells:
% (review) r0150b,r0213c
% (gr) r0148a, r0156a, r0158a, r0162b, r0168b, r0170a, r0208d, r0220a
% (nr) r0148a, r0156a, r0158a, r0162b, r0168b, r0169b, r0170a, r0223a
% bad gr but good rev pred: r0148a, r0156a, r0158a, r0162b, r0168b, 
%                           r0170a, r0208d, r0220a
% bad nr but good rev pred: r0148a, r0158a, r0162b, r0168b, r0170a


cellids={'r0148A','r0150b','r0156A','r0158A',...
         'r0162B','r0164C','r0166C','r0168B',...
         'r0170A','r0206B','r0208D','r0210A',...
         'r0211A','r0212B','r0215B','r0217B',...
         'r0219B','r0220A','r0221A','r0223A',...
         'r0225C'};

if 0
   %cellquality=[0 1 0 0 0 1 1 0 0 1 0 1 1 1 1 1 1 0 1 0 1];
   cellquality=[1 1 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 1];
   
   goodcells=find(cellquality);
   [xc,cid]=batchcomp([86;86],[2;3],1);
   [xc2,cid2]=batchcomp([86;86],[3;4],1);
   matchids=[];
   for ii=1:length(cid2),
      if sum(strcmp(cid2{ii},cellids)) > 0,
         matchids=[matchids ii];
      end
   end
   plotcomp(xc(matchids,1),xc(matchids,2),'gr','nr');
   plotcomp(xc2(matchids,1),xc2(matchids,2),'nr','rev');
end

hypredxc=[]; hyprederr=[];
for cellidx=1:length(cellids),
   [txc, terr]=xcpredtest(cellids{cellidx},86,1);
   hypredxc=[hypredxc txc];
   hyprederr=[hyprederr terr];
end
grpredxc=[]; grprederr=[];
for cellidx=1:length(cellids),
   [txc, terr]=xcpredtest(cellids{cellidx},87,2);
   grpredxc=[grpredxc txc];
   grprederr=[grprederr terr];
end
nrpredxc=[]; nrprederr=[];
for cellidx=1:length(cellids),
   [txc, terr]=xcpredtest(cellids{cellidx},88,2);
   nrpredxc=[nrpredxc txc];
   nrprederr=[nrprederr terr];
end
%poshypredxc=[]; poshyprederr=[];
%for cellidx=1:length(cellids),
%   [txc, terr]=xcpredtest(cellids{cellidx},83,1);
%   poshypredxc=[poshypredxc txc];
%   poshyprederr=[poshyprederr terr];
%end

cellquality=[2 1 2 2 2 1 1 2 2 1 2 1 1 1 1 1 1 2 1 2 1];

subplot(2,2,1);
set1=hypredxc(2,:)';
set2=hypredxc(3,:)';
e1=hyprederr(2,:)';
e2=hyprederr(3,:)';
set1hi=find(set1-e1 > set2+e2);
set2hi=find(set1+e1 < set2-e2);
sigidx=ones(size(set1))*2;
sigidx(set1hi)=3;
sigidx(set2hi)=1;
minx=min([-0.2 min(set1) min(set2)]);
plotcomp(set1,set2,'gr hybr','nr hybr',[minx 1 minx 1],sigidx);

subplot(2,2,2);
set1=hypredxc(3,:)';
set2=hypredxc(4,:)';
e1=hyprederr(3,:)';
e2=hyprederr(4,:)';
set1hi=find(set1-e1 > set2+e2);
set2hi=find(set1+e1 < set2-e2);
sigidx=ones(size(set1))*2;
sigidx(set1hi)=3;
sigidx(set2hi)=1;
minx=min([-0.2 min(set1) min(set2)]);
plotcomp(set1,set2,'nr hybr','rev',[minx 1 minx 1],sigidx);

subplot(2,2,3);
set1=grpredxc(3,:)';
set2=nrpredxc(4,:)';
e1=grprederr(3,:)';
e2=nrprederr(4,:)';
set1hi=find(set1-e1 > set2+e2);
set2hi=find(set1+e1 < set2-e2);
sigidx=ones(size(set1))*2;
sigidx(set1hi)=3;
sigidx(set2hi)=1;
minx=min([-0.2 min(set1) min(set2)]);
plotcomp(set1,set2,'gr','nr',[minx 1 minx 1],sigidx);


subplot(2,2,4);
bar(mean([grpredxc(3,:)' nrpredxc(4,:)' ...
          hypredxc(2,:)' hypredxc(3,:)' hypredxc(4,:)']));
ht=xticks(1:5,{'gr','nr','grhyb','nrhyb','rev'},...
          'FontSize',10,'Rotation',45,'HorizontalAlignment','right');
axis square
axis ([0 6 0 0.7]);

colormap(gray);

