% function [xc,celllist,xcerr]=batchcomp(batches,nloutidx,cnfidx,xcflag,celllist);
%
% xc(:,1) is first set of pred scores
% xc(:,2) is second
% celllist is a cellarray, where each entry gives the cell whose
% data is in that corresponding row of xc
% 
% xcflag: 0 (default) standard cc computed over time (predxc)
%         1 cc over fixation-averaged responses (predfix)
%         2 cc corrected for validataion noise (predinf)
%
% created SVD 12/02
%
function [xc,celllist,xcerr]=multibatchcomp(batches,nloutidx,cnfidx,xcflag,celllist);

% callmaster : exp X in-model X cell
% predxc: exp X in-model X cell X out-model X cnf

dbopen;

batchcount=length(batches);

bname={};
for ii=1:length(batches),
   bname{ii}=dbget('sBatch',batches(ii),'name');
end
if ~exist('nloutidx','var'),
   nloutidx=ones(size(batches)).*2;
elseif length(nloutidx)==1,
   nloutidx=ones(size(batches)).*nloutidx;
end
if ~exist('cnfidx','var'),
   cnfidx=ones(size(batches));
elseif length(cnfidx)==1,
   cnfidx=ones(size(batches)).*cnfidx;
end
if ~exist('celllist','var'),
   usemasterlist=0;
else
   usemasterlist=1;
end
if ~exist('newfig','var'),
   newfig=1;
end
if ~exist('xcflag','var'),
   xcflag=0;
end
if ~exist('bcellid','var'),
   bcellid=0;
end

predxcB={};
cellidB={};
prederrB={};
for batidx=1:batchcount,
   [predxcB{batidx},cellidB{batidx},prederrB{batidx}]= ...
       dbresv1(batches(batidx),xcflag);
   if ~usemasterlist,
      if batidx==1,
         celllist=cellidB{batidx};
      else
         celllist=intersect(celllist,cellidB{batidx});
      end
   end
end

cellcount=length(celllist);
if ~cellcount,
   disp('no matches across batches!');
   xc=[];
   return;
end

xc=zeros(cellcount,batchcount);
xcerr=zeros(cellcount,batchcount);
for batidx=1:batchcount,
   [c,ia,ib]=intersect(cellidB{batidx},celllist);
   xc(ib,batidx)=predxcB{batidx}(cnfidx(batidx),nloutidx(batidx),1,ia);
   xcerr(ib,batidx)=prederrB{batidx}(cnfidx(batidx),nloutidx(batidx),1,ia);
end

return




% load pred results for each batch
res=loadbigdata(batches,nloutidx,cnfidx);
if xcflag==1,
   predxc=res.predfix;
   prederr=res.predfixerr;
elseif xcflag==2,
   predxc=res.predinf;
   prederr=res.prederr;
elseif xcflag==0 || xcflag==3,
   predxc=res.predxc;
   prederr=res.prederr;
elseif xcflag==4
   predxc=res.expxc;
   prederr=res.prederr;
end
predp=res.predp;
celllist=res.celllist;

%badcellidx=find(strcmp(celllist,'r0150B') | strcmp(celllist,'93G83A') ...
%    | strcmp(celllist,'models'));
%predxc(badcellidx,:)=nan;

if 0,
MINLEN=500;
badcellidx=find(res.totlen(:,1)<MINLEN);
predxc(badcellidx,1,:,:,:)=nan;
end

PMIN=0.01;
expstr={'Rev','GR','NR'};
instr={'pix','pow','psf'};
cnfstr=expstr;
outstr={'lin','thr','ethr','full'};

% review-only model comparison
if newfig,
   figure;
end

n1=sprintf('%s (onl=%d cidx=%d)',bname{1},nloutidx(1),cnfidx(1));
n2=sprintf('%s (onl=%d cidx=%d)',bname{2},nloutidx(2),cnfidx(2));

set1=predxc(:,1);
set2=predxc(:,2);
goodidx=find(~isnan(set1) & ~isnan(set2));
set1=set1(goodidx);
set2=set2(goodidx);
p1=predp(goodidx,1);
p2=predp(goodidx,2);
e1=prederr(goodidx,1);
e2=prederr(goodidx,2);

set1hi=find(set1-e1 > set2+e2);
set2hi=find(set1+e1 < set2-e2);
sigidx=ones(size(set1))*2;
sigidx(set1hi)=3;
sigidx(set2hi)=1;

%sigidx=2-(set1(goodidx)+e1(goodidx) < set2(goodidx)-e2(goodidx));
%sigidx=2-(p1<PMIN | p2<PMIN);

minx=min([-0.2 min(set1) min(set2)]);

[p,m]=plotcomp(set1,set2,n1,n2,[minx 1.0 minx 1.0],sigidx,'mean');

if bcellid,
   hold on
   for ii=1:length(goodidx),
      text(set1(ii)-0.01,set2(ii)+0.01,celllist{goodidx(ii)},...
           'FontSize',8,'HorizontalAlignment','right');
   end
   hold off
end

if nargout>0,
   xc=[set1 set2];
else
   sres={'<',' ','>'};
   for ii=1:length(goodidx),
      
      fprintf('%-7s %6.3f  %s  %6.3f\n',celllist{goodidx(ii)},...
              set1(ii),sres{sigidx(ii)},set2(ii));
   end
end

celllist={celllist{goodidx}};

fprintf('%s v %s:\n',n1,n2);
fprintf('mean  : %.2f v %.2f (p<%.3f) %d > %d < %d / %d \n',...
        mean(set1),mean(set2),p(end),sum(sigidx==3),...
        sum(sigidx==2),sum(sigidx==1),length(sigidx))
fprintf('l2norm: %.2f v %.2f \n',sqrt(mean(set1.^2)),sqrt(mean(set2.^2)))

