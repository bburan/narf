% function [p,m,s]=randpairtest(x1,x2,n,tail,stat);
%
% n - number of randomizations (default 100)
% tail - 0 - two tailed (default)
%       -1 - test mean(x2-x1) < 0  ie, x2<x1
%        1 - test mean(x2-x1) > 0  ie, x2>x1
% stat - 'mean' (default) or 'median'
%
function [p,m,s]=randpairtest(x1,x2,n,tail,stat);

if length(x1)~=length(x2),
   disp('x1 and x2 must be same length!');
   return
end
if ~exist('stat','var'),
   stat='mean';
end
if strcmpi(stat,'median'),
   statnum=1;
else
   statnum=0;
end

if length(x1)==0,
   m=0;
   p=1;
   return
end
if length(x1)==1,
   m=x2-x1;
   p=1;
   return
end

if ~exist('n'),
   n=100;
end
if ~exist('tail'),
   tail=0;
end

dcount=length(x1);
mdiff=zeros(dcount,1);
if sum(abs(x2))==0,
   x2=x1;
   x1=zeros(size(x2));
end

for ii=1:n,
   bflip=round(rand(dcount,1));
   if statnum==0,
      mdiff(ii)=mean([x2(find(bflip)) -   x1(find(bflip)) ; ...
                      x1(find(1-bflip)) - x2(find(1-bflip))]);   
   else
      mdiff(ii)=median([x2(find(bflip)) -   x1(find(bflip)) ; ...
                        x1(find(1-bflip)) - x2(find(1-bflip))]);   
   end
end

if statnum==0,
   m=mean(x2-x1);
else
   m=median(x2-x1);
end
spos=round(n*0.66);
if tail==0,
   mlist=sort(-abs([mdiff; m]));
   bmin=max(find(mlist<=-abs(m)));
   s=abs(mlist(end-spos));
elseif tail==-1, % mean x2 < mean x1
   mlist=sort([mdiff; m]);
   bmin=max(find(mlist<=m));
   s=abs(mlist(end-spos));
elseif tail==1,
   mlist=sort(-[mdiff; m]);
   bmin=min(find(mlist>=-m));
   s=abs(mlist(spos));
end

s=std(mdiff);
p=bmin/(n+1);

%keyboard