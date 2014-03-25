% function [p,m,s]=randttest(x1,x2,n,tail);
%
% n - number of randomizations (default 100)
% tail - 0 - two tailed (default)
%       -1 - test mean(x2) < mean(x1)  ie, x2<x1
%        1 - test mean(x2) > mean(x1)  ie, x2>x1
%
function [p,m,s]=randttest(x1,x2,n,tail);

if length(x1)==0 | length(x2)==0,
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

xall=[x1(:);x2(:)];
dcount1=length(x1(:));
dcount2=length(x2(:));

mdiff=zeros(n,1);

for ii=1:n,
   xshuff=shuffle(xall);
   
   mdiff(ii)=mean(xshuff(dcount1+1:end))-mean(xshuff(1:dcount1));
end

m=mean(x2(:))-mean(x1(:));
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