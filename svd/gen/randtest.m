function p=randtest(x,m,s,n,diststr);

if ~exist('m'),
   m=0;
end
if ~exist('s'),
   s=1;
end
if ~exist('n'),
   n=100;
end
if ~exist('diststr'),
   diststr='gaussian';
end
dcount=length(x);

if ~strcmp(diststr,'uniform'),
   disp('randtest: only diststr==''uniform'' currently valid');
   return
end

if strcmp(diststr,'uniform'),
   a=m-s/2;
   b=m+s/2;
   
   mlist=zeros(n,1);
   for ii=1:n,
      y=rand(dcount,1)*s+a;
      mlist(ii)=mean(y);
   end
   
   mlist=sort(mlist);
   bmin=max([find(mlist<mean(x)); 0]);
   p=(bmin+1)/n;
end
