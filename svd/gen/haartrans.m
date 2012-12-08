function wc=haartrans(f,minj,maxj);

n=length(f);
J=floor(log2(n));

if maxj>J,
   maxj=J;
end
if minj<1,
   minj=1;
end

wc=zeros(n,maxj-minj+1);

for j=minj:maxj,
   g=haar(0:2^j-1,0,j);
   rt=conv(f,flipud(g));
   %keyboard
   %[sum(rt)*1000 sum(rt(length(rt)-length(g)+1:length(rt))) sum(rt(1:length(g)-1)) length(g) length(rt) length(rt)-length(g)*2+1]
   %rt(1:length(g)-1)=0;
   rt(length(rt)-length(g)+1:length(rt))=0;
   %wc(:,j-minj+1)=rt(1:length(rt)-length(g)+1);
   wc(:,j-minj+1)=rt(length(g):length(rt));
   %sum(wc(:,j-minj+1))
end

   

