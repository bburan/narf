function d=lsext(s)

NLCHAR=10;

r=ls(s);
ri=double(r);
edges=find(ri==10);

d=cell(length(edges),1);

lastedge=0;
for ii=1:length(edges),
   d{ii}=r(lastedge+1:edges(ii)-1);
   lastedge=edges(ii);
end
