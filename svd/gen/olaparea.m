% function a=olaparea(dist,olap,r1,r2)
%
% support for venn.m
%
% created SVD 1/2/04
%
function a=olaparea(dist,olap,r1,r2)

d1=dist.*r1 ./(r1+r2);
d2=dist-d1;

h=sqrt(r1.^2-d1.^2);
theta=acos(d1./r1);

%h d1 d2]

a1=h*d1;
a2=h*d2;

if dist>(r1+r2),
   a=a1+a2;
else
   s1=theta * r1.^2;
   s2=theta * r2.^2;

   a=(s1+s2)-(a1+a2);
end

return

