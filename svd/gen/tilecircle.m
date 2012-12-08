% function im=tilecircle(r,h);
%
% parameters:
% r=radius of circle
% h=height of triangles
%
% returns:
% im=image with each 45-deg triangle inside the circle tiled with a different number
%
% Simplest if r and h are both even numbers and that r be an integer
% multiple of h
%
function im=tilecircle(r,h);

d=r*2;
im=zeros(d,d);

w=h;

[xx,yy]=meshgrid(1:d,1:d);

cim=zeros(d,d);
cim(sqrt((xx-r-0.5).^2 + (yy-r-0.5).^2)<r)=1;

rows=ceil(d./h);
trianglesperrow=ceil(d./w);

tridx=0;
counter=0;
for rr=1:rows,
    r0=(rr-1).*h+1;
    r1=rr.*h;
    for tt=1:trianglesperrow,
        
        c0=(tt-1).*w+1;
        c1=tt.*w;
        
        if ((tt+rr)./2)~=round((tt+rr)./2),
            b0=tt.*w;
            s0=-1;
            ff1=find(xx>=r0 & xx<=r1 & yy>=c0 & yy<=c1 & yy-b0<=s0*(xx-r0)-double(xx-r0>=h./2));
            ff2=find(xx>=r0 & xx<=r1 & yy>=c0 & yy<=c1 & yy-b0>s0*(xx-r0)-double(xx-r0>=h./2));
        else
            b0=(tt-1).*w;
            s0=1;
            ff1=find(xx>=r0 & xx<=r1 & yy>=c0 & yy<=c1 & yy-b0<=s0*(xx-r0)+1*double(xx-r0<h./2));
            ff2=find(xx>=r0 & xx<=r1 & yy>=c0 & yy<=c1 & yy-b0>s0*(xx-r0)+1*double(xx-r0<h./2));
        end
        
        if sum(cim(ff1)==0)==0,
            counter=counter+1;
            im(ff1)=counter;
        end
        if sum(cim(ff2)==0)==0,
            counter=counter+1;
            im(ff2)=counter;
        end
    end
end

figure;
imagesc(im);