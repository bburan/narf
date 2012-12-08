% function zout=edjoinsmooth(z1,z2,curvecut,n);
function zout=edjoinsmooth(z1,z2,curvecut,n);

if ~exist('curvecut','var'),
   curvecut=0.3;
end
if ~exist('n','var'),
   n=21;
end

ppidx=min(find(abs(z1-z1(end))<curvecut));
if ~isempty(ppidx),
   z1=z1(1:(ppidx-1));
end
ppidx=max(find(abs(z2-z2(1))<curvecut));
if ~isempty(ppidx),
   z2=z2((ppidx+1):end);
end


if abs(z1(end)-z2(1))<1e-10,
   zout=[z1;z2];
   return;
end

t1=z1(end)-z1(end-1);
t2=z2(2)-z2(1);

t1=t1./norm(t1);
t2=t2./norm(t2);

p1=exp(i*pi/2)*t1;
p2=exp(-i*pi/2)*t2;

t=(linspace(0,1,n))';
zjoin=(2.*t.^3-3.*t.^2+1).*z1(end)+(-2.*t.^3+3.*t.^2).*z2(1) + ...
      (t.^3 - 2.*t.^2 + t).*t1 + (t.^3-t.^2).*t2;

if 0,
   figure(3);
   plot([z1(end) z2(1)],'o');
   hold on
   plot([z1(end) z1(end)+t1]);
   plot([z2(1) z2(1)+t2]);
   plot(z1,'g');
   plot(z2,'g');
   plot(zjoin,'r');
   hold off
end

zout=[z1;zjoin;z2];

