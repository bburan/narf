% function [m,g,l,p]=randdiffs(x,y,N);
%
% cxy is the cross covariance of x and y.  then x is shuffled
% multiple times to get a measure of expected error cross
% covariance, exy, for 95% of shuffled trials.
%
% maxlag-number of bins +/- from 0 in cross-cov
% N-number of shuffles (default 100)
%
% created: SVD 10/17/04 (ripped off randxcov)
%
function [m,g,l,p]=randdiffs(x,y,N);

if ~exist('N','var'),
   N=100;
end

gidx=find(~isnan(x+y));
x=x(gidx);
y=y(gidx);

L=length(x);
epsilon=1e-100;

mx=mean(x);
my=mean(y);
m=mean(y)-mean(x);
g=0;
p=ones(1,3);

if sum(abs(x-mx))<=epsilon | sum(abs(y-my))<=epsilon,
   return
end


mr=(x+y)./2;
dr=(y-x)./2;

bootcount=20;
m=zeros(bootcount,1);
g=zeros(bootcount,1);

yest=zeros(L,1);

for bootidx=1:bootcount,
   
   bidx=[1:round((bootidx-1)/bootcount*L) (round((bootidx)/bootcount*L)+1):L];
   vidx=[(round((bootidx-1)/bootcount*L)+1):round((bootidx)/bootcount*L)];
   
   mpg=polyfit(mr(bidx),dr(bidx),1);
   m(bootidx)=mpg(2);
   g(bootidx)=mpg(1);
   
   yest(vidx)=(1+g(bootidx))./(1-g(bootidx)) .* x(vidx) + 2 * m(bootidx);
   
end

me=std(m).*sqrt(bootcount);
m=mean(m);
ge=std(g).*sqrt(bootcount);
g=mean(g);

mmax=max([yest(:);y(:)]);
mm=linspace(0,mmax,100);

outidx=find(y>yest+sqrt(yest).*2 | y<yest-sqrt(yest).*2);

plot(yest,y,'.');
hold on
plot(yest(outidx),y(outidx),'r.');

plot(mm,mm,'k-');
%tt=repmat(mm,[200 1]);
%pp=sort(random('poiss',tt));
%plot(mm,pp(round(size(pp,1).*0.05),:),'k--');
%plot(mm,pp(round(size(pp,1).*0.95),:),'k--');
plot(mm,mm+sqrt(mm).*2,'k--');
plot(mm,mm-sqrt(mm).*2,'k--');
hold off

axis equal
axis square


l= length(noiseidx)./L





xrange=0:100;
murange=0:50;

poissprob=zeros(length(xrange),length(murange));
for ii=1:length(xrange),
   poissprob(ii,:)=murange.^xrange(ii) .* exp(-murange) ./ factorial(xrange(ii));
end

psum=sum(poissprob,2);
invpoiss=poissprob./repmat(psum,[1 length(murange)]);

psame=zeros(size(y));
yest(find(yest<0))=0;
for ii=1:L,
   psame(ii)=sum(min([invpoiss(round(y(ii)+1),:);invpoiss(round(yest(ii)+1),:)]));
end  




[xx,lambda]=meshgrid(0:50,0:50);

pxgivenl=lambda.^xx .* exp(-lambda);
psum=sum(pxgivenl,2);
pxgivenl=pxgivenl./repmat(psum,[1 size(xx,2)]);

psum=sum(pxgivenl,1);
plgivenx=pxgivenl./repmat(psum,[size(xx,1) 1]);






% true defs:
% x= mr-dr
% y= mr+dr
%
% fit dr= mr*g + m
%
% model for x:  x= mr - (m+mr*g)
% model for y:  y= mr + (m+mr*g)

yest=(1+g)./(1-g) .* x + 2 * m;

mmax=max([yest(:);y(:)]);
mm=linspace(0,mmax,100);

plot(yest,y,'.');
hold on
plot([0 max([yest(:);y(:)])],[0 max([yest(:);y(:)])],'k-');
plot(mm,mm+sqrt(mm).*3,'k--');
plot(mm,mm-sqrt(mm).*3,'k--');
hold off


mm=linspace(0,max(mr),100);

outidx=find(abs(dr)>sqrt(mr)*2);

plot(y-2*(m+mr*g),x,'.');
hold on
plot(mm,mm-(m+mm*g)+sqrt(mm-(m+mm*g))*3,'k--');
plot(mm,mm-(m+mm*g)-sqrt(mm-(m+mm*g))*3,'k--');
hold off

if 0

plot(mr,dr,'.');
hold on
plot([min(mr) max(mr)],[min(mr) max(mr)].*mpg(1)+mpg(2),'k-');
%plot(mm,(m+mm*g)+sqrt(mm+(m+mm*g))*3,'k--');
%plot(mm,(m+mm*g)-sqrt(mm-(m+mm*g))*3,'k--');
plot(mm, sqrt(mm + (m+mm*g))*2,'k--');
plot(mm,-sqrt(mm - (m+mm*g))*2,'k--');
%plot(mr(outidx),dr(outidx),'r.');
hold off

fracount=length(outidx)./L
end

if 0

nxhi=find(x>mr+(m+mr*g)+sqrt(mr+(m+mr*g))*3)
nyhi=find(y>mr-(m+mr*g)+sqrt(mr-(m+mr*g))*3)

plot(x,y,'.');
hold on
plot(mm+(m+mm*g)+sqrt(mm+(m+mm*g))*3,mm,'k--');
plot(mm,mm-(m+mm*g)+sqrt(mm-(m+mm*g))*3,'k--');
plot(x([nxhi;nyhi]),y([nxhi;nyhi]),'r.');
hold off
axis equal;
axis square

end


l=0;



return

% test mean diffs
rm=zeros(N,1);
for ii=1:N,
   switchidx=(rand(L,1)>0.5);
   rx=[x(switchidx==0);y(switchidx==1)];
   ry=[y(switchidx==0);x(switchidx==1)];
   rm(ii)=mean(ry)-mean(rx);
end

rm=sort(abs([rm;m]));
p(1)=sum(abs(m)<=rm)./(N+1);

if m<0,
   x1=x;
   y1=y-m;
else
   x1=x+m;
   y1=y;
end

g=(x1'*y1)/(x1'*x1);

rg=zeros(N,1);
for ii=1:N,
   switchidx=(rand(L,1)>0.5);
   rx=[x1(switchidx==0);y1(switchidx==1)];
   ry=[y1(switchidx==0);x1(switchidx==1)];
   rg(ii)=(rx'*ry)/(rx'*rx);
end
   
hist(rg);

rg=sort(abs(log([rg;g])));
p(2)=sum(abs(log(g))<=rg)./(N+1);

x2=x1.*g;
y2=y1;

if 1,
   %hist(rg);
else
   
   mr=sortrows([x1+y1 x1 y1 x2 y2]);
   plot(gsmooth(mr(:,2:3)),'--');
   hold on
   plot(gsmooth(mr(:,4:5)));
   hold off
   mean([x2 y2])
end

l=xcorr(x2,y2,0,'coeff');

rl=zeros(N,1);
for ii=1:N,
   switchidx=(rand(L,1)>0.5);
   rx=[x2(switchidx==0);y2(switchidx==1)];
   ry=[y2(switchidx==0);x2(switchidx==1)];
   rl(ii)=xcorr(rx,ry,0,'coeff');
end

if 0
hist(rl);
hold on
plot(l,N/10,'x');
hold off
end

rl=sort([rl;l]);
p(3)=sum(l>=rl)./(N+1);





