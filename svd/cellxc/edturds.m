function edturds

global SUBIDX
SUBIDX=[];

brange=[90 135 180 225 270];
N=41;

% number of convexities, cidx=1 most convex, cidx=4 most concave
cvcount=4;

xset=zeros(N,cvcount,length(brange));
yset=zeros(N,cvcount,length(brange));
erotset=ones(size(xset))*i;
srotset=ones(size(xset))*i;
brotset=ones(size(xset))*i;

figure(1);
for bidx=1:length(brange);
   b=brange(bidx);
   subplot(1,length(brange),bidx);
   
   if 0,
      a2=cos(b*pi/180./2);
      b2=sin(b*pi/180./2);
      a1=a2;
      b1=b2;
   else
      a2=cos(b*pi/180   );
      b2=sin(b*pi/180   );
      a1=1;
      b1=0;
   end
   
   if b<=180,
      [xset(:,1,bidx),yset(:,1,bidx)]=arc(0,0,0,brange(bidx),1,N);
      [xset(:,2,bidx),yset(:,2,bidx)]=arc12(a1, b1,a2,b2, 2,N);
      [xset(:,3,bidx),yset(:,3,bidx)]=arc12(a1,b1,a2, b2,-2,N);
      [xset(:,4,bidx),yset(:,4,bidx)]=arc12(a1,b1,a2, b2,-1,N);
   else
      [xset(:,1,bidx),yset(:,1,bidx)]=arc(0,0,0,brange(bidx),1,N);
      [xset(:,2,bidx),yset(:,2,bidx)]=arc12(a1, b1,a2,b2, 1,N);
      [xset(:,3,bidx),yset(:,3,bidx)]=arc12(a1,b1,a2, b2, 2,N);
      [xset(:,4,bidx),yset(:,4,bidx)]=arc12(a1,b1,a2, b2,-2,N);
   end
   
   cla
   plot(xset(:,1,bidx),yset(:,1,bidx));
   hold on
   plot(xset(:,2,bidx),yset(:,2,bidx),'r');
   plot(xset(:,3,bidx),yset(:,3,bidx),'g');
   plot(xset(:,4,bidx),yset(:,4,bidx),'c');
   hold off
   axis([-1 1 -1 1]);
   axis square
   title(sprintf('ang=%.0f',b));
   
   erotset(:,:,bidx)=exp(-i*pi/36)*(xset(:,:,bidx)-1+i*yset(:,:,bidx))+1;
   srotset(:,:,bidx)=exp(-i*b*pi/180)*(xset(:,:,bidx)+i*yset(:,:,bidx));
   srotset(:,:,bidx)=exp(i*pi/36)*(srotset(:,:,bidx)-1)+1;
   srotset(:,:,bidx)=exp(i*b*pi/180)*srotset(:,:,bidx);
   
   brotset(:,:,bidx)=exp(-i*b*pi/180)*erotset(:,:,bidx);
   brotset(:,:,bidx)=exp(i*pi/36)*(brotset(:,:,bidx)-1)+1;
   brotset(:,:,bidx)=exp(i*b*pi/180)*brotset(:,:,bidx);
end

nrotset=xset+i.*yset;

figure(2);
clf
fullpage portrait

%
% assemble 2-projection stimuli
%
%brange=[90 135 180 225 270];

%sep=90: use bidx=1, 5
rotamt=[90] *pi/180;
b1=1; cv1=3;
b2=5; cv2=2;

z1=xset(:,cv1,b1) + i * yset(:,cv1,b1);
z2=exp(i.*rotamt) * (xset(:,cv2,b2) + i * yset(:,cv2,b2));
plotnext([z1;z2],8);

zout=[srotset(:,cv1,b1);exp(i.*rotamt) * erotset(:,cv2,b2)];
m=round(length(zout)/2);
zout=edjoinsmooth(zout(m:end),zout(1:(m-1)));
plotnext(zout,8);

zout=edjoinsmooth(erotset(:,cv1,b1),exp(i.*rotamt) * srotset(:,cv2,b2));
plotnext(zout,8);

zout=edjoinsmooth(brotset(:,cv1,b1),exp(i.*rotamt) * brotset(:,cv2,b2));
m=round(length(zout)/2);
zout=edjoinsmooth(zout(m:end),zout(1:(m-1)));
plotnext(zout,8);

rotamt=[90]*pi/180;
b1=1; cv1=4;
b2=5; cv2=1;
z1=xset(:,cv1,b1) + i * yset(:,cv1,b1);
z2=exp(i.*rotamt) * (xset(:,cv2,b2) + i * yset(:,cv2,b2));
%plotnext([z1;z2]);

zout=edjoinsmooth(brotset(:,cv1,b1),exp(i.*rotamt) * brotset(:,cv2,b2),0.4);
m=round(length(zout)/2);
zout=edjoinsmooth(zout(m:end),zout(1:(m-1)));
plotnext(zout,8);

%sep=135: use bidx=2,4 
rotamt=[135] *pi/180;
b1=2; cv1=3;
b2=4; cv2=2;
z1=xset(:,cv1,b1) + i * yset(:,cv1,b1);
z2=exp(i.*rotamt) * (xset(:,cv2,b2) + i * yset(:,cv2,b2));
plotnext([z1;z2],8);

zout=[srotset(:,cv1,b1);exp(i.*rotamt) * erotset(:,cv2,b2)];
m=round(length(zout)/2);
zout=edjoinsmooth(zout(m:end),zout(1:(m-1)));
plotnext(zout,8);

zout=edjoinsmooth(erotset(:,cv1,b1),exp(i.*rotamt) * srotset(:,cv2,b2));
plotnext(zout,8);

zout=edjoinsmooth(z1,z2);
m=round(length(zout)/2);
zout=edjoinsmooth(zout(m:end),zout(1:(m-1)));
plotnext(zout,8);

%sep=180: use bidx=3,3
rotamt=[180] *pi/180;
b1=3; cv1=2;
b2=3; cv2=2;
z1=xset(:,cv1,b1) + i * yset(:,cv1,b1);
z2=exp(i.*rotamt) * (xset(:,cv2,b2) + i * yset(:,cv2,b2));
plotnext([z1;z2],4);

zout=[srotset(:,cv1,b1);exp(i.*rotamt) * erotset(:,cv2,b2)];
m=round(length(zout)/2);
zout=edjoinsmooth(zout(m:end),zout(1:(m-1)));
plotnext(zout,8);

%zout=edjoinsmooth(erotset(:,cv1,b1),exp(i.*rotamt) * srotset(:,cv2,b2));
%plotnext(zout);

zout=edjoinsmooth(brotset(:,cv1,b1),exp(i.*rotamt) * brotset(:,cv2,b2));
m=round(length(zout)/2);
zout=edjoinsmooth(zout(m:end),zout(1:(m-1)));
plotnext(zout,4);


%
% assemble 3-projection stimuli
%
%brange=[90 135 180 225 270];

%sep=90,180: use bidx=1,1,3 
rotamt=[90 180] *pi/180;
b1=1; cv1=4;
b2=1; cv2=4;
b3=3; cv3=2;
z1=xset(:,cv1,b1) + i * yset(:,cv1,b1);
z2=exp(i.*rotamt(1)) * (xset(:,cv2,b2) + i * yset(:,cv2,b2));
z3=exp(i.*rotamt(2)) * (xset(:,cv3,b3) + i * yset(:,cv3,b3));
plotnext([z1;z2;z3],8);

zout=edjoinsmooth(srotset(:,cv1,b1),exp(i.*rotamt(1)) * nrotset(:,cv2,b2),0.5);
zout=edjoinsmooth(zout,exp(i.*rotamt(2)) * erotset(:,cv3,b3),0.5);
m=round(length(zout)/2);
zout=edjoinsmooth(zout(m:end),zout(1:(m-1)),0.5);
plotnext(zout,8);

zout=edjoinsmooth(erotset(:,cv1,b1),exp(i.*rotamt(1)) * srotset(:,cv2,b2),0.5);
zout=edjoinsmooth(zout,exp(i.*rotamt(2)) * nrotset(:,cv3,b3),0.5);
m=round(length(zout)/2);
zout=edjoinsmooth(zout(m:end),zout(1:(m-1)),0.5);
plotnext(zout,8);

zout=edjoinsmooth(brotset(:,cv1,b1),exp(i.*rotamt(1)) * srotset(:,cv2,b2),0.5);
zout=edjoinsmooth(zout,exp(i.*rotamt(2)) * erotset(:,cv3,b3),0.5);
m=round(length(zout)/2);
zout=edjoinsmooth(zout(m:end),zout(1:(m-1)),0.5);
plotnext(zout,8);

zout=edjoinsmooth(nrotset(:,cv1,b1),exp(i.*rotamt(1)) * erotset(:,cv2,b2),0.5);
zout=edjoinsmooth(zout,exp(i.*rotamt(2)) * srotset(:,cv3,b3),0.5);
m=round(length(zout)/2);
zout=edjoinsmooth(zout(m:end),zout(1:(m-1)),0.5);
plotnext(zout,8);

zout=edjoinsmooth(srotset(:,cv1,b1),exp(i.*rotamt(1)) * erotset(:,cv2,b2),0.5);
zout=edjoinsmooth(zout,exp(i.*rotamt(2)) * brotset(:,cv3,b3),0.5);
m=round(length(zout)/2);
zout=edjoinsmooth(zout(m:end),zout(1:(m-1)),0.5);
plotnext(zout,8);

zout=edjoinsmooth(erotset(:,cv1,b1),exp(i.*rotamt(1)) * brotset(:,cv2,b2),0.5);
zout=edjoinsmooth(zout,exp(i.*rotamt(2)) * srotset(:,cv3,b3),0.5);
m=round(length(zout)/2);
zout=edjoinsmooth(zout(m:end),zout(1:(m-1)),0.5);
plotnext(zout,8);

zout=edjoinsmooth(brotset(:,cv1,b1),exp(i.*rotamt(1)) * brotset(:,cv2,b2),0.5);
zout=edjoinsmooth(zout,exp(i.*rotamt(2)) * brotset(:,cv3,b3),0.5);
m=round(length(zout)/2);
zout=edjoinsmooth(zout(m:end),zout(1:(m-1)),0.5);
plotnext(zout,8);

%
%sep=90,180: use bidx=1,1,3 
%
rotamt=[90 180] *pi/180;
b1=1; cv1=1;
b2=1; cv2=4;
b3=3; cv3=2;
z1=xset(:,cv1,b1) + i * yset(:,cv1,b1);
z2=exp(i.*rotamt(1)) * (xset(:,cv2,b2) + i * yset(:,cv2,b2));
z3=exp(i.*rotamt(2)) * (xset(:,cv3,b3) + i * yset(:,cv3,b3));
%plotnext([z1;z2;z3]);

zout=edjoinsmooth(brotset(:,cv1,b1),exp(i.*rotamt(1)) * srotset(:,cv2,b2));
zout=edjoinsmooth(zout,exp(i.*rotamt(2)) * erotset(:,cv3,b3));
m=round(length(zout)/2);
zout=edjoinsmooth(zout(m:end),zout(1:(m-1)));
plotnext(zout,8);

zout=edjoinsmooth(brotset(:,cv1,b1),exp(i.*rotamt(1)) * brotset(:,cv2,b2));
zout=edjoinsmooth(zout,exp(i.*rotamt(2)) * brotset(:,cv3,b3));
m=round(length(zout)/2);
zout=edjoinsmooth(zout(m:end),zout(1:(m-1)));
plotnext(zout,8);

%
%sep=90,180: use bidx=1,1,3 
%
rotamt=[90 180] *pi/180;
b1=1; cv1=4;
b2=1; cv2=1;
b3=3; cv3=2;
z1=xset(:,cv1,b1) + i * yset(:,cv1,b1);
z2=exp(i.*rotamt(1)) * (xset(:,cv2,b2) + i * yset(:,cv2,b2));
z3=exp(i.*rotamt(2)) * (xset(:,cv3,b3) + i * yset(:,cv3,b3));
%plotnext([z1;z2;z3]);

zout=edjoinsmooth(erotset(:,cv1,b1),exp(i.*rotamt(1)) * brotset(:,cv2,b2));
zout=edjoinsmooth(zout,exp(i.*rotamt(2)) * srotset(:,cv3,b3));
m=round(length(zout)/2);
zout=edjoinsmooth(zout(m:end),zout(1:(m-1)));
plotnext(zout,8);

zout=edjoinsmooth(brotset(:,cv1,b1),exp(i.*rotamt(1)) * brotset(:,cv2,b2));
zout=edjoinsmooth(zout,exp(i.*rotamt(2)) * brotset(:,cv3,b3));
m=round(length(zout)/2);
zout=edjoinsmooth(zout(m:end),zout(1:(m-1)));
plotnext(zout,8);

rotamt=[90 180] *pi/180;
b1=1; cv1=4;
b2=1; cv2=4;
b3=3; cv3=1;
z1=xset(:,cv1,b1) + i * yset(:,cv1,b1);
z2=exp(i.*rotamt(1)) * (xset(:,cv2,b2) + i * yset(:,cv2,b2));
z3=exp(i.*rotamt(2)) * (xset(:,cv3,b3) + i * yset(:,cv3,b3));
%plotnext([z1;z2;z3]);

zout=edjoinsmooth(srotset(:,cv1,b1),exp(i.*rotamt(1)) * erotset(:,cv2,b2),0.5);
zout=edjoinsmooth(zout,exp(i.*rotamt(2)) * brotset(:,cv3,b3));
m=round(length(zout)/2);
zout=edjoinsmooth(zout(m:end),zout(1:(m-1)));
plotnext(zout,8);

zout=edjoinsmooth(brotset(:,cv1,b1),exp(i.*rotamt(1)) * brotset(:,cv2,b2),0.5);
zout=edjoinsmooth(zout,exp(i.*rotamt(2)) * brotset(:,cv3,b3));
m=round(length(zout)/2);
zout=edjoinsmooth(zout(m:end),zout(1:(m-1)));
plotnext(zout,8);

if 0
   zpre=[brotset(:,cv1,b1);
         exp(i.*rotamt(1)) * brotset(:,cv2,b2);
         exp(i.*rotamt(2)) * brotset(:,cv3,b3)];
   hold on
   plot(zpre,'g');
   hold off
end

%
%sep=90,135: use bidx=1,2,2 
%
rotamt=[90 135+90] *pi/180;
b1=1; cv1=4;
b2=2; cv2=3;
b3=2; cv3=3;
z1=xset(:,cv1,b1) + i * yset(:,cv1,b1);
z2=exp(i.*rotamt(1)) * (xset(:,cv2,b2) + i * yset(:,cv2,b2));
z3=exp(i.*rotamt(2)) * (xset(:,cv3,b3) + i * yset(:,cv3,b3));
plotnext([z1;z2;z3],8);

zout=edjoinsmooth(srotset(:,cv1,b1),exp(i.*rotamt(1)) * nrotset(:,cv2,b2),0.5);
zout=edjoinsmooth(zout,exp(i.*rotamt(2)) * erotset(:,cv3,b3),0.5);
m=round(length(zout)/2);
zout=edjoinsmooth(zout(m:end),zout(1:(m-1)),0.5);
plotnext(zout,8);

zout=edjoinsmooth(erotset(:,cv1,b1),exp(i.*rotamt(1)) * srotset(:,cv2,b2),0.5);
zout=edjoinsmooth(zout,exp(i.*rotamt(2)) * nrotset(:,cv3,b3),0.5);
m=round(length(zout)/2);
zout=edjoinsmooth(zout(m:end),zout(1:(m-1)),0.5);
plotnext(zout,8);

zout=edjoinsmooth(brotset(:,cv1,b1),exp(i.*rotamt(1)) * srotset(:,cv2,b2),0.5);
zout=edjoinsmooth(zout,exp(i.*rotamt(2)) * erotset(:,cv3,b3),0.5);
m=round(length(zout)/2);
zout=edjoinsmooth(zout(m:end),zout(1:(m-1)),0.5);
plotnext(zout,8);

zout=edjoinsmooth(nrotset(:,cv1,b1),exp(i.*rotamt(1)) * erotset(:,cv2,b2),0.5);
zout=edjoinsmooth(zout,exp(i.*rotamt(2)) * srotset(:,cv3,b3),0.5);
m=round(length(zout)/2);
zout=edjoinsmooth(zout(m:end),zout(1:(m-1)),0.5);
plotnext(zout,8);

zout=edjoinsmooth(srotset(:,cv1,b1),exp(i.*rotamt(1)) * erotset(:,cv2,b2),0.5);
zout=edjoinsmooth(zout,exp(i.*rotamt(2)) * brotset(:,cv3,b3),0.5);
m=round(length(zout)/2);
zout=edjoinsmooth(zout(m:end),zout(1:(m-1)),0.5);
plotnext(zout,8);

zout=edjoinsmooth(erotset(:,cv1,b1),exp(i.*rotamt(1)) * brotset(:,cv2,b2),0.5);
zout=edjoinsmooth(zout,exp(i.*rotamt(2)) * srotset(:,cv3,b3),0.5);
m=round(length(zout)/2);
zout=edjoinsmooth(zout(m:end),zout(1:(m-1)),0.5);
plotnext(zout,8);

zout=edjoinsmooth(brotset(:,cv1,b1),exp(i.*rotamt(1)) * brotset(:,cv2,b2),0.5);
zout=edjoinsmooth(zout,exp(i.*rotamt(2)) * brotset(:,cv3,b3),0.5);
m=round(length(zout)/2);
zout=edjoinsmooth(zout(m:end),zout(1:(m-1)),0.5);
plotnext(zout,8);

rotamt=[90 135+90] *pi/180;
b1=1; cv1=2;
b2=2; cv2=3;
b3=2; cv3=3;
z1=xset(:,cv1,b1) + i * yset(:,cv1,b1);
z2=exp(i.*rotamt(1)) * (xset(:,cv2,b2) + i * yset(:,cv2,b2));
z3=exp(i.*rotamt(2)) * (xset(:,cv3,b3) + i * yset(:,cv3,b3));
%plotnext([z1;z2;z3]);

zout=edjoinsmooth(brotset(:,cv1,b1),exp(i.*rotamt(1)) * srotset(:,cv2,b2));
zout=edjoinsmooth(zout,exp(i.*rotamt(2)) * erotset(:,cv3,b3));
m=round(length(zout)/2);
zout=edjoinsmooth(zout(m:end),zout(1:(m-1)));
plotnext(zout,8);

zout=edjoinsmooth(brotset(:,cv1,b1),exp(i.*rotamt(1)) * brotset(:,cv2,b2));
zout=edjoinsmooth(zout,exp(i.*rotamt(2)) * brotset(:,cv3,b3));
m=round(length(zout)/2);
zout=edjoinsmooth(zout(m:end),zout(1:(m-1)));
plotnext(zout,8);

%
%sep=90,90,90: use bidx=1,1,1,1 
%
rotamt=[90 180 270] *pi/180;
b1=1; cv1=4;
b2=1; cv2=4;
b3=1; cv3=4;
b4=1; cv4=4;
z1=xset(:,cv1,b1) + i * yset(:,cv1,b1);
z2=exp(i.*rotamt(1)) * (xset(:,cv2,b2) + i * yset(:,cv2,b2));
z3=exp(i.*rotamt(2)) * (xset(:,cv3,b3) + i * yset(:,cv3,b3));
z4=exp(i.*rotamt(3)) * (xset(:,cv4,b4) + i * yset(:,cv4,b4));
plotnext([z1;z2;z3;z4],2);

zout=edjoinsmooth(srotset(:,cv1,b1),exp(i.*rotamt(1)) * nrotset(:,cv2,b2),0.5);
zout=edjoinsmooth(zout,exp(i.*rotamt(2)) * nrotset(:,cv3,b3),0.5);
zout=edjoinsmooth(zout,exp(i.*rotamt(3)) * erotset(:,cv4,b4),0.5);
m=round(length(zout)/2);
zout=edjoinsmooth(zout(m:end),zout(1:(m-1)),0.5);
plotnext(zout,8);

zout=edjoinsmooth(srotset(:,cv1,b1),exp(i.*rotamt(1)) * erotset(:,cv2,b2),0.5);
zout=edjoinsmooth(zout,exp(i.*rotamt(2)) * srotset(:,cv3,b3),0.5);
zout=edjoinsmooth(zout,exp(i.*rotamt(3)) * erotset(:,cv4,b4),0.5);
m=round(length(zout)/2);
zout=edjoinsmooth(zout(m:end),zout(1:(m-1)),0.5);
plotnext(zout,4);

zout=edjoinsmooth(erotset(:,cv1,b1),exp(i.*rotamt(1)) * brotset(:,cv2,b2),0.5);
zout=edjoinsmooth(zout,exp(i.*rotamt(2)) * srotset(:,cv3,b3),0.5);
zout=edjoinsmooth(zout,exp(i.*rotamt(3)) * nrotset(:,cv4,b4),0.5);
m=round(length(zout)/2);
zout=edjoinsmooth(zout(m:end),zout(1:(m-1)));
plotnext(zout,8);

zout=edjoinsmooth(erotset(:,cv1,b1),exp(i.*rotamt(1)) * brotset(:,cv2,b2),0.5);
zout=edjoinsmooth(zout,exp(i.*rotamt(2)) * brotset(:,cv3,b3),0.5);
zout=edjoinsmooth(zout,exp(i.*rotamt(3)) * srotset(:,cv4,b4),0.5);
m=round(length(zout)/2);
zout=edjoinsmooth(zout(m:end),zout(1:(m-1)),0.5);
plotnext(zout,8);

zout=edjoinsmooth(brotset(:,cv1,b1),exp(i.*rotamt(1)) * brotset(:,cv2,b2),0.5);
zout=edjoinsmooth(zout,exp(i.*rotamt(2)) * brotset(:,cv3,b3),0.5);
zout=edjoinsmooth(zout,exp(i.*rotamt(3)) * brotset(:,cv4,b4),0.5);
m=round(length(zout)/2);
zout=edjoinsmooth(zout(m:end),zout(1:(m-1)),0.5);
plotnext(zout,2);

%
%sep=90,90,90: use bidx=1,1,1,1 
%
rotamt=[90 180 270] *pi/180;
b1=1; cv1=4;
b2=1; cv2=4;
b3=1; cv3=4;
b4=1; cv4=1;
z1=xset(:,cv1,b1) + i * yset(:,cv1,b1);
z2=exp(i.*rotamt(1)) * (xset(:,cv2,b2) + i * yset(:,cv2,b2));
z3=exp(i.*rotamt(2)) * (xset(:,cv3,b3) + i * yset(:,cv3,b3));
z4=exp(i.*rotamt(3)) * (xset(:,cv4,b4) + i * yset(:,cv4,b4));
%plotnext([z1;z2;z3;z4]);

zout=edjoinsmooth(srotset(:,cv1,b1),exp(i.*rotamt(1)) * nrotset(:,cv2,b2));
zout=edjoinsmooth(zout,exp(i.*rotamt(2)) * erotset(:,cv3,b3));
zout=edjoinsmooth(zout,exp(i.*rotamt(3)) * brotset(:,cv4,b4));
m=round(length(zout)/2);
zout=edjoinsmooth(zout(m:end),zout(1:(m-1)));
plotnext(zout,8);

zout=edjoinsmooth(srotset(:,cv1,b1),exp(i.*rotamt(1)) * erotset(:,cv2,b2));
zout=edjoinsmooth(zout,exp(i.*rotamt(2)) * brotset(:,cv3,b3));
zout=edjoinsmooth(zout,exp(i.*rotamt(3)) * brotset(:,cv4,b4));
m=round(length(zout)/2);
zout=edjoinsmooth(zout(m:end),zout(1:(m-1)));
plotnext(zout,8);

zout=edjoinsmooth(brotset(:,cv1,b1),exp(i.*rotamt(1)) * srotset(:,cv2,b2));
zout=edjoinsmooth(zout,exp(i.*rotamt(2)) * erotset(:,cv3,b3));
zout=edjoinsmooth(zout,exp(i.*rotamt(3)) * brotset(:,cv4,b4));
m=round(length(zout)/2);
zout=edjoinsmooth(zout(m:end),zout(1:(m-1)));
plotnext(zout,8);

zout=edjoinsmooth(brotset(:,cv1,b1),exp(i.*rotamt(1)) * brotset(:,cv2,b2));
zout=edjoinsmooth(zout,exp(i.*rotamt(2)) * brotset(:,cv3,b3));
zout=edjoinsmooth(zout,exp(i.*rotamt(3)) * brotset(:,cv4,b4));
m=round(length(zout)/2);
zout=edjoinsmooth(zout(m:end),zout(1:(m-1)));
plotnext(zout,8);


%
%sep=90,90,90: use bidx=1,1,1,1 
%
rotamt=[90 180 270] *pi/180;
b1=1; cv1=4;
b2=1; cv2=1;
b3=1; cv3=4;
b4=1; cv4=1;
z1=xset(:,cv1,b1) + i * yset(:,cv1,b1);
z2=exp(i.*rotamt(1)) * (xset(:,cv2,b2) + i * yset(:,cv2,b2));
z3=exp(i.*rotamt(2)) * (xset(:,cv3,b3) + i * yset(:,cv3,b3));
z4=exp(i.*rotamt(3)) * (xset(:,cv4,b4) + i * yset(:,cv4,b4));
%plotnext([z1;z2;z3;z4]);

zout=edjoinsmooth(brotset(:,cv1,b1),exp(i.*rotamt(1)) * brotset(:,cv2,b2));
zout=edjoinsmooth(zout,exp(i.*rotamt(2)) * brotset(:,cv3,b3));
zout=edjoinsmooth(zout,exp(i.*rotamt(3)) * brotset(:,cv4,b4));
m=round(length(zout)/2);
zout=edjoinsmooth(zout(m:end),zout(1:(m-1)));
plotnext(zout,4);

return


function plotnext(z,rotcount);

global SUBIDX

if ~exist('rotcount','var'),
   rotcount=1;
end

for rr=1:rotcount,
   if isempty(SUBIDX),
      SUBIDX=1;
   else
      SUBIDX=SUBIDX+1;
   end
   
   subplot(23,16,SUBIDX);
   cla
   patch(real([z;z(1)]),imag([z;z(1)]),'k');
   %plot([z;z(1)],'k','LineWidth',2);
   axis([-1.2 1.2 -1.2 1.2]); 
   axis square
   set(gca,'XTickLabel',[],'YTickLabel',[]);
   
   z=exp(i.*pi/4) * z;
   
end
return


% function zout=edjoinsmooth(z1,z2,curvecut,n);
function zout=edjoinsmooth(z1,z2,curvecut,n);

if ~exist('curvecut','var'),
   curvecut=0.3;
end
if ~exist('n','var'),
   n=21;
end

if abs(z1(end)-z2(1))<1e-3,
   zout=[z1;z2];
   return;
end

ppidx=min(find(abs(z1-z1(end))<curvecut));
z1=z1(1:(ppidx-1));
ppidx=max(find(abs(z2-z2(1))<curvecut));
z2=z2((ppidx+1):end);


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


