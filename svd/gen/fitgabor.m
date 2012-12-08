% function [beta,tfout]=fitgabor(tf);
%
%  tf=xmax X ymax frame to fit
%
%  beta(1)=or
%  beta(2)=orw
%  beta(3)=sf
%  beta(4)=sfw
%  beta(5)=amp;
%  beta(6)=offset;
%
function [beta,tfout]=fitgabor(tf);

xmax=size(tf,1);
ymax=size(tf,2);

x0=round((xmax+1)/2);
xr=(1-x0):(xmax-x0);
yr=xr;

%tft=reshape(tf,xmax*ymax,1).^2;
%tftscale=max(tft);
%tft=tft./tftscale;
tft=reshape(tf,xmax*ymax,1);

[xx,yy]=meshgrid(1:xmax,1:ymax);
x=[reshape(xx,xmax*ymax,1)  reshape(yy,xmax*ymax,1)];

good=0;
while good==0,
   figure(1);
   clf
   imagesc(xr,yr,tf);
   colormap(hot);
   axis image
   %axis off

   beta0=[0.01,1,2,2,1,0.1];
   disp('Click on center of orientation tuning peak.');
   [wx,wy]=ginput(1);
   beta0(1)=-atan(wx/wy)*180/pi;
   beta0(3)=sqrt(wx^2+wy^2);

   disp('Click on rotational edge of tuning peak.');
   [wx,wy]=ginput(1);
   or2=-atan(wx/wy)*180/pi;
   beta0(2)=abs(beta0(1)-or2);
   if beta0(2)>90,
      beta0(2)=180-beta0(2);
   end
   if beta0(2)<0.1,
      beta0(2)=1;
   end
   if beta0(2)>45,
      beta0(2)=44;
   end
   beta0(5)=max(tft(:));
   beta0(6)=beta0(5)/50;
   
   disp('Click on radial edge of tuning peak.');
   [wx,wy]=ginput(1);
   beta0(4)=abs(beta0(3)-sqrt(wx^2+wy^2));
   
%   beta0(1)=input('Orientation (degrees cw from horizontal)? ');
%   beta0(2)=input('Orientation bw (degrees)? ');
%   beta0(3)=input('SF (cycles/2CRF)? ');
%   beta0(4)=input('SF bw? ');
   tf0=reshape(gaborsf(beta0),16,16);
   
   fprintf('Init: Or: %.3f Orw: %.3f SF: %.3f SFw: %.3f A: %.3f Of: %.3f\n',...
		beta0(1),beta0(2),beta0(3),beta0(4),beta0(5),beta0(6));
   %beta=nlinfit(x,tft,'gaborsf',beta0);
   
   betalb=[-Inf  5 0.5 0.2 0   0  ];
   betaub=[ Inf 50 7   3   Inf Inf];
   [beta,resnorm]=lsqcurvefit('gaborsf',beta0,x,tft,betalb,betaub);
   resnorm
   
   %beta(5)=beta(5).*tftscale;
   tfout=reshape(gaborsf(beta),16,16);
   
   clf
   subplot(2,2,1);
   imagesc(xr,yr,tf);
   title('Original kernel');
   axis image
   subplot(2,2,2);
   imagesc(xr,yr,tf0);
   title('Initial guess');
   axis image
   subplot(2,2,3);
   imagesc(xr,yr,tfout);
   title('Gabor fit');
   axis image
   subplot(2,2,4);
   imagesc(xr,yr,abs(tf-tfout),[0,max(max(tf))]);
   title('Difference');
   axis image
   colormap(hot);
   
   fprintf('Final: Or: %.3f Orw: %.3f SF: %.3f SFw: %.3f  A: %.3f Of: %.3f\n',...
		beta(1),beta(2),beta(3),beta(4),beta(5),beta(6));
   fiterror=mserr(tfout,tf);   
   fprintf('Fit error: %.3f\n',fiterror);
   
   good=input('0 to redo: ');
   if isempty(good),
      good=1;
   end
end
