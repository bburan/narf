function edstim();

orange=linspace(0,315,8);
shseprange=linspace(45,180,4);
smseprange=linspace(135,45,3);
ocount=length(orange);
sepcount=length(shseprange)+length(smseprange);

iconside=128;
polcount=3;
frames=zeros(iconside,iconside,ocount,sepcount,polcount);
thick=0.2;

for ii=1:ocount,
   for jj=1:sepcount,
      if jj<=length(shseprange),
         ob=180+orange(ii);
         osep=shseprange(jj)/2;
         o1=180+orange(ii)+osep;
         o2=180+orange(ii)-osep;
         
         z1=exp(i*o1/180*pi).*2;
         z2=exp(i*o2/180*pi).*2;
         z1a=exp(i*o1/180*pi).*2+thick./sin(osep./180*pi).*exp(i*ob/180*pi);
         z2a=exp(i*o2/180*pi).*2+thick./sin(osep./180*pi).*exp(i*ob/180*pi);
         if shseprange(jj)==180,
            offset=exp(i*ob/180*pi).*0.15;
         else
            offset=exp(i*ob/180*pi).*(0.25./sin(osep./180*pi));
         end
         
         zout=[z1;0;z2];
         zout0=[zout;flipud(z2a);thick./sin(osep./180*pi).*exp(i*ob/180*pi);...
                flipud(z1a)]-offset;
         zout=zout-offset;
         if 0,
            figure(2);
            clf
            hold off
            plot(zout0);
            hold on;
            plot(zout,'r');
            hold off
            axis equal
            
            keyboard
         end
         
         [x,y]=arc(0,0,o2,o1,2.5);
         zout=[zout;x+i*y]-offset;
         
         [x,y]=arc(0,0,o1,o2,2.5);
         zout2=[flipud(x+i*y);z1;0;z2]-offset;
     else
         ob=180+orange(ii);
         osep=shseprange(jj-length(shseprange))/2
         o1=ob+smseprange(jj-length(shseprange))/2;
         o2=ob-smseprange(jj-length(shseprange))/2;
         
         z1=exp(i*o1/180*pi).*[2; 0.1]+exp(i*(ob+90)/180*pi).*(0.4);
         z2=exp(i*o2/180*pi).*[0.1; 2]+exp(i*(ob-90)/180*pi).*(0.4);
         z1a=exp(i*o1/180*pi).*[2; 0.25]+exp(i*(ob+90)/180*pi).*(0.4)-...
             exp(i*(o1+90)/180*pi)*thick;
         z2a=exp(i*o2/180*pi).*[0.25; 2]+exp(i*(ob-90)/180*pi).*(0.4)-...
             exp(i*(o2-90)/180*pi)*thick;
         
         offset=exp(i*ob/180*pi).*0.15;
         
         zout=edjoinsmooth(z1,z2,0,20)-offset;
         zin=edjoinsmooth(z1a,z2a,0,20)-offset;
         
         zout0=[zout;flipud(zin)];
         if 0,
            figure(2);
            clf
            hold off
            plot(zout);
            hold on;
            plot(z1,'r');
            plot(z2,'r');
            hold off
            axis equal
            keyboard
         end
         
         [x,y]=arc(0,0,angle(z1(1))*180/pi,angle(z2(end))*180/pi,2.5);
         zout2=[flipud(x+i*y); zout];
         
         [x,y]=arc12(real(z2(end)),imag(z2(end)),real(z1(1)),imag(z1(1)),2.5);
         zout=[zout; x+i*y];
      end
      
      figure(1);
      pos=get(1,'Position');
      set(1,'Position',[pos(1) pos(2) iconside iconside]);
      
      clf
      subplot('position',[0 0 1 1]);
      patch([-1.1; -1.1; 1.1; 1.1; -1.1],...
            [-1.1; 1.1; 1.1; -1.1; -1.1],[0 0 0]);
      hold on
      patch(real(zout0),imag(zout0),[1 1 1]);
      %plot([zout],'w','LineWidth',12);
      hold off
      axis([-1.1 1.1 -1.1 1.1]);
      axis square
      axis off
      drawnow
      
      f=getframe(1);
      frames(:,:,ii,jj,1)=f.cdata(:,:,1);
      
      clf
      subplot('position',[0 0 1 1]);
      patch([-1.1; -1.1; 1.1; 1.1; -1.1],...
            [-1.1; 1.1; 1.1; -1.1; -1.1],[0 0 0]);
      hold on
      patch(real(zout),imag(zout),[1 1 1]);
      hold off
      axis([-1.1 1.1 -1.1 1.1]);
      axis square
      axis off
      drawnow
      
      f=getframe(1);
      frames(:,:,ii,jj,2)=f.cdata(:,:,1);
      
      clf
      subplot('position',[0 0 1 1]);
      patch([-1.1; -1.1; 1.1; 1.1; -1.1],...
            [-1.1; 1.1; 1.1; -1.1; -1.1],[0 0 0]);
      hold on
      patch(real(zout2),imag(zout2),[1 1 1]);
      hold off
      axis([-1.1 1.1 -1.1 1.1]);
      axis square
      axis off
      drawnow
      
      f=getframe(1);
      frames(:,:,ii,jj,3)=f.cdata(:,:,1);
   end
end


ff=reshape(frames,iconside,iconside,ocount*sepcount*polcount);

% normalize each frame to have the same contrast
for ii=1:size(ff,3),
   tf=ff(:,:,ii);
   tf=(tf-mean(tf(:)))./std(tf(:)) .*20+127;
   ff(:,:,ii)=tf;
end

!rm /auto/k5/david/sample/stim/edstim.imsm*
writeimfile(ff,'/auto/k5/david/sample/stim/edstim.imsm');

ffs=loadimfile('/auto/k5/david/sample/stim/edstim.imsm',0,0,20,0,20,0,1);
ffs=(ffs-mean(ffs(:)))./std(ffs(:));
ffs=ffs*21+40;

offset=56*0;
ffg1=reshape(ffs(:,:,offset+(1:56)),20*20,8,7);
ffg1=permute(ffg1,[1 3 2]);

figure(2);
clf
showkern(ffg1,'space');
fullpage landscape


keyboard


pffs=movpower(ffs,0,0,1,0.5,0);

ppg1=reshape(pffs(:,offset+(1:56)),size(pffs,1),8,7);
ppg1=permute(ppg1,[1 3 2]);
showkern(ppg1,'pfftgr');

