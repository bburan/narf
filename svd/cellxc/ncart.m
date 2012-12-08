

indexfile='/auto/k5/hayden/ncart/index';
pixfn='/tmp/ncart.imsm';

if 1,
   repcount=1;
   
   [indexfn,indexpath]=basename(indexfile);
   [framecount,iconside]=pypestim2imsmraw(indexpath,indexfn,pixfn,repcount);
end

framecount=imfileinfo(pixfn);

loadidx=1:4:framecount;
mov=loadimframes(pixfn,loadidx,20,0,20,0,1);

% subtract off background
mov=mov-mov(1);

framecount=size(mov,3);
fmov=movpower(mov,0,0,1,1,0);
mov=reshape(mov,20*20,framecount);


showstim(mov(:,1:100),'space',[20 20],10,10,'ncart space domain');
set(gcf,'PaperOrientation','portrait','PaperPosition',[0.25 0.25 8 10.5]);

showstim(fmov(:,1:100),'pfft',[200],10,10,'ncart fp domain');
set(gcf,'PaperOrientation','portrait','PaperPosition',[0.25 0.25 8 10.5]);

mstim=mean(fmov,2);
fmov=fmov-repmat(mstim,[1 framecount]);

sa=fmov*fmov';
[u,s,v]=svd(sa);

showstim(u(:,1:100),'pfft',[200],10,10,'ncart eigs fp domain');
set(gcf,'PaperOrientation','portrait','PaperPosition',[0.25 0.25 8 10.5]);
