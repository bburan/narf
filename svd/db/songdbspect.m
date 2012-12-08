chronux_addpath;

dbopen;
sql=['SELECT id from dMusic WHERE not(crap)',...
     ' AND not(filemissing) AND not(dup) AND playsec>8'];
songdata=mysql(sql);
songids=cat(1,songdata.id);

runcount=500;
runstep=floor(length(songids)./(runcount+1));
runidx=round((1:runcount).*runstep-runstep./2);

sdata=[];
edata=[];
songlabel={};
figure(1);
clf;
for ii=1:runcount,
   songid=songids(runidx(ii));
   
   sql=['SELECT * FROM dMusic WHERE id=',num2str(songid)];
   songdata=mysql(sql);
   songlabel{ii}=[songdata(1).artist '-' songdata(1).title];
   songfile=songdata(1).file;
   fprintf('song %d: %s\n',ii,songfile);
   
   % load song
   %[Y,FS,NBITS,OPTS] = mp3read(FILE,N,MONO,DOWNSAMP,DELAY)
   if songdata.playsec>360
      try,
         [y,fs,nbits,opts]=mp3read(songfile,16000000);
      catch
         [y,fs,nbits,opts]=mp3read(songfile);
      end
   else
      [y,fs,nbits,opts]=mp3read(songfile);
   end   
   
   % average over stereo channels
   y=mean(y,2);
   
   % downsample to Fs
   Fs=8000;
   y=resample(y,Fs,opts.fmt.nSamplesPerSec);
   
   reps=floor(length(y)./(Fs.*2));
   replen=floor(length(y)./reps);
   x=reshape(y(1:(reps.*replen)),replen,reps);
   x=x(:,round(linspace(1,reps,10)));
   
   enrate=round(Fs./40);
   env=rconv2(abs(y),ones(enrate,1));
   env=env((enrate./2):enrate:end);
   
   params=[];
   params.Fs=Fs;
   params.fpass=[20 Fs./2];
   params.trialave=1;
   [S,f]=mtspectrumc(x,params);
   
   params=[];
   params.Fs=40;
   params.fpass=[0.5 20];
   params.trialave=1;
   [Senv,fenv]=mtspectrumc(env,params);
   
   fbins=50;
   fbinsize=round(length(S)./fbins);
   S=rconv2(S,ones(fbinsize,1));
   S=S(round(fbinsize./2):fbinsize:end);
   f=f(round(fbinsize./2):fbinsize:end);
   
   ebins=50;
   ebinsize=round(length(Senv)./ebins);
   Senv=rconv2(Senv,ones(ebinsize,1));
   Senv=Senv(round(ebinsize./2):ebinsize:end);
   fenv=fenv(round(ebinsize./2):ebinsize:end);
   
   sfigure(1);
   if runcount>50,
      clf
      subplot(2,1,1);
      semilogy(f,S);
      [pp,bb,ee]=fileparts(songfile);
      ht=title(bb);
      set(ht,'Interpreter','none');
      
      subplot(2,1,2);
      semilogy(fenv,Senv,'r');
      title('envelope spectrum');
      set(gca,'Xticklabel',[]);
   else
      
      subplot(runcount./5,5,ii);
      semilogy(S.*10000);
      [pp,bb,ee]=fileparts(songfile);
      ht=title(bb);
      set(ht,'Interpreter','none');
      %set(gca,'Xticklabel',[]);
      
      hold on;
      semilogy(Senv,'r');
      set(gca,'Xticklabel',[]);
      hold off
   end
   
   edata=[edata Senv];
   sdata=[sdata S];
   
   if ii>2,
      mean_power=mean(sdata,1);
      
      ss=std(sdata);
      sdata0=sdata./repmat(ss,[size(sdata,1) 1]);
      s_mean=mean(sdata0,2);
      sdata0=sdata0-repmat(s_mean,[1 size(sdata0,2)]);
      
      ee=std(edata);
      edata0=edata./repmat(ee,[size(edata,1) 1]);
      e_mean=mean(edata0,2);
      edata0=edata0-repmat(e_mean,[1 size(edata0,2)]);
      
      A=[mean_power;sdata0;edata0];
      [u,s,v]=svd(A);
      
      sfigure(2);
      clf
      subplot(2,2,1);
      plot(f,s_mean);
      title('average spectrum');
      subplot(2,2,3);
      plot(fenv,e_mean);
      title('average rate spectrum');
      subplot(2,2,2);
      plot(u(:,1:3));
      title('1st three pcs');
      
      sfigure(3);
      clf
      plot(v(:,1),v(:,2),'.');
      for jj=1:ii,
         ht=text(v(jj,1)+0.002,v(jj,2),songlabel{jj});
      end
      
   end
   
   drawnow

end

songidused=songids(runidx)
save /auto/p1/svd/code/db/songdb_spectdata.mat edata sdata s_mean ...
   e_mean u s v songlabel songidused
