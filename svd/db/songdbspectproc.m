chronux_addpath;

load /auto/p1/svd/code/db/songdb_spectdata.mat

pmtx=(s(:,1:101)^(-1))*u';


dbopen;
sql=['SELECT id from dMusic INNER JOIN dEigProj',...
     ' ON dMusic.id=dEigProj.songid WHERE playsec>8 and dEigProj.st1=0',...
     ' ORDER BY dMusic.id'];
songdata=mysql(sql);
songids=cat(1,songdata.id);
fprintf('%d songs without spectro-temporal fingerprint ...\n',...
        length(songids));
%figure(1);
%clf;
for ii=1:length(songids),
   songid=songids(ii);
   
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
         try,
            [y,fs,nbits,opts]=mp3read(songfile);
         catch
            disp('error loading mp3, skipping');
            y=[];
         end

      end
   else
      try,
         [y,fs,nbits,opts]=mp3read(songfile);
      catch
         disp('error loading mp3, skipping');
         y=[];
      end
   end   
   
   if ~isempty(y),
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
   
   mean_power=mean(S,1);
   
   ss=std(S);
   sdata0=S./repmat(ss,[size(sdata,1) 1]);
   sdata0=sdata0-repmat(s_mean,[1 size(sdata0,2)]);
   
   ee=std(Senv);
   edata0=Senv./repmat(ee,[size(edata,1) 1]);
   edata0=edata0-repmat(e_mean,[1 size(edata0,2)]);
   
   A=[mean_power;sdata0;edata0];
   
   stproj=pmtx*A;
   %A=U*S*V'
   %V'=s^(-1)*u'*A;
   
   sql=sprintf(['UPDATE dEigProj set st1=%.6f,st2=%.6f,st3=%.6f,',...
                'st4=%.6f,st5=%.6f WHERE songid=%d'],...
               stproj(1:5),songid)
   mysql(sql);
   
   end % if ~isempty(y)

end
