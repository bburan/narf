% function [tunedata,fitdata]=xctune(outfiles/cellid,batch)
%
function [tunedata,fitdata]=xctune(cellid,batch)

dbopen;

if nargout>0,
   tunedata=[];
   fitdata=[];
end

if iscell(cellid),
   resfiles=cellid;
   batchcount=length(resfiles);
elseif strcmp(cellid(end-3:end),'.mat') | ...
      strcmp(cellid(end-2:end),'.gz'),
   resfiles={cellid};
   batchcount=1;
else
   goodbatch=zeros(1,length(batch));
   batchcount=0;
   resfiles={};
   for ii=1:length(batch),
      sql=['SELECT * from sRunData WHERE cellid="',cellid,'"',...
           ' AND batch=',num2str(batch(ii))];
      trd=mysql(sql);
      if ~isempty(trd),
         goodbatch(ii)=1;
         batchcount=batchcount+1;
         resfiles{batchcount}=[trd.respath,trd.resfile,'.gz'];
         rundata(batchcount)=trd;
      end
   end
   batch=batch(find(goodbatch));
   batchcount=length(batch);
end

fprintf('xctune.m:  %d files:\n',batchcount);

for batchidx=1:batchcount
   fprintf(' %s\n',resfiles{batchidx});
   if strcmp(resfiles{batchidx}(end-2:end),'.gz'),
      zload(resfiles{batchidx});
   else
      load(resfiles{batchidx});
   end
   
   nlidx=params.nlidxsave;
   if nlidx>length(strf),
      nlidx=1;
   end
   fprintf('chose to show nlidx=%d\n',nlidx);
   
   tfitdata=[];
   for bootidx=1:params.fitboot
      tunedata(bootidx,batchidx)=kern2tune(strf(nlidx,bootidx,1));
      
      tfitdata(bootidx).orfit=fitflatcos(tunedata(bootidx,batchidx).obins,...
                                         tunedata(bootidx,batchidx).por(:,1));
      tfitdata(bootidx).sffit=fitgauss1d(log2(tunedata(bootidx,batchidx).sfbins),...
                                         tunedata(bootidx,batchidx).psf(:,1));
      tfitdata(bootidx).timefit=...
          fitexpdecaydiff(tunedata(bootidx,batchidx).tbins,...
                          tunedata(bootidx,batchidx).ptime(:,1));
      
   end
   
   fitdata(batchidx).orfit=mean(cat(2,tfitdata.orfit),2)';
   fitdata(batchidx).orfiterr=std(cat(2,tfitdata.orfit),1,2)' .* ...
       sqrt(params.fitboot);
   fitdata(batchidx).sffit=mean(cat(2,tfitdata.sffit),2)';
   fitdata(batchidx).sffiterr=std(cat(2,tfitdata.sffit),1,2)' .* ...
       sqrt(params.fitboot);
   fitdata(batchidx).timefit=mean(cat(2,tfitdata.timefit),2)';
   fitdata(batchidx).timefiterr=std(cat(2,tfitdata.timefit),1,2)' .* ...
       sqrt(params.fitboot);
   
end

obincount=length(tunedata(1,batchidx).obins);

for batchidx=1:batchcount
   figure(1);
   
   subplot(4,batchcount,batchidx);
   mor=cat(3,tunedata(:,batchidx).por);
   eor=std(mor(:,1,:),1,3) .* sqrt(params.fitboot-1);
   mor=mean(mor(:,1,:),3);
   errorbar(tunedata(1,batchidx).obins,mor,eor);
   
   hold on
   plot(tunedata(1,batchidx).obins,...
        flatcos(fitdata(batchidx).orfit,tunedata(1,batchidx).obins),'r--');
   hold off
   
   title(sprintf('%s bat %d or',cellid,batch(batchidx)));
   
   subplot(4,batchcount,batchidx+batchcount);
   msf=cat(3,tunedata(:,batchidx).psf) ;
   esf=std(msf(:,1,:),1,3) .* sqrt(params.fitboot-1);
   msf=mean(msf(:,1,:),3);
   errorbar(tunedata(1,batchidx).sfbins,msf,esf);
   title(sprintf('%s bat %d sf',cellid,batch(batchidx)));
   
   hold on
   plot(tunedata(1,batchidx).sfbins,...
        gauss1(fitdata(batchidx).sffit,log2(tunedata(1,batchidx).sfbins)),'r--');
   hold off
   
   subplot(4,batchcount,batchidx+batchcount*2);
   mtime=cat(3,tunedata(:,batchidx).ptime) ;
   etime=std(mtime(:,1,:),1,3) .* sqrt(params.fitboot-1);
   mtime=mean(mtime(:,1,:),3);
   errorbar(tunedata(1,batchidx).tbins,mtime,etime);
   title(sprintf('%s bat %d time',cellid,batch(batchidx)));
   
   hold on
   plot(tunedata(1,batchidx).tbins,...
        expdecaydiff(fitdata(batchidx).timefit,tunedata(1,batchidx).tbins),'r--');
   hold off
end

for batchidx=1:batchcount
   figure(2);
   
   subplot(4,batchcount,batchidx);
   mor=cat(3,tunedata(:,batchidx).por);
   imagesc(squeeze(mor(:,1,:)));
   title(sprintf('%s bat %d or',cellid,batch(batchidx)));
   
   subplot(4,batchcount,batchidx+batchcount);
   msf=cat(3,tunedata(:,batchidx).psf) ;
   imagesc(squeeze(msf(:,1,:)));
   title(sprintf('%s bat %d sf',cellid,batch(batchidx)));
   
   subplot(4,batchcount,batchidx+batchcount*2);
   mtime=cat(3,tunedata(:,batchidx).ptime) ;
   imagesc(squeeze(mtime(:,1,:)));
   title(sprintf('%s bat %d time',cellid,batch(batchidx)));
   
end

fprintf('TUNING SUMMARY:\n');
fprintf('%5s %5s %5s  %5s %5s  %5s %5s\n',...
        'BATCH','OR','ORBW','SF','SFBW','LAT','ARAT');
for batchidx=1:batchcount,
   fprintf('%5d ',batch(batchidx));
   fprintf('%5.1f %5.1f  %5.2f %5.2f  %5.1f %5.2f\n',...
           fitdata(batchidx).orfit(1:2),fitdata(batchidx).sffit(1:2),...
           fitdata(batchidx).timefit(1),...
           (fitdata(batchidx).timefit(3)-fitdata(batchidx).timefit(6)) ./ ...
           fitdata(batchidx).timefit(3));
end

keyboard








