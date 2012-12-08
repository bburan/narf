% function [tunedata,fitdata]=xcfittune(cellid,batch)
%
function [tunedata,fitdata]=xcfittune(cellid,batch)

dbopen;

if nargout>0,
   tunedata=[];
   fitdata=[];
end

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
      kernfiles{batchcount}=[trd.respath,trd.kernfile];
      rundata(batchcount)=trd;
   end
end
batch=batch(find(goodbatch));
batchcount=length(batch);

fprintf('xcfittune.m:  %d file(s):\n',batchcount);

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
   fprintf('using nlidx=%d\n',nlidx);
   
   tfitdata=[];
   
   for bootidx=1:params.fitboot
      fprintf('bootidx=%d:\n',bootidx);
      if bootidx==1 & batchidx==1,
         tunedata=kern2tune(strf(nlidx,bootidx,1));
      else
         tunedata(bootidx,batchidx)=kern2tune(strf(nlidx,bootidx,1));
      end
      
      tfitdata(bootidx).orfit=fitflatcos(tunedata(bootidx,batchidx).obins,...
                                         tunedata(bootidx,batchidx).por(:,1));
      tfitdata(bootidx).sffit=fitgauss1d(...
         log2(tunedata(bootidx,batchidx).sfbins),...
         tunedata(bootidx,batchidx).psf(:,1));
      tfitdata(bootidx).timefit=...
          fitexpdecaydiff(tunedata(bootidx,batchidx).tbins,...
                          tunedata(bootidx,batchidx).ptime(:,1));
      
      h=strf(nlidx,bootidx,1).h;
      hsum=cumsum(h);
      ttime=squeeze(mean(hsum,1));
      if max(ttime)>0,
         ttmax=min(find(ttime==max(ttime)));
      else
         ttmax=min(find(abs(ttime)==max(abs(ttime))));
      end
      
      Xmax=sqrt(size(h,1)*2);
      [cfilt,cfiltconj,cmask,cmaskconj]=gencfilt(Xmax,Xmax);
      orsf=zeros(Xmax,Xmax);
      orsf(cfilt)=h(:,ttmax);
      orsf(cfiltconj)=h(:,ttmax);
      
      [beta,beta0]=fitgaussfp(orsf);
      tfitdata(bootidx).orsffit=beta;
      
      orsf(cfilt)=h(:,end);
      orsf(cfiltconj)=h(:,end);
      tfitdata(bootidx).orsffitsust=fitgaussfp(orsf);
      
      dbsetqueue;
      
      %keyboard
   end
   
   fitdata(batchidx).orfit=mean(cat(2,tfitdata.orfit),2)';
   to=cat(2,tfitdata.orfit);
   fitdata(batchidx).orfiterr=std(to,1,2)' .* ...
       sqrt(params.fitboot);
   if fitdata(batchidx).orfiterr(1)>150,
      to(1,:)=mod(to(1,:)+90,180);
      fitdata(batchidx).orfiterr(1)=std(to(1,:),1)'.*sqrt(params.fitboot);
   end
   
   fitdata(batchidx).orsffit=mean(cat(2,tfitdata.orsffit),2)';
   to=cat(2,tfitdata.orsffit);
   fitdata(batchidx).orsffiterr=std(to,1,2)' .* ...
       sqrt(params.fitboot);
   if fitdata(batchidx).orsffiterr(1)>pi/6*5,
      to(1,:)=mod(to(1,:)+pi/2,pi);
      fitdata(batchidx).orsffit(1)=mean(to(1,:));
      fitdata(batchidx).orsffiterr(1)=std(to(1,:),1)'.*sqrt(params.fitboot);
   end
   
   fitdata(batchidx).orsffitsust=mean(cat(2,tfitdata.orsffitsust),2)';
   to=cat(2,tfitdata.orsffitsust);
   fitdata(batchidx).orsffitsusterr=std(to,1,2)' .* ...
       sqrt(params.fitboot);
   if fitdata(batchidx).orsffitsusterr(1)>pi/6*5,
      to(1,:)=mod(to(1,:)+pi/2,pi);
      fitdata(batchidx).orsffitsust(1)=mean(to(1,:));
      fitdata(batchidx).orsffitsusterr(1)=...
          std(to(1,:),1)'.*sqrt(params.fitboot);
   end
   
   fitdata(batchidx).sffit=mean(cat(2,tfitdata.sffit),2)';
   fitdata(batchidx).sffiterr=std(cat(2,tfitdata.sffit),1,2)' .* ...
       sqrt(params.fitboot);
   
   fitdata(batchidx).timefit=mean(cat(2,tfitdata.timefit),2)';
   fitdata(batchidx).timefiterr=std(cat(2,tfitdata.timefit),1,2)' .* ...
       sqrt(params.fitboot);
   
   fprintf('Saving fit data to kernfile %s...\n',kernfiles{batchidx});
   save(kernfiles{batchidx},'fitdata','tunedata','params','strf');
   
end


