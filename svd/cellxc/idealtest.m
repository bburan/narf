
dbopen;
batchidx=11;

sql=['SELECT * FROM sRunData where batch=',num2str(batchidx)];
rundata=mysql(sql);

for runidx=1:1,  % length(rundata),
   fprintf('Cell id: %s\n',rundata(runidx).cellid);
   z=zload([rundata(runidx).respath,rundata(runidx).kernfile,'.gz']);
   
   
   latidx=6;
   snrthresh=10;
   
   snr=abs(z.eigH(:,latidx,1)./z.eigHerr(:,latidx));
   sigidx=find(snr>snrthresh);
   
   figure(1);
   clf
   
   for hh=1:length(sigidx),
      [n,x]=hist(z.eigH(sigidx(hh),latidx,z.attuse+1:end),20);
      
   end
   
   return
   
   patchcount=size(z.eigpatches,2);
   targcount=size(z.teigpatches,2);
   attcount=z.attcount;
   attuse=z.attuse;
   spacelim=z.spacelim;
   latidx=5;
   eigrange=15; % could be 2:spacelim;
   eigshow=15;
   
   mp=repmat(mean(z.eigpatches,2),[1 patchcount]);
   sp=repmat(std(z.eigpatches,0,2),[1 patchcount]);
   V=(z.eigpatches-mp) * (z.eigpatches-mp)';
   
   d=zeros(size(z.teigpatches));
   
   for tidx=1:targcount,
      % compute Z score for each eigenvector for each target. preserve
      % sign to remember whether patch is d sd's below or above the mean
      d(:,tidx)=(z.teigpatches(:,tidx)-mp(:,1))./sp(:,1);
   end
   
   xc=zeros(targcount,attcount,spacelim);
   p=zeros(targcount,attuse,spacelim);
   
   for eigidx=eigrange,
      for tidx=2:targcount,
         for attidx=2:attcount,
            cc=xcorr(z.eigH(1:eigidx,latidx,attidx),d(1:eigidx,tidx),...
                     0,'coeff');
            xc(tidx,attidx,eigidx)=cc;
         end
         
         for a2idx=2:attuse,
            sxc=sort(xc(tidx,[a2idx attuse+1:end],eigidx));
            p(tidx,a2idx,eigidx)=1-max(find(xc(tidx,a2idx,eigidx)>=sxc)-1)./...
                length(sxc);
         end
         
         
      end
   end
   
   imagesc(squeeze(p(2:5,2:5,eigshow)),[0 1]); % ,[-1 1]
   %imagesc(squeeze(xc(2:5,2:5,eigshow))); % ,[-1 1]
   xlabel('attidx');
   ylabel('target');
   colorbar
   
   p(2:5,2:5,eigshow)
   
   pause;
end


