% function res=rescalesum(batchid,cnfidx[=1],reload[=0])
%
% display spatial tuning for normalized strfs and characterize
% inhibitory tuning
%
% created SVD 6/21/04 - hacked from kvares.m
%
function res=rescalesum(batchid,cnfidx,reload);

dbopen;
sql=['SELECT * FROM sRunData',...
     ' WHERE batch=',num2str(batchid),...
     ' AND cellid<"r0306"',...
     ' ORDER BY cellid'];

rundata=mysql(sql);

if length(rundata)==0,
   fprintf('no runs found for batch %d.\n',batchid);
   return
end

batchdata=dbget('sBatch','',batchid);

RESPATH='/auto/k5/david/tmp/rescale/';
bigresfile=sprintf('%sbatch%d.mat',RESPATH,batchid);
if ~exist('reload','var'),
   reload=0;
end
if ~reload & ~exist(bigresfile,'file'),
   reload=1;
end

if ~exist('cnfidx','var'),
   cnfidx=1;
end
if ismember(batchid,[78 100 103]),
   grestbatch=78;
elseif batchid==101,
   grestbatch=81;
elseif batchid==102,
   grestbatch=80;
else
   disp('what is grbatch?');
   keyboard
end


if reload,
   clear res
   
   for ii=1:length(rundata),
      res(ii).cellid=rundata(ii).cellid;
      resfile=[rundata(ii).respath,rundata(ii).resfile,'.gz'];
      
      if ~exist(resfile,'file'),
         fprintf('%s not found. skipping\n',resfile);
         res(ii).bad=1;
      elseif 0 & sum(strcmp(res(ii).cellid,skipcells))>0,
         fprintf('cell %s marked to skip\n',res(ii).cellid);
         res(ii).bad=1;
      else         
         res(ii).bad=0;
         fprintf('processing strf for cell %s\n',rundata(ii).cellid);
         r=zload(resfile);
         
         nlidx=r.params.nlidxsave;
         res(ii).strf=r.strf(nlidx);
         
         if ~isfield(res(ii).strf,'hspace'),
            h=res(ii).strf.h;
            ttime=sum(h,1)';
            ttimepos=sum(h.*(h>0),1)';
            ttimeneg=sum(h.*(h<0),1)';
      
            cumttime=cumsum(ttime);
            tmax=min(find(cumttime==max(cumttime)));
            if tmax==1,
               tmax=min(find(ttimepos==max(ttimepos)));
            end
            
            %fprintf('tmax=%d\n',tmax);
            tspacemaxpos=sum(h(:,1:tmax),2);
            tspaceall=sum(h,2);
            
            [u,s,v]=svd(h);
            if tspaceall'*u(:,1) > 0;
               ttimesep=v(:,1);
               tspacesep=u(:,1);
            else
               ttimesep=-v(:,1);
               tspacesep=-u(:,1);
            end
            teigvals=diag(s);
            
            res(ii).strf.hspace=tspacesep;
         end
         
         res(ii).hspace=res(ii).strf.hspace;
         res(ii).predxc=r.predres(cnfidx).predxc(nlidx);
         res(ii).prederr=r.predres(cnfidx).prederr(nlidx);
         
         sql=['SELECT * FROM sRunData WHERE cellid="',res(ii).cellid,'"',...
              ' AND batch=',num2str(grestbatch)];
         grrundata=mysql(sql);
         grresfile=[grrundata(1).respath,grrundata(1).resfile,'.gz'];
         z=zload(grresfile);
         nlidx=z.params.nlidxsave;
         grstrf=z.strf(nlidx);
         grkernfmt=grstrf.parms.kernfmt;
         res(ii).hspaceraw=grstrf.hspace;
         
         grtunedata=kern2tune(grstrf);
         res(ii).fitor=grtunedata.fitor(1,1);
         res(ii).fitsf=2^grtunedata.fitsf(1,1);
         
         nlidx=z.params.nlidxsave;
         if grestbatch==78,
            res(ii).grpredxc=z.predres(3).predxc(nlidx);
            res(ii).grprederr=z.predres(3).prederr(nlidx);
         elseif grestbatch==81,
            res(ii).grpredxc=z.predres(2).predxc(nlidx);
            res(ii).grprederr=z.predres(2).prederr(nlidx);
         elseif grestbatch==80,
            res(ii).grpredxc=z.predres(1).predxc(nlidx);
            res(ii).grprederr=z.predres(1).prederr(nlidx);
         end
      end
   end
   
   fprintf('saving to %s\n',bigresfile);
   
   clear r z grstrf hs
   save(bigresfile);
   
else 
   % load previous res data
   fprintf('loading %s\n',bigresfile);
   load(bigresfile);
end

EXAMPLES=0;

if batchid==102,
   % with sketchy review cells:
   skipcells={'93G73A','93G92C','93G98C','93G99A','r0023A',...
              'r0023B','r0024A','r0029A','r0034A',...
              'r0036A','r0036B','r0037A','r0038A','r0049A','r0055B',...
              'r0061A','r0065A','r0072A','r0076B','r0077A',...
              'r0143C','r0158A',...
              'r0162B','r0166B','r0217B','r0259','r0282','r0295',...
              'modelhis','modelhic'};
elseif EXAMPLES,
   skipcells={...
      'e0006','e0012','r0154A','r0158A',...
      'r0170A','r0208D','r0220A','r0230B','r0232A',....
      'r0257','r0275','r0279',...
      'r0034A','r0036A','r0077A','r0081A','r0097A',...
      'r0162B','r0166B','r0168B','r0169B','r0217B',...
      'r0259','r0282','r0295','r0300','r0301'...
      'e0002','e0003',...
      'modelhis','modelhic',...
             };
   
else
   
   skipcells={'r0034A','r0036A','r0077A','r0081A','r0097A',...
              'r0162B','r0166B','r0168B','r0169B','r0217B',...
              'r0259','r0282','r0295','r0300',...
              'e0003',...
              'modelhis','modelhic',...
             };
   %'r0275',
   % these cells don't scale well to 4 cyc/crf vertical
   %              'r0156A','r0162B','r0169B','r0170A','r0208D','r0225C',...
   %              'r0230B','r0230C','r0279','e0006'};
end %

for ii=1:length(rundata),
   if sum(strcmp(res(ii).cellid,skipcells))>0,
      fprintf('cell %s marked to skip\n',res(ii).cellid);
      res(ii).bad=1;
   end
end

kernfmt=res(1).strf.parms.kernfmt;
if ~exist('grkernfmt','var'),
   grkernfmt='pfft';
end

ff=find(~cat(1,res.bad));
predxc=cat(1,res(ff).predxc);
prederr=cat(1,res(ff).prederr);
grpredxc=cat(1,res(ff).grpredxc);
grprederr=cat(1,res(ff).grprederr);

if 1,
   ffkeep=(predxc>prederr);
   %ffkeep=(grpredxc>grprederr);
   ff=ff(ffkeep);
   predxc=predxc(ffkeep);
   grpredxc=grpredxc(ffkeep);
end

fitor=cat(1,res(ff).fitor);
fitsf=cat(1,res(ff).fitsf);

%[psort,tt]=sort(-predxc);
%psort=-psort;
%ff=ff(tt);

hs=cat(2,res(ff).hspace);
hsraw=cat(2,res(ff).hspaceraw);
for ii=1:size(hs,2),
   %hs(:,ii)=hs(:,ii)./norm(hs(:,ii));
   %hsraw(:,ii)=hsraw(:,ii)./norm(hsraw(:,ii));
   hs(:,ii)=hs(:,ii)./std(hs(:,ii));
   hsraw(:,ii)=hsraw(:,ii)./std(hsraw(:,ii));
   
   if hs(:,ii)'*mean(hs,2)<0,
      hs(:,ii)=-hs(:,ii);
   end
end
hs0=hs;



powunbiased=cat(2,res(ff).strf);
powunbiased=cat(2,powunbiased.powunbiased);

hm=mean(hs,2);
hmn=mean(hs.*(hs<0),2);
hmp=mean(hs.*(hs>0),2);

ll=ceil(size(hs,2)/8).*8;
hs(:,size(hs,2)+1:ll)=0;
hsraw(:,size(hsraw,2)+1:ll)=0;
powunbiased(:,size(powunbiased,2)+1:ll)=0;

hs=reshape(hs,size(hs,1),8,ll/8);
hsraw=reshape(hsraw,size(hsraw,1),8,ll/8);
powm=mean(powunbiased,2);
powunbiased=reshape(powunbiased,size(powunbiased,1),8,ll/8);

if strcmp(kernfmt,'pfft'),
   iconside=sqrt(size(hs,1)*2).*[1 1];
else
   iconside=sqrt(size(hs,1)./2).*[1 1];
end

figure(1);
showkern(hs,kernfmt);
colormap(redblue);
fullpage('portrait')

for llidx=1:length(ff),
   subplot(ll/8,8,llidx);
   title(sprintf('%s: %.2f',res(ff(llidx)).cellid,predxc(llidx)));
   xlabel('');
   
   if ~EXAMPLES,
      hold on
      a=axis;
      mx=(a(2)-a(1))/2+1;
      my=(a(4)-a(3))/2+1;
      plot([mx-4 mx+4],[my my],'kx');
      hold off
   end
end

figure(2);
showkern(hsraw,grkernfmt);
colormap(redblue);
fullpage('portrait')

for llidx=1:length(ff),
   subplot(ll/8,8,llidx);
   title(sprintf('%s: %.2f',res(ff(llidx)).cellid,grpredxc(llidx)));
   xlabel('');
end

Xmax=iconside(1);
spacecount=size(hs,1);
phasecount=spacecount./(Xmax.^2./2);
pixcount=Xmax*Xmax;
tk=cat(3,hm./norm(hm),hmp./norm(hmp),...
       hmn./norm(hmn),hmp./norm(hmp)+hmn./norm(hmn),powm);
kcount=size(tk,3);
[cfilt,cfiltconj,cmask,cmaskconj]=gencfilt(Xmax,Xmax);
thsum=zeros(pixcount,2,kcount);
for phaseidx=1:phasecount,
   spacerange=(1:(pixcount/2))+(phaseidx-1).*(pixcount/2);
   thsum(cfilt,1,:)=thsum(cfilt,1,:)+tk(spacerange,:,:);
   thsum(cfiltconj,1,:)=thsum(cfiltconj,1,:)+tk(spacerange,:,:);
end

tth=reshape(thsum(:,1,:),[Xmax Xmax*kcount]);
tth=ifft(conj(fft(tth)));
thsum(:,2,:)=mean([thsum(:,1,:) reshape(tth,prod(iconside),1,kcount)],2);

% optionally set dc gain to 1 cyc/rf gain so that it doesn't
% dominate color range.. .but misleading to do this
%xc=round((Xmax+1)/2);
%ctridx=sub2ind([Xmax Xmax],[xc xc],[xc-1 xc]);
%thsum(ctridx(2),1,3)=thsum(ctridx(1),1,3);

figure(3);
titles={'mean','mean pos','mean neg','mean pos+neg'};
showkern(thsum,'space',iconside,titles);
colormap(redblue(1.1));
fullpage('portrait')

for ii=1:kcount*2,
   subplot(kcount,2,ii);
   hold on
   a=axis;
   mx=(a(2)-a(1))/2+1;
   my=(a(4)-a(3))/2+1;
   plot([mx-4 mx+4],[my my],'kx');
   hold off
end

colorbar

if EXAMPLES
   figure(4);
   ii=find(strcmp({rundata.cellid},'r0210A'));
   ii=find(ii==ff);
   
   hplusminus=[hsraw(:,ii), hsraw(:,ii).*(hsraw(:,ii)>0), ...
               hsraw(:,ii).*(hsraw(:,ii)<0)];
   showkern(hplusminus,kernfmt);
   fullpage('portrait');
   colormap(redblue);
   
   figure(5);
   clf
   
   subplot(3,1,1);
   plot(hs(37,:));
   title('dc vs cells');
   
   subplot(3,1,2);
   plot(hs(38,:));
   title('low vert vs cells');
   
   subplot(3,1,3);
   plot(hs(32,:));
   title('horiz vs cells');
   
   
   
else
   


   addpath /auto/k1/david/code/toolbox/kmeans
   kcount=size(hs0,2);
   
   FLIP=1;
   if FLIP
      [cfilt,cfiltconj,cmask,cmaskconj]=gencfilt(Xmax,Xmax);
      hs1=zeros(pixcount,kcount);
      for phaseidx=1:phasecount,
         spacerange=(1:(pixcount/2))+(phaseidx-1).*(pixcount/2);
         hs1(cfilt,:)=hs1(cfilt,:)+hs0(spacerange,:).*(hs0(spacerange,:)<0);
      end
      hs1(cfiltconj,:)=hs1(cfilt,:);
      
      tth=reshape(hs1(:,:),[Xmax Xmax*kcount]);
      tth=ifft(conj(fft(tth)));
      hs1=(hs1 + reshape(tth,prod(iconside),kcount))./2;
      kkernfmt='space';
      
      hs1=hs1(cfilt,:);
      kkernfmt=kernfmt;
   else
      hs1=hs0.*(hs0<0);
      kkernfmt=kernfmt;
   end
   
   clustercount=6;
   [z, c] = kmeans(hs1,clustercount);
   
   %for ii=1:clustercount,
   %   z(:,ii)=z(:,ii)./abs(min(z(:,ii)));
   %end
   
   figure(4);
   z=cat(2,z,zeros(size(z,1),1));
   showkern(z,kkernfmt);
   subplot(1,size(z,2),size(z,2));
   hist(c,size(z,2)-1);
   colormap(redblue);
   
   
   keyboard
end


%


