
batchid=82;
RESPATH='/auto/k5/david/tmp/kvaparms/';
resfile=sprintf('%sdmsattsum.batch%d.mat',RESPATH,batchid);

load(resfile);

PTHRESH=0.05;

predxc=cat(3,res.predxc);
randxc=cat(3,res.randxc);
pxc=cat(2,res.pxc)';
origlocalsig=pxc(:,2:3)<PTHRESH;

fullxc=cat(3,res.fullxc);
fullp=cat(3,res.fullp);

ftsig=squeeze(fullp(2:4,2,:)<PTHRESH)';
spsig=squeeze(fullp(2:4,3,:)<PTHRESH)';

ftfracs=(squeeze(fullxc(:,2,:))').^2;
spfracs=(squeeze(fullxc(:,3,:))').^2;
ftfracs=ftfracs./repmat(ftfracs(:,end),[1 size(ftfracs,2)]);
spfracs=spfracs./repmat(spfracs(:,end),[1 size(spfracs,2)]);
ftfracs(ftfracs>1)=1;
spfracs(spfracs>1)=1;

ftd=diff(ftfracs(:,[1 3 4]),[],2);
spd=diff(spfracs(:,[1 3 4]),[],2);
ftd(ftd<0)=0;
spd(spd<0)=0;

catidx=ones(length(ftsig),1)*3;
catidx(find(ftsig(:,2)))=2;
catidx(find(ftsig(:,3)))=1;

nv_idx=find(predxc(1,1,:)<randxc(1,1,:).*2);
catidx(nv_idx)=3;

anymodidx=find(sum([ftsig(:,2:3) spsig(:,2:3)],2));
mean([ftfracs(anymodidx,[1 3 4]) spfracs(anymodidx,[1 3 4])])

globaldc=cat(3,res.globaldc);
globalgain=cat(3,res.globalgain);
dc=squeeze(mean(globaldc(2:3,:,:),2));
gg=squeeze(mean(globalgain(2:3,:,:),2));

resp=cat(3,res.rank);
resperr=cat(3,res.resperr);
mresp=cat(3,res.mresp);
strfresp=cat(3,res.strfresp);

figure;

mm=zeros(length(res),2);
rr=zeros(2,2,length(res));
dd=zeros(length(res),1);
for ii=1:length(res),
   
   mm(ii,1)=abs(resp(1,2,ii)-resp(2,2,ii)) > ...
            (resperr(1,2,ii)+resperr(2,2,ii)).*2;
   mm(ii,2)=abs(resp(1,3,ii)-resp(2,3,ii)) > ...
            (resperr(1,3,ii)+resperr(2,3,ii)).*2;
   
   if mm(ii,1) | mm(ii,2),
      cc='r.';
   else
      cc='b.';
   end
   hold on
   if res(ii).resp(1,1)<res(ii).resp(2,1)
      rr(1,1,ii)=resp(1,3,ii); % worse out
      rr(2,1,ii)=resp(1,2,ii); % worse in
      rr(1,2,ii)=resp(2,2,ii); % better out
      rr(2,2,ii)=resp(2,3,ii); % better in
      dd(ii)=mresp(1,3,ii)-mresp(1,2,ii);
   else
      rr(1,1,ii)=resp(2,2,ii); % worse out
      rr(2,1,ii)=resp(2,3,ii); % worse in
      rr(1,2,ii)=resp(1,3,ii); % better out
      rr(2,2,ii)=resp(1,2,ii); % better in
      dd(ii)=mresp(1,2,ii)-mresp(1,3,ii);
   end
   rr(:,:,ii)=rr(:,:,ii)./mresp(1,1,ii);
   if catidx(ii)<3
      plot(rr(:,1,ii),rr(:,2,ii));
   end
end
hold off

catidx2=zeros(size(catidx));
catidx2(find(rr(1,1,:)<1 & rr(2,1,:)<1))=2;
catidx2(find(rr(1,1,:)<1 & rr(2,1,:)>=1))=2;
catidx2(find(rr(1,1,:)>=1 & rr(2,1,:)<1))=2;
catidx2(find(rr(1,1,:)>=1 & rr(2,1,:)>=1))=1;
hist(catidx+catidx2*3)


mean(rr(:,:,catidx2==1),3)
std(rr(:,:,catidx2==1),0,3)./sqrt(sum(catidx2==1))
mean(rr(:,:,catidx2==2),3)
std(rr(:,:,catidx2==2),0,3)./sqrt(sum(catidx2==2))
mean(rr(:,:,catidx2==3),3)
std(rr(:,:,catidx2==3),0,3)./sqrt(sum(catidx2==3))

return


for ii=1:length(res),
   plot((res(ii).resp(1,2)+res(ii).resperr(1,2).*[-1 1]).*60,...
        res(ii).resp(1,3).*[1 1].*60,'r-');
   hold on
   
   plot(res(ii).resp(1,2)*[1 1].*60,...
        (res(ii).resp(1,3)+res(ii).resperr(1,3).*[-1 1]).*60,'r-');
   plot((res(ii).resp(2,2)+res(ii).resperr(2,2).*[-1 1]).*60,...
        res(ii).resp(2,3).*[1 1].*60,'k-');
   plot(res(ii).resp(2,2)*[1 1].*60,...
        (res(ii).resp(2,3)+res(ii).resperr(2,3).*[-1 1]).*60,'k-');
end
hold off


for ii=1:length(res),
   plot((res(ii).resp(1,3)-res(ii).mresp(3)).*[1 1],...
        res(ii).resp(1,2)-res(ii).mresp(2)+res(ii).resperr(1,2).*[-1 1],'-');
   hold on
   plot(res(ii).resp(1,3)-res(ii).mresp(3)+res(ii).resperr(1,3).*[-1 1],...
        (res(ii).resp(1,2)-res(ii).mresp(2)).*[1 1],'-');

   plot((res(ii).resp(2,2)-res(ii).mresp(2)).*[1 1],...
        res(ii).resp(2,3)-res(ii).mresp(3)+res(ii).resperr(2,3).*[-1 1],'-');
   plot(res(ii).resp(2,2)-res(ii).mresp(2)+res(ii).resperr(2,2).*[-1 1],...
        (res(ii).resp(2,3)-res(ii).mresp(3)).*[1 1],'-');
end

plot([-0.5 1],[-0.5 1],'k--');
hold off


        
