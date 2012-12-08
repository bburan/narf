% function r=tnlres(cellid,batch)
%
% load results from resfile in sRunData and display fits, pred
% results
%
% r=0 if no entries found in db, =1 otherwise
%
function r=tnlres(runidx,batch)

if ~exist('batch','var'),
   batch=80;
   fprintf('assuming batch=80\n');
end

dbopen;
if isnumeric(runidx),
   sql=['SELECT * from sRunData WHERE id=',num2str(runidx)];
   rundata=mysql(sql);
   cellid=rundata.cellid;
else
   cellid=runidx;
   goodbatch=zeros(1,length(batch));
   sql=['SELECT * from sRunData WHERE cellid="',cellid,'"',...
        ' AND batch=',num2str(batch)];
   rundata=mysql(sql);
end

if length(rundata)==0,
   disp('no entries found in db!');
   if nargout>0,
      r=0;
   end
   return
end

outfile=[rundata(1).respath,rundata(1).resfile,'.gz'];
fprintf('outfile: %s\n',outfile);

if strcmp(outfile(end-2:end),'.gz'),
   zload(outfile);
else
   load(outfile);
end

nlrange=[1 4 5];
nshowcount=length(nlrange);
colcount=max([nshowcount 2]);

figure(1);
clf

h=strf(1).h;
h=repmat(h(:,1:(-params.maxlag(1)+1)),[1 1 4]);
h(:,:,2:end)=0;
showkern(h,params.kernfmt);

subplot(4,colcount,1);
plot(mresp);
hold on
plot([transstartidx transstartidx],[min(mresp) max(mresp)],'r-');
plot([transendidx transendidx],[min(mresp) max(mresp)],'r-');
hold off
title(sprintf('%s pfth',params.cellid));

subplot(4,colcount,2);

% quick plot of trans vs sust response
nzidx=find(~isnan(transresp) & ~isnan(sresp) & transresp>0);
slope=mean(sresp(nzidx).*transresp(nzidx))./mean(transresp(nzidx).^2);

scatter(transresp,sresp,'.');
hold on
plot([0 max(transresp)],[0 max(transresp)].*slope,'k--');
hold off

xlabel('trans spikes');
ylabel('sust spikes');
axis equal


mmax=max(max(laggain(:,1,modeloptsfs(nlrange(1)),nlrange)));
mmin=min(min(laggain(:,1,modeloptsfs(nlrange(1)),nlrange)));
for ii=1:length(nlrange),
   subplot(4,colcount,ii+colcount);
   cla
   
   nlidx=nlrange(ii);
   %plot(squeeze(laggain(:,1,:,nlidx)),'k:');
   plot(laggain(:,1,modeloptsfs(nlidx),nlidx),'b-');

   axis([1 fixlen mmin mmax]);
   title(sprintf('%s variable thresh',nlstr{nlidx}));
   ylabel('threshold');
end

mmax=max(max(max(laggain(:,2:3,modeloptsfs(nlrange(1)),nlrange))));
mmin=min(min(min(laggain(:,2:3,modeloptsfs(nlrange(1)),nlrange))));
for ii=1:length(nlrange),
   subplot(4,colcount,ii+colcount*2);
   cla
   
   nlidx=nlrange(ii);
   %plot(squeeze(laggain(:,2,:,nlidx)),'k:');
   plot(laggain(:,2,modeloptsfs(nlidx),nlidx),'b-');
   hold on
   plot(laggain(:,3,modeloptsfs(nlidx),nlidx),'r-');
   hold off
   
   axis([1 fixlen mmin mmax]);
   title(sprintf('%s gain',nlstr{nlidx}));
   ylabel('gain');
end

subplot(4,colcount,10);
mm=max([xc(:)]);
imagesc(squeeze(sum(xc,1)),[0 mm]);
hold on
for nlidx=1:nlcount,
   plot(nlidx,modeloptsfs(nlidx),'x');
end
hold off
title(sprintf('exp pred'));
ylabel('reg parm');
xlabel('output nlidx');
colorbar

subplot(4,colcount,11);
mm=max([allpredxc(:)]);
imagesc(squeeze(sum(allpredxc(:,:,:),1)),[0 mm]);
hold on
for nlidx=1:nlcount,
   plot(nlidx,modeloptsfs(nlidx),'x');
end
hold off
title(sprintf('val pred'));
xlabel('output nlidx');
colorbar

subplot(4,colcount,12);
mm=max([predpower(:)]);
imagesc(predpower,[0 mm]);
hold on
plot(showsig1,showidx1,'x');
%plot(lambda,ones(size(lambda)).*std(aresp(nnidx)),'k');
%plot(lambda,mean(predpower2,2),'r');
%plot(lambda(showidx2),mean(predpower2(showidx2,:)),'kx');
%plot(lambda,mean(vpredpower2,2),'g');
%plot(lambda(showidx3),mean(vpredpower2(showidx3,:)),'kx');
hold off
colorbar


if exist('predfix','var'),
   [modeloptsfs'; squeeze(sum(xc(:,modeloptsfs(3),:),1))'; predxc; predfix]
else
   [modeloptsfs'; squeeze(sum(xc(:,modeloptsfs(3),:),1))'; predxc]
end

colormap(gray);


% subtract mean from fit and val stimuli
tstim=stim(nnidx,:)-repmat(mSall',length(nnidx),1);
tcstim=cdata.stim(vnnidx,:)-repmat(mSall',length(vnnidx),1);

fitresp=sum(respbak(nnidx,:),2);
valresp=sum(vresp(vnnidx,:),2);


h=strf(1).h(:,-params.maxlag(1)+1);

r1=tstim*h;
dcgparms=fitdcgain(r1,fitresp);

h=h.*dcgparms(2);

hp=h.*(h>0);
hn=-h.*(h<0);

%project stim onto postive and negative components of the kernel
r1p=tstim*hp;
r1n=tstim*hn;
v1p=tcstim*hp;
v1n=tcstim*hn;

z1=-mSall'*hp;
z2=-mSall'*hn;

pmin=min([r1p;v1p]);
nmin=min([r1n;v1n]);

%r1p=log(r1p-pmin-0.5);
%r1n=log(r1n-nmin-0.5);
%v1p=log(v1p-pmin-0.5);
%v1n=log(v1n-nmin-0.5);

fitpred1=sum(finalrpred(:,:,1),2);
valpred1=sum(finalvpred(:,:,1),2);
fitpred2=sum(finalrpred(:,:,5),2);
valpred2=sum(finalvpred(:,:,5),2);

posfiterr1=(fitpred1-fitresp).*(fitpred1>fitresp);
negfiterr1=-(fitpred1-fitresp).*(fitpred1<fitresp);
posvalerr1=(valpred1-valresp).*(valpred1>valresp);
negvalerr1=-(valpred1-valresp).*(valpred1<valresp);
posfiterr2=(fitpred2-fitresp).*(fitpred2>fitresp);
negfiterr2=-(fitpred2-fitresp).*(fitpred2<fitresp);
posvalerr2=(valpred2-valresp).*(valpred2>valresp);
negvalerr2=-(valpred2-valresp).*(valpred2<valresp);

if 0,
   aa1=(r1p-pmin)+(r1n-nmin);
   bb1=((r1p-pmin)-(r1n-nmin));
   aa2=(v1p-pmin)+(v1n-nmin);
   bb2=((v1p-pmin)-(v1n-nmin));
   
   zz=[z1-pmin z2-nmin];
elseif 0,
   pmean=mean([r1p; v1p]);
   nmean=mean([r1n; v1n]);

   xx1=[r1p-pmean r1n-nmean; v1p-pmean v1n-nmean];
   n=size(xx1,1);
   
   C=(xx1'*xx1)./n;
   
   yy1=xx1 * C^(-0.5);
   zz=[z1-pmean z2-nmean] * C^(-0.5);
   
   aa1=yy1(1:length(r1p),1);
   bb1=yy1(1:length(r1p),2);
   aa2=yy1(length(r1p)+1:end,1);
   bb2=yy1(length(r1p)+1:end,2);
else
   aa1=r1p;
   bb1=r1n;
   aa2=v1p;
   bb2=v1n;
   zz=[z1 z2 ]
end

scf=5;

figure(2);
clf


xi=linspace(min(aa1),max(aa1),80);
yi=linspace(min(bb1),max(bb1),80)';
zi = griddata(aa1,bb1,fitresp,xi,yi);

if 1,
   g=[.1 .15 .2 .25 .2 .15 .1 ];
   g=g'*g;
   g=g./sqrt(sum(g(:).^2));
   
   zi=conv2(zi,g,'same');
end

imagesc(xi,yi,zi);
%surf(xi,yi,zi);
xlabel('pos component');
ylabel('neg component');
title(sprintf('%s projection onto pos/neg vs obs resp',cellid));
axis xy
colorbar




if 0 & ismember(batch,[80]),
   disp('skipping pred display');
   return
end

figure(3);

% load validation stimulus
sidx=times(3).start;
eidx=times(3).stop;
mov=loadimfix(params.stimfiles{times(3).fileidx},sidx,eidx,params.cellid);

if size(mov,3)>15,
   fprintf('%d frames loaded.\n',size(mov,3));
   %xx=input('restricting time range to 15 frames. first id [1]? ');
   fprintf('restricting time range to 15 frames.\n');
   xx=1;
   if isempty(xx),
      xx=1;
   end
   
   timerange=(1:15)+xx;
else
   timerange=1:size(mov,3);
end

mov=mov(:,:,timerange);

% normalize each frame for cleaner display
mbase=mov(1,1,1);
mov=mov-mbase;
mmax=255-mbase;
mmin=-mbase;
for ii=1:length(timerange),
   tmax=max(max(mov(:,:,ii).*(mov(:,:,ii)>0)));
   tmin=min(min(mov(:,:,ii).*(mov(:,:,ii)<0)));
   if (tmax/mmax)<(tmin/mmin) & tmin<0,
      mov(:,:,ii)=round(mov(:,:,ii).*(mmin/tmin));
   else
      mov(:,:,ii)=round(mov(:,:,ii).*(mmax/tmax));
   end
   
   %[tmax tmin mmax mmin max(max(mov(:,:,ii))) min(min(mov(:,:,ii)))]
end
mov=mov+mbase;

% stretch out validation response and preds into PSTH vectors
actualresp=vresp(timerange,:)';
seppred=allvpred(timerange,:,modeloptsfs(nlrange(1)),nlrange(1))';
inseppred=allvpred(timerange,:,modeloptsfs(nlrange(2)),nlrange(2))';

ticks=zeros(size(actualresp));
ticks(1,:)=1;

goodidx=find(~isnan(actualresp(:)));
actualresp=actualresp(goodidx);
seppred=seppred(goodidx);
inseppred=inseppred(goodidx);
ticks=ticks(goodidx);
ticks=find(ticks).*14-7;

subplot(3,1,1);
imagesc(mov(:,:));
hold on
for ii=1:size(mov,3),
   plot(((ii-1).*size(mov,1)+1) .* [1 1],...
        [size(mov,1) size(mov,1).*1.1],'k-');
end
hold off
colormap(gray);
axis image
axis off
title(sprintf('cell %s review validation stimulus',params.cellid));

ymax=max([actualresp;seppred;inseppred]).*1.1;
ymin=-ymax/10;

hs=subplot(3,1,2);
tbins=(0:(length(actualresp)-1)).*14;
plot(tbins,actualresp,'k--');
hold on
plot(tbins,seppred,'r');
for ii=1:length(ticks),
   plot([ticks(ii) ticks(ii)],[ymin ymin/2],'k-');
end
hold off
axis([0 tbins(end) ymin ymax]);
set(hs,'XTickMode','manual','XTick',[]);
title(sprintf('%s model pred, cc=%.3f',nlstr{nlrange(1)},predxc(nlrange(1))));

hs=subplot(3,1,3);
plot(tbins,actualresp,'k--');
hold on
plot(tbins,seppred,'b:');
plot(tbins,inseppred,'r');
for ii=1:length(ticks),
   plot([ticks(ii) ticks(ii)],[ymin ymin/2],'k-');
end
hold off
axis([0 tbins(end) ymin ymax]);
set(hs,'XTickMode','manual','XTick',[]);
title(sprintf('%s model pred, cc=%.3f',nlstr{nlrange(2)},predxc(nlrange(2))));
xlabel('time (ms)');





return

