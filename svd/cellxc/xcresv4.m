% function r=xcresv4(cellid,batch);
function r=xcresv4(cellid,batch);

dbopen;

sql=['SELECT * from sRunData WHERE cellid="',cellid,'"',...
     ' AND batch=',num2str(batch)];
rundata=mysql(sql);

if length(rundata)==0,
   disp('no entries found in db!');
   if nargout>0,
      r=0;
   end
   return
end

outfile=[rundata.respath,rundata.resfile,'.gz'];
fprintf('xcresv4.m: %s\n',outfile);

zload(outfile);

nlidx=params.nlidxsave;
nlidx=1;
if nlidx>length(strf),
   nlidx=1;
end
fprintf('nlidx=%d\n',nlidx);

hf=strf(nlidx).h;

if length(iconside)==1,
   iconside=[20 round(iconside/20)];
end

stephf=cumsum(hf,2);
steptime=sum(stephf,1);
maxidx=min(find(steptime==max(steptime)));
kern=stephf(:,maxidx);

% examples of kernel at different cutoffs
sampcount=6;
sampidx=[round(linspace(1,params.sfscount,sampcount))];

if size(mH,2)==diff(params.maxlag)+1,
   smm=mH(:, (1:size(hf,2))-params.maxlag(1), sampidx,1);
   sms=eH(:, (1:size(hf,2))-params.maxlag(1), sampidx,1) .* sigrange(1);
else
   smm=mH(:, :, sampidx,1);
   sms=eH(:, :, sampidx,1) .* sigrange(1);
end

smd=abs(smm)./(sms+(sms==0));
if ~params.shrinkage, % old--shrinkage filter
   tsf=smm.*(smd>1);
else
   smd=(1-smd.^(-2));
   smd=smd.*(smd>0);
   smd(find(isnan(smd)))=0;
   tsf=smm.*smd;
end

for ii=1:sampcount,
   tl=lambda(sampidx(ii));
   titles{ii}=sprintf('%s: STRF sample %.2f NLxc=%0.3f',...
                      params.cellid,tl,xc(sampidx(ii),1));
end

figure(1);
clf
showkern(tsf,params.kernfmt,iconside,titles,1,16);

figure(2);
clf
subplot(2,1,1);
imagesc(xc(:,:,nlidx,1,1));
hold on
plot(strf(nlidx,end,1).parms.sigfit,strf(nlidx,end,1).parms.sfsfit,'x');
hold off
colormap(hot);
colorbar;
title('fit xc');

figure(3);
clf
if isnumeric(params.predbatch),
   params.predbatch={params.predbatch};
end

if length(params.predbatch)>0,
   bstr=[];
   for ii=1:length(params.predbatch),
      if ~isnan(predxc(ii,nlidx,1,1)),
         bstr=[bstr, sprintf('%d: %.2f/',params.predbatch{ii},predxc(ii,nlidx,1,1))];
      end
   end
   bstr=bstr(1:end-1);
else
   bstr=sprintf('%.2f',predxc(1,nlidx,1,1));
end

sfsfit=strf(nlidx,1,1).parms.sfsfit;
sigfit=strf(nlidx,1,1).parms.sigfit;
titles={sprintf('%s Impulse response (fitxc=%.2f predxc=%s)',...
                params.cellid,xc(sfsfit,sigfit,nlidx,1,1),bstr),...
        sprintf('%s Step response',params.cellid)};

hlen=min([8 size(hf,2)]);
respcount=size(strf,3);
if respcount>1,
   
   hf=cat(3,strf(nlidx,1,:).h);
   for rr=1:respcount,
      if length(params.predbatch)>0,
         bstr=[];
         for ii=1:length(params.predbatch),
            if ~isnan(predxc(ii,1)),
               bstr=[bstr, sprintf('%d: %.2f/',params.predbatch{ii},...
                                   predxc(ii,nlidx,1,rr))];
            end
         end
         bstr=bstr(1:end-1);
      else
         bstr=sprintf('%.2f',predxc(1,nlidx,1,rr));
      end
      sfsfit=strf(nlidx,1,rr).parms.sfsfit;
      sigfit=strf(nlidx,1,rr).parms.sigfit;
      titles{rr}=sprintf('%s Impulse response (fitxc=%.2f predxc=%s)',...
                params.cellid,xc(sfsfit,sigfit,nlidx,1,rr),bstr);
      
   end
   
   showkern(cat(3,hf(:,1:hlen,:)),params.kernfmt,iconside,titles,1,16);
   
else
   showkern(cat(3,hf(:,1:hlen),stephf(:,1:hlen)),params.kernfmt,...
            iconside,titles,1,16);
end

cresp=[];
for fidx=1:length(params.respfiles),
   resploadparms=params.resploadparms;
   resploadparms{2}=1;
   tresp=feval(params.resploadcmd,params.respfiles{fidx},...
               resploadparms{:});
   if ~isempty(params.respfiltercmd),
      tresp=feval(params.respfiltercmd,tresp,params.respfilterparms{:});
   end
   cresp=cat(1,cresp,tresp);
end

%keyboard

cresp=cresp(1:size(predmtx,1),:);

%cresp=cresp(times(3).start:times(3).stop,:);

pred=predmtx;
resp=cresp(:,1);

[ps,psidx]=sort(pred);
bincount=5;
pbedges=round(linspace(1,length(pred)+1,bincount+1));

figure(4);
clf

bincount=8;
pstr={'c.','k.','b.','r.','g.'};
pstr2={'c-','k-','b-','r-','g-'};
for attidx=2:size(cresp,2),
   agoodidx=find(~isnan(cresp(:,attidx)));
   agoodidx=agoodidx(find(~isnan(pred(agoodidx))));
   
   if length(agoodidx)>100,
      scatter(pred(agoodidx(1:100)),cresp(agoodidx(1:100),attidx), ...
              pstr{attidx});
   else
      scatter(pred(agoodidx),cresp(agoodidx,attidx),pstr{attidx});
   end
   xc=xcov(pred(agoodidx),cresp(agoodidx),0,'coeff');
   fprintf('attidx=%d : %d xc=%.3f\n',...
           attidx,length(agoodidx),xc);
   hold on
   
   [ps,psidx]=sort(pred(agoodidx));
   pbedges=round(linspace(1,length(agoodidx)+1,bincount+1));
   
   pp=zeros(bincount,1);
   rr=zeros(bincount,1);
   for bbidx=1:bincount,
      bmatches=psidx(pbedges(bbidx):(pbedges(bbidx+1)-1));
      
      pp(bbidx)=mean(pred(agoodidx(bmatches)));
      rr(bbidx)=mean(resp(agoodidx(bmatches)));
   end
   
   plot(pp,rr,pstr2{attidx});
   
end
hold off
legend('1','2','3','4');