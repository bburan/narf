
% load all that crap
v1predsum

set1=squeeze(predxc(:,[1 3 2],nlinidx,[1 3 2],nloutidx));
goodidx=find(sum(isnan(set1(:,:,1)),2)==0 & ...
             ~strcmp(celllist,'R150B'));
goodcount=length(goodidx);

domcount=2; % not 3 any more
nfft=100;
cmtx=zeros(nfft/2+1,expcount,domcount,cnfcount,goodcount);
ppsth=zeros(expcount,domcount,cnfcount,goodcount);
ptot=zeros(expcount,domcount,cnfcount,goodcount);
psig=zeros(expcount,domcount,cnfcount,goodcount);
pnoise=zeros(expcount,domcount,cnfcount,goodcount);
mse=zeros(expcount,domcount,cnfcount,goodcount);
nmse=zeros(expcount,domcount,cnfcount,goodcount);
xc=zeros(expcount,domcount,cnfcount,goodcount);
xc1=zeros(expcount,domcount,cnfcount,goodcount);
xcasymp=zeros(expcount,domcount,cnfcount,goodcount);

cnfidx=1;    % rev pred
predidx=5;   % rectified output model
for cellidx=1:goodcount,
   cellid=celllist{goodidx(cellidx)};
   for expidx=1:1, % 1:expcount,
      for domidx=2:2, % 1:domcount,
         sql=['SELECT * from sRunData WHERE id=',...
              num2str(bigdata(expidx,domidx,goodidx(cellidx)).runid)];
         trd=mysql(sql);
         z=zload([trd.respath,trd.resfile,'.gz']);
         
         for cnfidx=1:cnfcount,
            
            [cellfiledata,times,batchdata]=...
                cellfiletimes(cellid,z.predbatch{cnfidx});
            
            respfile=[cellfiledata(times(3).fileidx).path,...
                      cellfiledata(times(3).fileidx).respfile];
            ract=respload(respfile,...
                          cellfiledata(times(3).fileidx).respvarname,...
                          cellfiledata(times(3).fileidx).respfiletype,1);
            ract=ract(times(3).start:times(3).stop,2:end);
            rnz=sum(~isnan(ract),2);
            rnz(rnz==0)=inf;
            rcount=min(rnz);
            r=zeros(length(rnz),rcount);
            for ll=1:length(rnz),
               kk=find(~isnan(ract(ll,:)));
               if length(kk)>=rcount,
                  r(ll,:)=ract(ll,kk(1:rcount));
               end
            end
            p=mean(r,2);
            
            ppsth(expidx,domidx,cnfidx,cellidx)=var(p,1);
            ptot(expidx,domidx,cnfidx,cellidx)=mean(var(r,1));
            if rcount>1,
               psig(expidx,domidx,cnfidx,cellidx)=(rcount.*var(p(:),1) - ...
                         ptot(expidx,domidx,cnfidx,cellidx))./(rcount-1);
               pnoise(expidx,domidx,cnfidx,cellidx)=...
                   ptot(expidx,domidx,cnfidx,cellidx)-...
                   psig(expidx,domidx,cnfidx,cellidx);
            end
            
            % extract observed and pred response from results file
            tpred=z.mod_psth{cnfidx}(:,1,1,predidx);
            tact=z.act_resp{cnfidx}(:,1);
            bidx=find(~isnan(tpred) & ~isnan(tact));
            if sum(abs(tpred(1:13)))==0,
               bidx=bidx(find(bidx>13));
            end
            tpred=tpred(bidx);
            tact=tact(bidx);
            
            % compute coherence between obs and pred response
            [cmtx(:,expidx,domidx,cnfidx,cellidx),fs]=...
                cohere(tpred,tact,nfft,72,[],nfft./4);
            if sum(isnan(cmtx(:,expidx,domidx,cnfidx,cellidx)))>0,
               keyboard;
            end
            
            % compute mse and corr coeff between obs and pred response
            tpred=(tpred-mean(tpred))./std(tpred).*std(tact);
            tact=(tact-mean(tact));
            mse(expidx,domidx,cnfidx,cellidx)=var(tpred-tact,1);
            nmse(expidx,domidx,cnfidx,cellidx)=var(tpred-tact,1)./(2*var(tact,1));            
            if size(r,2)>1,
               [m,s,n]=predtrialcurve(tpred,r);
               if length(n)>2,
                  [pp,ss]=polyfit(1./n(2:end),1./m(2:end),1);
               else
                  [pp,ss]=polyfit(1./n,1./m,1);
               end
               fprintf('xc=%.3f asymp xc=%.3f\n',m(end),1./pp(2));
               xc1(expidx,domidx,cnfidx,cellidx)=m(1);
               xc(expidx,domidx,cnfidx,cellidx)=m(end);
               xcasymp(expidx,domidx,cnfidx,cellidx)=1./pp(2);
            else
               disp('n=1!');
               xc(expidx,domidx,cnfidx,cellidx)=...
                   xcov(tpred,tact,0,'coeff');
               xc1(expidx,domidx,cnfidx,cellidx)=xc(expidx,domidx,cnfidx,cellidx);
               xcasymp(expidx,domidx,cnfidx,cellidx)=...
                   xc(expidx,domidx,cnfidx,cellidx);
            end
         end
      end
   end
end

cm=mean(cmtx,5);

figure(3);
clf
plot(fs,squeeze(cm(:,:,2,1)));
legend('Nat','DNat','DGrat');

%[squeeze(ppsth(1,2,1,:)) squeeze(psig(1,2,1,:)) ...
% squeeze(mse(1,2,1,:)) ...
% (squeeze(ppsth(1,2,1,:))-squeeze(mse(1,2,1,:)))./squeeze(psig(1,2,1,:)) ...
% squeeze(xc(1,2,1,:))]

[squeeze(xc1(1,2,1,:)) squeeze(xc(1,2,1,:)) squeeze(xcasymp(1,2,1,:)) ]

