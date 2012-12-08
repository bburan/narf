% function res=kvaparms(batchid,reload[=0],recalc[=0])
%
% created SVD 6/21/04 - hacked from kvares.m
% modified SVD 3/20/05 - added tuning curve analysis in recalc section
%
function res=kvaparms(batchid,reload,recalc);
dbopen;
sql=['SELECT * FROM sRunData',...
     ' WHERE batch=',num2str(batchid),...
     ' AND not(cellid like "m0000")',...
     ' AND not(cellid like "v0170")',...
     ' AND not(cellid like "v0172")',...
     ' ORDER BY cellid'];
rundata=mysql(sql);

if length(rundata)==0,
   fprintf('no runs found for batch %d.\n',batchid);
   return
end

batchdata=dbget('sBatch',batchid);

RESPATH='/home/vdata/tmp/kvaparms/';

resfile=sprintf('%sbatch%d.mat',RESPATH,batchid);
if ~exist('reload','var'),
   reload=0;
end
if ~exist('recalc','var'),
   recalc=0;
end
if ~reload & ~exist(resfile,'file'),
   reload=1;
end

if reload,
   res.cellid={};
   res.attxc=nan.*zeros(length(rundata),2,4);
   
   for ii=1:length(rundata),
      res.cellid{ii}=rundata(ii).cellid;
      
      if ~exist([rundata(ii).respath,rundata(ii).kernfile,'.gz'],'file'),
         fprintf('%s not found. skipping\n',...
                 [rundata(ii).respath,rundata(ii).kernfile,'.gz']);
      elseif strcmp(batchdata.matcmd,'xcdms'),
         
         %dms summary
         if ii==1,
            res.resp=zeros(2,3,length(rundata));
            res.mresp=zeros(3,1,length(rundata));
            res.rank=zeros(2,3,length(rundata));
            res.strfresp=zeros(2,3,length(rundata));
            res.strfrank=zeros(2,3,length(rundata));
            res.strfmean=zeros(3,1,length(rundata));
            res.predxc=zeros(3,4,length(rundata));
            res.pxc=zeros(4,1,length(rundata));
            res.randxc=zeros(1,1,length(rundata));
            res.dcxc=zeros(4,1,length(rundata));
            res.locxc=zeros(4,1,length(rundata));
            res.onexc=zeros(4,2,length(rundata));
            res.targdist=zeros(1,1,length(rundata));
         end
         
         tres=dmsres(rundata(ii).cellid,batchid);
         res.cellid{ii}=rundata(ii).cellid;
         res.resp(:,:,ii)=tres.resp;
         res.mresp(:,:,ii)=tres.mresp(:);
         res.rank(:,:,ii)=tres.rank;
         res.strfresp(:,:,ii)=tres.strfresp;
         res.strfrank(:,:,ii)=tres.strfrank;
         res.strfmean(:,:,ii)=tres.strfmean;
         res.predxc(:,:,ii)=tres.predxc;
         res.pxc(:,:,ii)=tres.pxc;
         res.randxc(:,:,ii)=tres.randxc;
         res.dcxc(:,:,ii)=tres.dcxc;
         res.locxc(:,:,ii)=tres.locxc;
         res.onexc(:,:,ii)=tres.onexc;
         res.targdist(:,:,ii)=tres.targdist;
         
      elseif strcmp(batchdata.matcmd,'kernvsatt2'),
         
         %fvvs summary
         tres=kvares(rundata(ii).cellid,batchid);
         
         if ii==1,
            attcount=length(tres.mresp);
            res.mresp=zeros(attcount,1,length(rundata)).*nan;
            res.ncount=zeros(attcount,1,length(rundata)).*nan;
            res.dcerr=zeros([size(tres.dcerr) length(rundata)]).*nan;
            res.strfresp=zeros([size(tres.strfresp),length(rundata)]).*nan;
            res.strfrespnodc=zeros([size(tres.strfresp),length(rundata)]).*nan;
            res.strfresppos=zeros([size(tres.strfresp),length(rundata)]).*nan;
            res.strfrespneg=zeros([size(tres.strfresp),length(rundata)]).*nan;
            res.strfrank=zeros(attcount-1,attcount,length(rundata)).*nan;
            res.strfmean=zeros(attcount,1,length(rundata)).*nan;
            res.predxc=zeros([size(tres.predxc),length(rundata)]).*nan;
            res.pxc=zeros([size(tres.pxc),length(rundata)]).*nan;
            res.randxc=zeros(1,1,length(rundata));
            res.predxccross=zeros([size(tres.predxccross),...
                    length(rundata)]).*nan;
            res.predxccrossrand=zeros([size(tres.predxccrossrand),...
                    length(rundata)]).*nan;
            
            res.targsim=zeros([size(tres.targsim),length(rundata)]).*nan;
            res.prefsim=zeros([size(tres.prefsim),length(rundata)]).*nan;
            res.predsim=zeros([size(tres.predsim),length(rundata)]).*nan;
            res.dcsim=zeros([size(tres.dcsim),length(rundata)]).*nan;
            res.predin=zeros([size(tres.predin),length(rundata)]).*nan;
            res.predout=zeros([size(tres.predout),length(rundata)]).*nan;
         end
         
         res.cellid{ii}=rundata(ii).cellid;
         if length(tres.pxc)==size(res.pxc,1) &...
                length(tres.mresp)==size(res.mresp,1),
            %disp('length(tres.pxc)==size(res.pxc,1): counting');
            res.mresp(:,:,ii)=tres.mresp(:);
            res.ncount(:,:,ii)=tres.ncount(:);
            res.dcerr(:,:,ii)=tres.dcerr;
            res.strfresp(:,:,ii)=tres.strfresp;
            res.strfrespnodc(:,:,ii)=tres.strfrespnodc;
            res.strfresppos(:,:,ii)=tres.strfresppos;
            res.strfrespneg(:,:,ii)=tres.strfrespneg;
            res.strfrank(:,:,ii)=tres.strfrank;
            res.strfmean(:,:,ii)=tres.strfmean;
            res.predxc(:,:,ii)=tres.predxc;
            res.pxc(:,:,ii)=tres.pxc;
            res.randxc(:,:,ii)=tres.randxc;
            res.predxccross(:,:,ii)=tres.predxccross;
            res.predxccrossrand(:,:,ii)=tres.predxccrossrand;

            res.targsim(:,:,ii)=tres.targsim;
            res.prefsim(:,:,ii)=tres.prefsim;
            res.predsim(:,:,ii)=tres.predsim;
            res.dcsim(:,:,ii)=tres.dcsim;
            res.predin(:,:,ii)=tres.predin;
            res.predout(:,:,ii)=tres.predout;
         end
      elseif 1,
         % in vs out prediction summary
         if ii==1,
            res.attxc=nan.*zeros(length(rundata),2,4);
         end
         
         z=zload([rundata(ii).respath,rundata(ii).kernfile,'.gz']);
         fprintf('cellid %s:\n',rundata(ii).cellid);
         if isfield(z,'nlusecount'),
            for nlidx=1:z.nlusecount,
               fprintf('nlidx=%d: ',z.nluse(nlidx));
               for attidx=2:z.attcount,
                  for att2=attidx+1:z.attcount,
                     p1=z.valattpreds(:,attidx,nlidx);
                     p2=z.valattpreds(:,att2,nlidx);
                     a1idx=find(~isnan(z.bresp(z.anyokidx,1,attidx)));
                     a2idx=find(~isnan(z.bresp(z.anyokidx,1,att2)));
                     
                     pin=xcov(z.bresp(z.anyokidx([a1idx;a2idx]),:,1),...
                              [p1(a1idx);p2(a2idx)],0,'coeff');
                     pout=xcov(z.bresp(z.anyokidx([a1idx;a2idx]),:,1),...
                               [p2(a1idx);p1(a2idx)],0,'coeff');
                     fprintf(' %6.3f/%6.3f\n',pin,pout);
                     
                     res.attxc(ii,1,nlidx)=pin;
                     res.attxc(ii,2,nlidx)=pout;
                  end
               end
            end
         end
      else
         kvares(rundata(ii).cellid,batchid);
         
         figure(2);
         drawnow;
         print -f2 -dpsc -Pgcolor
         
         figure(1);
         drawnow;
         print -f1 -dpsc -Pgcolor
      end
   end
   fprintf('saving %s\n',resfile);
   
   clear recalc
   
   save(resfile);
   recalc=1;
else
   % load previous res data
   fprintf('loading %s\n',resfile);
   load(resfile);
end


if strcmp(batchdata.matcmd,'xcdms'),
   % dms summary stuff
   
   %obsbaseresp=[squeeze(res.resp(1,1,:)); squeeze(res.resp(2,1,:))];
   obsbaseresp=[squeeze(res.resp(1,3,:)); squeeze(res.resp(2,2,:))];
   strfbaseresp=[squeeze(res.strfresp(1,1,:)); squeeze(res.strfresp(2,1,:))];
   obsattresp=[squeeze(res.resp(1,2,:)); squeeze(res.resp(2,3,:))];
   obsoutresp=[squeeze(res.resp(1,3,:)); squeeze(res.resp(2,2,:))];
   strfattresp=[squeeze(res.strfresp(1,2,:)); squeeze(res.strfresp(2,3,:))];
   strfoutresp=[squeeze(res.strfresp(1,3,:)); squeeze(res.strfresp(2,2,:))];
   
   %obsbaserank=[squeeze(res.rank(1,1,:)); squeeze(res.rank(2,1,:))];
   obsbaserank=[squeeze(res.rank(1,3,:)); squeeze(res.rank(2,2,:))];
   strfbaserank=[squeeze(res.strfrank(1,1,:)); squeeze(res.strfrank(2,1,:))];
   obsattrank=[squeeze(res.rank(1,2,:)); squeeze(res.rank(2,3,:))];
   strfattrank=[squeeze(res.strfrank(1,2,:)); squeeze(res.strfrank(2,3,:))];
   strfoutrank=[squeeze(res.strfrank(1,3,:)); squeeze(res.strfrank(2,2,:))];
   
   %obsbaserankdiff=squeeze(res.rank(1,1,:)-res.rank(2,1,:));
   %obsbaserankmin=squeeze(min([res.rank(1,1,:); res.rank(2,1,:)]));
   %obsbaserankmean=squeeze(mean([res.rank(1,1,:); res.rank(2,1,:)]));
   obsbaserankdiff=squeeze(res.rank(1,3,:)-res.rank(2,2,:));
   obsbaserankmin=squeeze(min([res.rank(1,3,:); res.rank(2,2,:)]));
   obsbaserankmean=squeeze(mean([res.rank(1,3,:); res.rank(2,2,:)]));
   obsattrankdiff=squeeze(res.rank(1,2,:)-res.rank(2,3,:));
   obsattrankmin=squeeze(min([res.rank(1,2,:); res.rank(2,3,:)]));
   obsattrankmean=squeeze(mean([res.rank(1,2,:); res.rank(2,3,:)]));
   
   strfbaserankdiff=squeeze(res.strfrank(1,1,:)-res.strfrank(2,1,:));
   strfbaserankmin=squeeze(min([res.strfrank(1,1,:); res.strfrank(2,1,:)]));
   strfbaserankminsign=squeeze(res.strfrank(1,1,:)<res.strfrank(2,1,:))*2-1;
   strfbaserankmean=squeeze(mean([res.strfrank(1,1,:); res.strfrank(2,1,:)]));
   strfattrankdiff=squeeze(res.strfrank(1,2,:)-res.strfrank(2,3,:));
   strfattrankmin=squeeze(min([res.strfrank(1,2,:); res.strfrank(2,3,:)]));
   strfattrankminsign=squeeze(res.strfrank(1,2,:)<res.strfrank(2,3,:))*2-1;
   strfattrankmean=squeeze(mean([res.strfrank(1,2,:); res.strfrank(2,3,:)]));
   
   obsmeandiff=squeeze((res.mresp(2,1,:)-res.mresp(3,1,:))./...
                       res.mresp(1,1,:));
   strfmeandiff=squeeze((res.strfmean(2,1,:)-res.strfmean(3,1,:))./...
                        res.strfmean(1,1,:));
   strfbasemean=squeeze(res.strfmean(1,1,:));
   strfattmean=squeeze(res.strfmean(2:3,1,:))';
   dcbase=squeeze(res.mresp(1,1,:));
   dcattchange=[obsmeandiff; -obsmeandiff];
   
   localpredxcste=squeeze(res.predxc(1,2,:)./res.predxc(3,2,:));
   localpredxc=squeeze(res.predxc(1,2,:));
   basepredxc=squeeze(res.predxc(1,1,:));
   
   dcpredxc=squeeze(res.dcxc(2,:,:));
   locpredxc0=squeeze(res.locxc(1,:,:));
   locpredxc=squeeze(res.locxc(2,:,:));
   targdist=squeeze(res.targdist);
   
   baseonexc=squeeze(res.onexc(1,:,:))';
   baseonexcdiff=baseonexc(:,1)-baseonexc(:,2);
   baseonexcdiffnorm=(baseonexc(:,1)-baseonexc(:,2))./locpredxc0;
   
   loconexc=squeeze(res.onexc(2,:,:))';
   loconexcdiff=loconexc(:,1)-loconexc(:,2);
   
   % cells with significant baseline predictions
   catidx=squeeze((res.predxc(1,1,:)<res.randxc.*2)+1);
   %catidx=squeeze((res.predxc(1,1,:)<0.15)+1);
   
   % cells with significant local modulation 
   catidx2=ones(size(localpredxc));
   catidx2(find(res.pxc(2,:)>=0.05))=2;
   c2idx=find([catidx2;catidx2]==1 & strfbaserank>0);
   
   figure(1);
   clf
   
   subplot(2,3,1);
   corrcomp(obsbaserank,strfbaserank,'act','pred base targ rnk',...
            [],[catidx;catidx]);
   
   % subst dcpredxc for abs(strfmeandiff)
   subplot(2,3,2);
   %corrcomp(obsattrank,strfattrank,'act','pred base targ rnk',...
   %         [],[catidx;catidx]);
   corrcomp((obsattrankdiff),(strfmeandiff),...
            'act att rank diff','dc shift',[],catidx);
   
   subplot(2,3,3);
   corrcomp((obsbaserankdiff),(strfmeandiff),...
            'act base rank diff','dc shift',[],catidx);
   
   subplot(2,3,4);
   corrcomp((strfbaserankdiff),(strfmeandiff),...
            'pred rank diff','dc shift',[],catidx);
   
   subplot(2,3,5);
   corrcomp(strfbaserankmin,(strfmeandiff).*strfbaserankminsign,...
            'min pred att','dc shift',[],catidx);
   
   subplot(2,3,6);
   corrcomp(strfbaserankmean,abs(strfmeandiff),...
            'mean pred att','dc shift',[],catidx);
   
   figure(2);
   clf
   
   subplot(3,3,1);
   corrcomp(strfbaserankmin,localpredxc,...
            'min pred base rank','local att cc',[],catidx);
   
   subplot(3,3,2);
   corrcomp((strfmeandiff),localpredxc,...
            'abs dc shift','local att cc',[],catidx);
   
   subplot(3,3,3);
   corrcomp((obsattrankmin),localpredxc,...
            'min obs att rank','local att cc',[],catidx);
   
   subplot(3,3,4);
   corrcomp(abs(obsbaserankdiff),basepredxc,...
            'obs rank diff','base cc',[],catidx);
   
   subplot(3,3,5);
   corrcomp(obsbaserankmin,basepredxc,...
            'min obs att rank','base cc',[],catidx);
   
   subplot(3,3,6);
   corrcomp(obsbaserankmean,basepredxc,...
            'mean obs att rank','base cc',[],catidx);
   
   subplot(3,3,7);
   corrcomp(abs(strfattrankdiff),localpredxc,...
            'pred rank diff','local att cc',[],catidx);
   
   subplot(3,3,8);
   corrcomp(strfattrankmin,localpredxc,...
            'min pred att rank','local att cc',[],catidx2);
   
   subplot(3,3,9);
   corrcomp(strfattrankmean,localpredxc,...
            'mean pred att rank','local att cc',[],catidx2);
   
   figure(3);
   clf
   
   subplot(2,3,1);
   % A and B combined... like ben's match analysis.
   pp=[squeeze(res.rank(1,3,:))+squeeze(res.rank(2,2,:))];
   rr=[squeeze(res.rank(1,2,:)-res.rank(1,3,:))+...
       squeeze(res.rank(2,3,:)-res.rank(2,2,:))];
   corrcomp(pp,rr,'base rank','att rank inc',[],catidx);
  
   subplot(2,3,2);
   % A and B sep.
   pp=[squeeze(res.rank(1,3,:)); squeeze(res.rank(2,2,:))];
   rr=[squeeze(res.rank(1,2,:)-res.rank(1,3,:));
       squeeze(res.rank(2,3,:)-res.rank(2,2,:))];
   corrcomp(pp,rr,'base rank','att rank inc',[],[catidx;catidx]);
   
   subplot(2,3,3);
   pp=[squeeze(res.strfrank(1,3,:))+ squeeze(res.strfrank(2,2,:))];
   rr=[squeeze(res.strfrank(1,2,:))-squeeze(res.strfrank(1,3,:))+...
       squeeze(res.strfrank(2,3,:))-squeeze(res.strfrank(2,2,:))];
   corrcomp(pp,rr,'pred rank','pred att rank inc',[],catidx);
   %subplot(2,3,3);
   %pp=[squeeze(res.strfrank(1,3,:)); squeeze(res.strfrank(2,2,:))];
   %rr=[squeeze(res.strfrank(1,2,:))-squeeze(res.strfrank(1,3,:));
   %    squeeze(res.strfrank(2,3,:))-squeeze(res.strfrank(2,2,:))];
   %corrcomp(pp,rr,'pred rank','pred att rank inc',[],[catidx;catidx]);
   
   subplot(2,3,4);
   % A and B combined... like ben's match analysis.
   pp=[squeeze(res.resp(1,3,:))+squeeze(res.resp(2,2,:))];
   rr=[squeeze(res.resp(1,2,:)-res.resp(1,3,:))+...
       squeeze(res.resp(2,3,:)-res.resp(2,2,:))];
   corrcomp(pp,rr,'base resp','att resp inc',[],[catidx]);
  
   subplot(2,3,5);
   % A and B sep.
   %pp=[squeeze(res.resp(1,3,:)); squeeze(res.resp(2,2,:))];
   %rr=[squeeze(res.resp(1,2,:)-res.resp(1,3,:));...
   %    squeeze(res.resp(2,3,:)-res.resp(2,2,:))];
   %corrcomp(pp,rr,'base resp','att resp inc',[],[catidx;catidx]);
   
   %pp=basepredxc;
   pp=localpredxc;
   %rr=[squeeze(res.strfresp(1,2,:)-res.strfresp(1,3,:))+...
   %    squeeze(res.strfresp(2,3,:)-res.strfresp(2,2,:))];
   %rr=[squeeze(res.resp(1,2,:)-res.resp(1,3,:))+...
   %    squeeze(res.resp(2,3,:)-res.resp(2,2,:))];
   rr=[squeeze(res.rank(1,2,:)-res.rank(1,3,:))+...
       squeeze(res.rank(2,3,:)-res.rank(2,2,:))];
   %rr=[squeeze(res.strfrank(1,2,:)-res.strfrank(1,3,:))+...
   %    squeeze(res.strfrank(2,3,:)-res.strfrank(2,2,:))];
   corrcomp(pp,rr,'localpredxc','att rank inc',[],[catidx]);
   
   
   subplot(2,3,6);
   pp=[squeeze(res.strfresp(1,3,:))+squeeze(res.strfresp(2,2,:))];
   rr=[squeeze(res.strfresp(1,2,:)-res.strfresp(1,3,:))+...
       squeeze(res.strfresp(2,3,:)-res.strfresp(2,2,:))];
   corrcomp(pp,rr,'base strfresp','att strfresp inc',[],catidx);
   %pp=[squeeze(res.strfresp(1,3,:)); squeeze(res.strfresp(2,2,:))];
   %rr=[squeeze(res.strfresp(1,2,:)-res.strfresp(1,3,:)); ...
   %    squeeze(res.strfresp(2,3,:)-res.strfresp(2,2,:))];
   %corrcomp(pp,rr,'base strfresp','att strfresp inc',[],[catidx;catidx]);
   
   figure(4);
   clf
   
   subplot(2,3,1);
   xx=sort([sort(baseonexc')']);
   plot(xx(:,1),'r');
   hold on
   plot(xx(:,2),'b');
   hold off
   
   subplot(2,3,2);
   xx=loconexc;
   xx(find(strfattrankdiff>0),:)=flipdim(xx(find(strfattrankdiff>0),:),2);
   
   scatter(xx(:,1),xx(:,2),'k.');
   hold on
   plot([0 0.7],[0 0.7],'k--');
   hold off
   
   subplot(2,3,3);
   %corrcomp([abs(strfattrankdiff);-abs(strfattrankdiff)],...
   %         [xx(:,2)-xx(:,1);xx(:,1)-xx(:,2)],...
   %         'strfattrankdiff','baseonexcdiff',[],[catidx;catidx]);
   corrcomp(strfattrankdiff,loconexc(:,1)-loconexc(:,2),...
            'strfattrankdiff','loconexcdiff',[],catidx);
   
   subplot(2,3,4);
   corrcomp(strfmeandiff,loconexc(:,1)-loconexc(:,2),...
            'strfmeandiff','loconexcdiff',[],catidx);
   
   subplot(2,3,5);
   corrcomp(strfattrankdiff,strfmeandiff,...
            'strfattrankdiff','strfmeandiff',[],catidx);
   
   subplot(2,3,6);
   corrcomp(strfattrankdiff,loconexcdiff,...
            'strfattrankdiff','loconexcdiff',[],catidx);
   
   
   
   disp('paused after display ');
   keyboard
   
else % ie, strcmp(batchdata.matcmd,'kernvsatt2')
   % fvvs summary stuff
   
   % predicted resp to targets
   basestrfresp=squeeze(res.strfresp(:,1,:))';
   basestrfrespnodc=squeeze(res.strfrespnodc(:,1,:))';
   basestrfresppos=squeeze(res.strfresppos(:,1,:))';
   basestrfrespneg=squeeze(res.strfrespneg(:,1,:))';
   attuse=size(basestrfresp,2);
   attstrfresp=zeros(size(basestrfresp));
   attstrfrespnodc=zeros(size(basestrfresp));
   attstrfresppos=zeros(size(basestrfresp));
   attstrfrespneg=zeros(size(basestrfresp));
   attstrfresppos0=zeros(size(basestrfresp));
   attstrfrespneg0=zeros(size(basestrfresp));
   outstrfresp=zeros(size(basestrfresp));
   for attidx=1:attuse,
      attstrfresp(:,attidx)=squeeze(res.strfresp(attidx,attidx+1,:));
      attstrfrespnodc(:,attidx)=squeeze(res.strfrespnodc(attidx,attidx+1,:));
      attstrfresppos(:,attidx)=squeeze(res.strfresppos(attidx,attidx+1,:));
      attstrfrespneg(:,attidx)=squeeze(res.strfrespneg(attidx,attidx+1,:));
      att2=attidx+2;
      if att2>attuse,
         att2=2;
      end
      attstrfresppos0(:,attidx)=squeeze(res.strfresppos(attidx,att2,:));
      attstrfrespneg0(:,attidx)=squeeze(res.strfrespneg(attidx,att2,:));
      attstrfrespneg(:,attidx)=squeeze(res.strfrespneg(attidx,attidx+1,:));
      outstrfresp(:,attidx)=...
          squeeze(mean(res.strfresp(attidx,[2:attidx attidx+2:end],:)));
   end
   
   % predicted rank of targets among all stim in dataset
   basestrfrank=squeeze(res.strfrank(:,1,:))';
   % diff betw single att cond and mean of other att conds
   basestrfrankdiff=zeros(size(basestrfrank));
   % also min/mean rank of targets across all conds
   basestrfrankmin=min(basestrfrank,[],2);
   basestrfrankmax=max(basestrfrank,[],2);
   basestrfrankmean=mean(basestrfrank,2);
   basestrfrankstd=std(basestrfrank,0,2);
   basestrfrankrange=basestrfrankmax-basestrfrankmin;
   
   for attidx=1:attuse,
      rm=mean(basestrfrank(:,[1:attidx-1 attidx+1:end]),2);
      basestrfrankdiff(:,attidx)=basestrfrank(:,attidx)-rm;
   end
   
   % same things for att kernels
   attstrfrank=zeros(size(basestrfrank));
   outstrfrank=zeros(size(basestrfrank));
   attstrfrankdiff=zeros(size(attstrfrank));
   for attidx=1:attuse,
      attstrfrank(:,attidx)=squeeze(res.strfrank(attidx,attidx+1,:));
      outstrfrank(:,attidx)=...
          squeeze(mean(res.strfrank(attidx,[2:attidx attidx+2:end],:)));
      
      rm=squeeze(mean(res.strfrank([1:attidx-1 attidx+1:end],attidx+1,:)));
      attstrfrankdiff(:,attidx)=attstrfrankdiff(:,attidx)-rm;
   end
   attstrfrankmean=mean(attstrfrank,2);
   attstrfrankmin=min(attstrfrank,[],2);
   attstrfrankmax=max(attstrfrank,[],2);
   
   % difference in mean fireing rate between att cond and avg
   % firing rate over other conds
   dcbase=squeeze(res.mresp(1,1,:));
   dc=squeeze(res.mresp(2:end,:,:))';
   dcerr=squeeze(res.dcerr)';
   
   % normalized by mean att-out rate
   dcdiff=zeros(size(dc));
   dcdiffnorm=zeros(size(dc));
   dcnorm=zeros(size(dc));
   
   for attidx=1:attcount-1,
      dm=mean(dc(:,[1:attidx-1 attidx+1:end]),2);
      dcdiff(:,attidx)=dc(:,attidx)-dm;
      dcdiffnorm(:,attidx)=dcdiff(:,attidx)./(dm+(dm==0));
      
      dm=mean(dc,2);
      dcnorm(:,attidx)=(dc(:,attidx)-dm)./dm;
   end
   
   % cells with significant baseline predictions
   %catidx=squeeze((res.predxc(1,1,:)<res.randxc)+1);
   catidx=squeeze((res.predxc(1,1,:)<res.randxc.*2)+1);
   %catidx=squeeze((res.predxc(1,1,:)<0.1)+1);
   
   baseonexc=squeeze(res.predxccross(1,2:end,:))';
   loconexc=squeeze(res.predxccross(end,2:end,:))';
   loconexc0=squeeze(res.predxccrossrand(end,2:end,:))';
   basepredxc=squeeze(res.predxc(1,1,:));
   locpredxc=squeeze(res.predxc(1,end,:));
   locpredxcrand=squeeze(res.predxc(2,end,:));
   localxcdiff=(locpredxc.^2-locpredxcrand.^2);
   localpredp=squeeze(res.pxc(end,1,:))';
   dcpredxc=squeeze(res.predxc(1,2,:));
   dcpredxcrand=squeeze(res.predxc(2,2,:));
   dcpredp=squeeze(res.pxc(2,1,:))';
   
   locpreddiff=locpredxc.*abs(locpredxc)-locpredxcrand.*abs(locpredxcrand);
   locpreddiff=sqrt(abs(locpreddiff)).*sign(locpreddiff);
   
   dcpreddiff=dcpredxc.*abs(dcpredxc)-dcpredxcrand.*abs(dcpredxcrand);
   dcpreddiff=sqrt(abs(dcpreddiff)).*sign(dcpreddiff);
   
   PTHRESH=0.05;
   
   catidx2=(dcpredp(:)>PTHRESH)+1;
   catidx2(isnan(dcpredp))=0;
   
   catidx3=(localpredp(:)>PTHRESH)+1;
   catidx3(isnan(localpredp))=0;
   
   catidx4=(localpredp(:)>PTHRESH & dcpredp(:)>PTHRESH)+1;
   catidx4(isnan(localpredp))=0;
   
   catidx5=(loconexc(:)<loconexc0(:))+1;
   catidx5(isnan(loconexc(:)))=0;
   
   cellcount=length(res.cellid);
   dccorrin=zeros(cellcount,1);
   predcorrin=zeros(cellcount,1);
   dccorr=zeros(cellcount,1);
   predcorr=zeros(cellcount,1);
   prefcorr=zeros(cellcount,1);
   matchtest=zeros(cellcount,1);
   matchtest4=zeros(cellcount,4).*nan;
   %for pidx=1:paircount,
   %   p1=pairidx(pidx,1);
   %   p2=pairidx(pidx,2);
   %   bstrfrank(:,pidx)=basestrfrank(:,p1)-basestrfrank(:,p2);
   %end
   
   % ts is target similarity between each pair of targets,
   % measured as xcorr in the linearized domain
   ts=res.targsim.*abs(res.targsim);
   %ts=res.targsim;
   for ii=1:cellcount,
      dccorr(ii)=xcov(ts(:,ii),res.dcsim(:,ii),0,'coeff');
      %dccorr(ii)=xcov(res.prefsim(:,ii),res.dcsim(:,ii),0,'coeff');
      predcorr(ii)=xcov(ts(:,ii),res.predsim(:,ii),0,'coeff');
      prefcorr(ii)=xcov(res.prefsim(:,ii).*abs(res.prefsim(:,ii)),...
                        res.predsim(:,ii),0,'coeff');
      
      dccorrin(ii)=xcov(basestrfresp(ii,:)',res.mresp(2:end,ii),0,'coeff');
      predcorrin(ii)=xcov(basestrfresp(ii,:)',...
                          attstrfresp(ii,:)'-basestrfresp(ii,:)',...
                          0,'coeff');
      %predcorrin(ii)=xcov(ts(:,ii),bstrfrank(ii,:)',0,'coeff');
      tt=find(abs(res.predin(:,ii)+res.predout(:,ii))>0);
      if length(tt)>0,
         matchtest(ii)=median((res.predin(tt,ii)-res.predout(tt,ii))./...
                            (res.predin(tt,ii)+res.predout(tt,ii)));
         matchtest4(ii,tt)=((res.predin(tt,ii)-res.predout(tt,ii))./...
                            (res.predin(tt,ii)+res.predout(tt,ii)))';
      end
   end
   histcenters=-7/8:0.25:7/8;
   
   % [1] to include only sig baseline cells, [1 2] any cells
   catset=ismember(catidx,[1 2 ]);
   predok=ones(size(squeeze(res.randxc)));
   %predok=squeeze(res.predxc(1,1,:)>res.randxc);
   modset=(catidx3==1 & predok);
   modset2=(catidx4==1 & predok);
   modset3=((catidx2==1 | catidx3==1) & predok);
   
   % change catidx>=1 to catidx==1 to include only sig baseline
   % cells
   
   figure(1);
   clf
   subplot(2,2,1);
   %[n,x]=hist(prefcorr(modset & catset),histcenters);
   histcenters2=histcenters./6;
   matchtest(find(matchtest>max(histcenters2)))=max(histcenters2);
   matchtest(find(matchtest<min(histcenters2)))=min(histcenters2);
   ttd=matchtest(modset & catset);
   [n,x]=hist(ttd(:),histcenters2);
   bar(x,n);
   axis square
   %title(sprintf('prefsim vs modsim(loc)=%.2f/%.2f',...
   %              nanmean(prefcorr(modset & catset)),...
   %              nanmedian(prefcorr(modset & catset))));
   title(sprintf('matchtest(loc)=%.2f/%.2f',...
                 nanmean(matchtest(modset & catset)),...
                 nanmedian(matchtest(modset & catset))));
   
   subplot(2,2,2);
   [n,x]=hist(predcorr(modset & catset),histcenters);
   bar(x,n);
   axis([-1 1 0 6]);
   axis square
   title(sprintf('targsim v modsim(loc)=%.2f (mn)/%.2f (med)',...
                 mean(predcorr(modset & catset)),...
                 median(predcorr(modset & catset))));
   
   subplot(2,2,3);
   %catidx=squeeze((res.predxc(1,5,:)<res.randxc.*1)+1);
   ncount=squeeze(res.ncount)';
   catidx=squeeze((ncount(:,1)<2500)+1);
   catset2=ismember(catidx,[1 ]);
   %histcomp(matchtest(modset & catset),[],'targsim','','attidx',[-1 1 0 6]);
   histcomp(matchtest(modset & catset2),...
            matchtest(modset & ~catset2 & catset),'targsim','','attidx',...
            [-0.15 0.15 0 10]);
   %[n,x]=hist(predcorrin(modset & catset),histcenters);
   %bar(x,n);
   %title(sprintf('predcorrin(loc)=%.2f/%.2f',...
   %              mean(predcorrin(modset & catset)),...
   %              median(predcorrin(modset & catset))));
   
   subplot(2,2,4);
   %catidx=squeeze((res.predxc(1,5,:)<res.randxc.*1)+1);
   ncount=squeeze(res.ncount)';
   catidx=squeeze((ncount(:,1)<2500)+1);
   catset2=ismember(catidx,[1 ]);
   %histcomp(predcorr(modset & catset),[],'targsim','','attidx',[-1 1 0 6]);
   histcomp(predcorr(modset & catset2),...
            predcorr(modset & ~catset2 & catset),'targsim','','attidx',...
            [-1 1 0 6]);
   
   %[n,x]=hist(dccorr(catidx4==1 & catset),histcenters);
   %bar(x,n);
   %title(sprintf('dccorr(dc)=%.2f/%.2f',...
   %              nanmean(dccorr(modset2 & catset)),...
   %              nanmedian(dccorr(modset2 & catset))));
   
   % dump details for local cells
   [crap,ttidx]=sortrows([catidx3 predcorr]);
   for ii=1:cellcount,
      ttt=ttidx(ii);
      if catidx3(ttt)==1,
         fprintf('%6s %3d %2d %.3f %.3f %.3f %.3f\n',...
                 res.cellid{ttt},ttt,catidx4(ttt),predcorr(ttt),...
                 prefcorr(ttt),dccorr(ttt),localxcdiff(ttt));
      end
   end
   set(gcf,'PaperPosition',[2 2 4 4]);
   
   colormap(gray);
   
   keyboard
   
   if 0,
      
      ts=(res.targsim(:,find(catidx3==1)));
      ts=ts.*abs(ts);
      ps=(res.predsim(:,find(catidx3==1)));
      figure(2);
      clf
      scatter(ts(:)+ps(:),ts(:)-ps(:))
      [xcov(ts(:),ps(:),0,'coeff')]
      
      cc=zeros(400,1);
      for ccidx=1:length(cc),
         rr=rand(100,4);
         pp=rand(100,4);
         
         pairidx=[1 2; 1 3; 1 4; 2 3; 2 4 ; 3 4];
         paircount=6;
         tt=zeros(paircount,1);
         tt2=zeros(paircount,1);
         for pidx=1:paircount,             
            p1=pairidx(pidx,1);         
            p2=pairidx(pidx,2);         
            tt(pidx)=xcorr(pp(:,p1),pp(:,p2),0,'coeff');
            tt2(pidx)=xcorr(rr(:,p1),rr(:,p2),0,'coeff');
         end
         cc(ccidx)=xcov(tt,tt2,0,'coeff');
      end
   end
   
   figure(1);
   clf
   
   ranksign=attstrfrank-0.05>basestrfrank;
   ranksign=attstrfrank+0.05<basestrfrank;
   bar([hist(sum(ranksign(catidx3==1,:)'),5);
        hist(sum(ranksign(catidx3==2,:)'),5)]');
   
   corrcomp(dc(:),...
            (attstrfresp(:)-basestrfresp(:)),...
            'baserank','attrank',[],repmat(catidx2,[4 1]));
   plotcomp(nanmean(basestrfrank')',nanmean(attstrfrank')',...
            'baserankdiff','attrankdiff',[],catidx3);
    
   
   
   subplot(3,3,1);
   corrcomp(attstrfrankmin,locpreddiff,...
            'min att rank','local att cc',[],catidx);
   
   subplot(3,3,2);
   corrcomp(attstrfrankmean,locpreddiff,...
            'mean att rank','local att cc',[],catidx);
   
   subplot(3,3,3);
   corrcomp(attstrfrankmax,locpreddiff,...
            'max att rank','local att cc',[],catidx);
   
   subplot(3,3,4);
   corrcomp(attstrfrankmin,dcxcdiff,...
            'min att rank','dc att cc',[],catidx);
   
   subplot(3,3,5);
   corrcomp(attstrfrankmean,dcxcdiff,...
            'mean att rank','dc att cc',[],catidx);
   
   subplot(3,3,6);
   corrcomp(attstrfrankmax,dcxcdiff,...
            'max att rank','dc att cc',[],catidx);
   
   subplot(3,3,4);
   corrcomp(attstrfrankmax-attstrfrankmin,dcxcdiff,...
            'obs rank diff','att cc',[],catidx);
   
   subplot(3,3,5);
   corrcomp(basestrfrankmin,basepredxc,...
            'min obs att rank','base cc',[],catidx);
   
   subplot(3,3,6);
   corrcomp(nanmax(dcnorm')'-nanmin(dcnorm')',locpreddiff,...
            'max dc shift','local att cc',[],catidx);
   
   subplot(3,3,7);
   corrcomp(attstrfrankmax-attstrfrankmin,locpreddiff,...
            'pred rank diff','local att cc',[],catidx);
   
   subplot(3,3,8);
   corrcomp(attstrfrankmin,locpreddiff,...
            'min pred att rank','local att cc',[],catidx);
   
   subplot(3,3,9);
   corrcomp(attstrfrankmean,locpreddiff,...
            'mean pred att rank','local att cc',[],catidx);
   
   
   
   
   
   
   
   figure(1);
   clf
   disp('processing rank vs predxc');
   
   subplot(2,3,1);
   xx=sort([sort(loconexc(catidx==1,:)')']);
   plot(xx);
   title('local xc sorted');
   
   xx=zeros(size(loconexc,1),2).*nan;
   bxx=zeros(size(baseonexc,1),2).*nan;
   brank=zeros(size(loconexc,1),2).*nan;
   arank=zeros(size(loconexc,1),2).*nan;
   dcextreme=zeros(size(loconexc,1),2).*nan;
   for ii=1:length(xx),
      
      % could use basestrfrank here, but attstrfrank appears better
      % could use baseonexc here, but loconexc appears better
      [ss,ssi]=sort(attstrfrank(ii,:));
      %[ss,ssi]=sort(basestrfrank(ii,:));
      
      if ~isnan(ss(ssi(2))),
         aa1=ssi(1);
         aa2=ssi(max(find(~isnan(ss))));
         %xx(ii,:)=loconexc(ii,[aa1 aa2]);
         xx(ii,:)=baseonexc(ii,[aa1 aa2]);
         brank(ii,:)=[basestrfrank(ii,[aa1 aa2])];
         arank(ii,:)=[attstrfrank(ii,[aa1 aa2])];
         %srank(ii,:)=[attstrfresp(ii,[aa1 aa2])];
         %srank(ii,:)=[basestrfresp(ii,[aa1 aa2])];
         %dcextreme(ii,:)=[dc(ii,[aa1 aa2])];
         dcextreme(ii,:)=[dcnorm(ii,[aa1 aa2])];
      end
   end
   ccname='baseonexcdiff';
   rankname='attstrfrankdiff';
   dcname='dcdiffnorm';
   
   subplot(2,3,2);
   plotcomp(xx(:,1),xx(:,2),'badtargcc','goodtargcc',[0 0.7 0 0.7],catidx);
   
   subplot(2,3,3);
   corrcomp((arank(:,2)-arank(:,1)),xx(:,2)-xx(:,1),...
            rankname,ccname,[],catidx);
   
   subplot(2,3,4);
   corrcomp(dcextreme(:,2)-dcextreme(:,1),xx(:,2)-xx(:,1),...
            dcname,ccname,[],catidx);
   
   subplot(2,3,5);
   corrcomp(arank(:,2)-arank(:,1),dcextreme(:,2)-dcextreme(:,1),...
            rankname,dcname,[],catidx);
   
   %subplot(2,3,6);
   %corrcomp((dcextreme(:).*1000),xx(:),...
   %         'dcone','xcone');
   
   subplot(2,3,6);
   corrcomp(brank(:,2)-brank(:,1),...
            (arank(:,2)-brank(:,2))-(arank(:,1)-brank(:,1)),...
            'baserankdiff','attrankdiff',[],catidx);
   
    
   figure(3);
   clf
      
   aa=(attstrfresp-repmat(dcbase,[1 attcount-1]))./...
      repmat(dcbase,[1 attcount-1]);
   bb=(basestrfresp-repmat(dcbase,[1 attcount-1]))./...
      repmat(dcbase,[1 attcount-1]);
   
   aa=attstrfrank;
   bb=basestrfrank;
   
   GOODTHRESH=nanmean(nanmean(aa(catidx3>0,:)));
   GOODTHRESH=0.00;
   GOODTHRESH=0.55;
   
   ccidx=zeros(size(catidx));
   ccidx(catidx3==1)=1;
   ccidx(catidx2==1 & catidx3~=1)=2;
   ccidx(catidx2~=1 & catidx3~=1)=3;
   %ccidx(ccidx==3 & catidx==2)=0;
   ccidx(catidx==2)=0;
   
   cc=[mean(sum(aa(ccidx==1,:)>GOODTHRESH,2)) ...
       std(sum(aa(ccidx==1,:)>GOODTHRESH,2))./sqrt(sum(ccidx==1)); ...
       mean(sum(aa(ccidx==2,:)>GOODTHRESH,2)) ...
       std(sum(aa(ccidx==2,:)>GOODTHRESH,2))./sqrt(sum(ccidx==2)); ...
       mean(sum(aa(ccidx==2,:)>GOODTHRESH,2)) ...
       std(sum(aa(ccidx==2,:)>GOODTHRESH,2))./sqrt(sum(ccidx==2));
       mean(sum(bb(ccidx==1,:)>GOODTHRESH,2)) ...
       std(sum(bb(ccidx==1,:)>GOODTHRESH,2))./sqrt(sum(ccidx==1)); ...
       mean(sum(bb(ccidx==2,:)>GOODTHRESH,2)) ...
       std(sum(bb(ccidx==2,:)>GOODTHRESH,2))./sqrt(sum(ccidx==2)); ...
       mean(sum(bb(ccidx==2,:)>GOODTHRESH,2)) ...
       std(sum(bb(ccidx==2,:)>GOODTHRESH,2))./sqrt(sum(ccidx==2))];
   dd=[nanmean(nanmin(aa(ccidx==1,:)')) ...
       nanstd(nanmin(aa(ccidx==1,:)'))./sqrt(sum(ccidx==1)); ...
       nanmean(nanmin(aa(ccidx==2,:)')) ...
       nanstd(nanmin(aa(ccidx==2,:)'))./sqrt(sum(ccidx==2));
       nanmean(nanmin(aa(ccidx==3,:)')) ...
       nanstd(nanmin(aa(ccidx==3,:)'))./sqrt(sum(ccidx==3));
       nanmean(nanmin(bb(ccidx==1,:)')) ...
       nanstd(nanmin(bb(ccidx==1,:)'))./sqrt(sum(ccidx==1)); ...
       nanmean(nanmin(bb(ccidx==2,:)')) ...
       nanstd(nanmin(bb(ccidx==2,:)'))./sqrt(sum(ccidx==2));
       nanmean(nanmin(bb(ccidx==3,:)')) ...
       nanstd(nanmin(bb(ccidx==3,:)'))./sqrt(sum(ccidx==3))];
   ee=[nanmean(nanmean(aa(ccidx==1,:)')) ...
       nanstd(nanmean(aa(ccidx==1,:)'))./sqrt(sum(ccidx==1)); ...
       nanmean(nanmean(aa(ccidx==2,:)')) ...
       nanstd(nanmean(aa(ccidx==2,:)'))./sqrt(sum(ccidx==2));
       nanmean(nanmean(aa(ccidx>1,:)')) ...
       nanstd(nanmean(aa(ccidx>1,:)'))./sqrt(sum(ccidx>1));
       nanmean(nanmean(bb(ccidx==1,:)')) ...
       nanstd(nanmean(bb(ccidx==1,:)'))./sqrt(sum(ccidx==1)); ...
       nanmean(nanmean(bb(ccidx==2,:)')) ...
       nanstd(nanmean(bb(ccidx==2,:)'))./sqrt(sum(ccidx==2));
       nanmean(nanmean(bb(ccidx>1,:)')) ...
       nanstd(nanmean(bb(ccidx>1,:)'))./sqrt(sum(ccidx>1))];
   
   subplot(2,3,1);
   errorbar(cc(:,1),cc(:,2),'k+');
   hold on
   bar(cc(:,1));
   hold off
   title('# targets over mean');
   
   subplot(2,3,2);
   errorbar(dd(:,1),dd(:,2),'k+');
   hold on
   bar(dd(:,1));
   hold off
   title('min target');
   
   subplot(2,3,3);
   errorbar(ee(:,1),ee(:,2),'k+');
   hold on
   bar(ee(:,1));
   hold off
   title('mean target');
   
   aa=attstrfrankmax-attstrfrankmin;
   drange=sort(aa(catidx3>0));
   ll=round(linspace(1,length(drange)+1,4));
   brange2=[drange(ll(1:end-1)); inf];
   locatt=histc(aa(catidx3==1),brange2);
   locnoatt=histc(aa(catidx3==2),brange2);
   locatt=locatt(1:end-1); locnoatt=locnoatt(1:end-1);
   [brange2(1:end-1) locatt locnoatt]
   locatt=locatt./(locatt+locnoatt);
   
   dcatt=histc(aa(catidx2==1),brange2);
   dcnoatt=histc(aa(catidx2==2),brange2);
   dcatt=dcatt(1:end-1); dcnoatt=dcnoatt(1:end-1);
   dcatt=dcatt./(dcatt+dcnoatt);
   
   subplot(2,3,4);
   bar(brange2(1:end-1),[locatt dcatt]);
   
   drange=sort(nanmean(aa(catidx3>0,:)'));
   ll=round(linspace(1,length(drange)+1,5));
   brange2=[drange(ll(1:end-1)) inf];
   
   locatt=histc(nanmean(aa(catidx3==1,:)'),brange2);
   locnoatt=histc(nanmean(aa(catidx3==2,:)'),brange2);
   [locatt; locnoatt]
   locatt=locatt./(locatt+locnoatt);
   
   dcatt=histc(nanmean(aa(catidx2==1,:)'),brange2);
   dcnoatt=histc(nanmean(aa(catidx2==2,:)'),brange2);
   dcatt=dcatt./(dcatt+dcnoatt);
   
   subplot(2,2,3);
   bar(brange2,[locatt; dcatt]');
   
   drange=sort(nanmin(aa(catidx3>0,:)'));
   ll=round(linspace(1,length(drange)+1,5));
   brange2=[drange(ll(1:end-1)) inf];
   
   locatt=histc(nanmin(aa(catidx3==1,:)'),brange2);
   locnoatt=histc(nanmin(aa(catidx3==2,:)'),brange2);
   [locatt; locnoatt]
   locatt=locatt./(locatt+locnoatt);
   
   dcatt=histc(nanmin(aa(catidx2==1,:)'),brange2);
   dcnoatt=histc(nanmin(aa(catidx2==2,:)'),brange2);
   dcatt=dcatt./(dcatt+dcnoatt);
   
   subplot(2,2,3);
   bar(brange2(1:end-1),[locatt(1:end-1); dcatt(1:end-1)]');
   
   drange=sort(nanmax(aa(catidx3>0,:)')-nanmin(aa(catidx3>0,:)'));
   ll=round(linspace(1,length(drange)+1,5));
   brange2=drange(ll(1:end-1));
   %brange2=0:4;
   locatt=hist(nanmax(aa(catidx3==1,:)')-nanmin(aa(catidx3==1,:)'),brange2);
   locnoatt=hist(nanmax(aa(catidx3==2,:)')-nanmin(aa(catidx3==2,:)'),brange2);
   [locatt; locnoatt]
   locatt=locatt./(locatt+locnoatt);
   
   dcatt=hist(nanmax(aa(catidx2==1,:)')-nanmin(aa(catidx2==1,:)'),brange2);
   dcnoatt=hist(nanmax(aa(catidx2==2,:)')-nanmin(aa(catidx2==2,:)'),brange2);
   [dcatt; dcnoatt]
   dcatt=dcatt./(dcatt+dcnoatt);
   
   subplot(2,2,4);
   bar(brange2,[locatt; dcatt]');
   
   
   
   disp('paused before finished to avoid crashing b/c buggy plot code');
   keyboard
   
   
   
   plotcomp(nanmean(basestrfrank')',nanmean(attstrfrank')',...
            'baserankdiff','attrankdiff',[],catidx3);
   
   
   plotcomp(brank(:,2)-brank(:,1),arank(:,2)-arank(:,1),...
            'baserankdiff','attrankdiff',[],catidx3);
   
   dcxcdiff=((dcpredxc.^2-dcpredxcrand.^2)./dcpredxc.^2);
   dcxcdiff(isnan(dcxcdiff))=0;
   localxcdiff=((locpredxc.^2-locpredxcrand.^2)./locpredxc.^2);
   localxcdiff(isnan(localxcdiff))=0;
   
   corrcomp(nanmedian(attstrfrank')',-dcxcdiff,...
            'baseminrank','localxcdiff',[],catidx2);
   corrcomp(brank(:,1),localxcdiff,...
            'baseminrank','dcxcdiff',[],catidx3);
   
   
   
   corrcomp(brank(:,1),localxcdiff,...
            'baserankdiff','localxcdiff',[],catidx3);
   corrcomp(min(basestrfrank,[],2),localxcdiff,...
            'attminrank','localxcdiff',[],catidx);
   corrcomp(mean(basestrfrank,2),localxcdiff,...
            'basemeanrank','localxcdiff',[],catidx);
   corrcomp(brank(:,2)-brank(:,1),dcxcdiff,...
            rankname,dcname,[],catidx);
   corrcomp(arank(:,2)-arank(:,1),max(dcnorm,[],2)-min(dcnorm,[],2),...
            rankname,dcname,[],catidx);
   corrcomp(arank(:,2)-arank(:,1),dcextreme(:,2)-dcextreme(:,1),...
            rankname,dcname,[],catidx2);
   corrcomp(brank(:,2)-brank(:,1),dcextreme(:,2)-dcextreme(:,1),...
            rankname,dcname,[],catidx2);
   
   
   
   subplot(2,2,1);
   aa=attstrfrank(:);
   bb=dc(:);
   gidx=find(~isnan(aa+bb) & ~isinf(bb));
   aa=aa(gidx);
   bb=bb(gidx);
   corrcomp(aa,bb,...
            'basestrfrank','dcerr',[]);
   
   subplot(2,2,2);
   cla
   trank=[];
   tdc=[];
   for ii=1:size(attstrfrank,1),
      tt=sortrows([attstrfrank(ii,:);dc(ii,:)]');
      %plot(tt(:,1),tt(:,2));
      trank=[trank tt(:,1)];
      tdc=[tdc tt(:,2)];
   end
   %hold off
   

   corrcomp(attstrfrank(:),dc(:),...
            'attstrfrank','dc',[],repmat(catidx2,4,1));
   hold on
    for ii=1:size(attstrfrank,1),
      tt=sortrows([attstrfrank(ii,:);dc(ii,:)]');
      plot(tt(:,1),tt(:,2));
   end
   hold off
   
   corrcomp(attstrfrank(:),dcnorm(:),...
            'attstrfrank','dc',[],repmat(catidx2,4,1));
   
   plotcomp(basestrfrank(:),attstrfrank(:),...
            'basestrfrank','attstrfrank',[],repmat(catidx3,4,1));
   plotcomp(basestrfrank(:),attstrfrank(:),...
            'basestrfrank','attstrfrank',[],catidx4);
   plotcomp(basestrfresp(:),attstrfresp(:),...
            'basestrfresp','attstrfresp',[],repmat(catidx3,4,1));
   corrcomp(basestrfrespnodc(:),attstrfrespnodc(:)-basestrfrespnodc(:),...
            'basestrfrespnodc','attstrfrespnodc',[],repmat(catidx3,4,1));
   
   corrcomp(basestrfresppos(:),attstrfresppos(:)-basestrfresppos(:),...
            'basestrfresppos','attstrfresppos',[],repmat(catidx3,4,1));
   corrcomp(basestrfrespneg(:),attstrfrespneg(:)-basestrfrespneg(:),...
            'basestrfrespneg','attstrfrespneg',[],repmat(catidx3,4,1));
   corrcomp(basestrfrespneg(:),attstrfrespneg0(:)-basestrfrespneg(:),...
            'basestrfrespneg','attstrfrespneg',[],repmat(catidx3,4,1));
   corrcomp(basestrfresppos(:),attstrfresppos0(:)-basestrfresppos(:),...
            'basestrfresppos','attstrfresppos0',[],repmat(catidx3,4,1));
   
   
   subplot(2,2,2);
   corrcomp(attstrfrank(:),dc(:),...
            'attstrfrank','dc',[],repmat(catidx2,4,1));
   corrcomp(basestrfrespnodc(:),attstrfrespnodc(:),...
            'basestrfnodc','attstrfnodcinc',[],repmat(catidx,4,1));
   corrcomp(basestrfresppos(:)+basestrfrespneg(:),attstrfresppos(:)-basestrfresppos(:)+attstrfrespneg(:)-basestrfrespneg(:),...
            'basestrfpos','attstrfposinc',[],repmat(catidx,4,1));
   
   
   corrcomp(attstrfresppos(:)-basestrfresppos(:),...
            attstrfrespneg(:)-basestrfrespneg(:),...
            'basestrfpos','basestrfneg',[],repmat(catidx,4,1));
   corrcomp(basestrfresppos(:),attstrfresppos(:)-basestrfresppos(:),...
            'basestrfpos','attstrfposinc',[],repmat(catidx,4,1));
   corrcomp(basestrfrespneg(:),attstrfrespneg(:)-basestrfrespneg(:),...
            'basestrfneg','attstrfneginc',[],repmat(catidx,4,1));
   corrcomp(basestrfresppos(:)+basestrfrespneg(:),basestrfrespnodc(:),...
            'basestrfneg','attstrfneginc',[],repmat(catidx,4,1));
   corrcomp(attstrfresp(:),dc(:),...
            'attstrfresp','dc',[],repmat(catidx2,4,1));
   
   
   
   figure(2);
   
   % dcbase vs basepredxc is correlated
   corrcomp(log10(dcbase.*1000),basepredxc,...
            'log10(dcbase)','basepredxc',[],catidx);
   % mean dc vs dcpreddiff not correlated
   % std(dcdiff,0,2) vs dcpreddiff is correlated

   % basestrfrankmean vs dcpreddiff weakly correlated
   %  (better than basestrfrank-min, -max, -range or -std)
   corrcomp(basestrfrankmean,dcpreddiff,...
            'basestrfrankmean','dcpreddiff',[],catidx);
   
   %  basestrfrankmin vs std(dc,0,2) weakly correlated
   corrcomp(basestrfrankmin,std(dc,0,2),...
            'basestrfrankmin','std(dc,0,2)',[],catidx);
   
   % basestrfrankmin vs basepredxc weakly correlated
   corrcomp(basestrfrankmin,basepredxc,...
            'basestrfrankmin','basepredxc',[],catidx);
   
   % basestrfrank vs basexcone reasonably correlated
   corrcomp(basestrfrank(:),baseonexc(:),...
            'basestrfrank','baseonexc',[],repmat(catidx,attcount-1,1));
   
   % basestrfrank vs dcdiffnorm not correlated
   % (same for dcnorm, dcdiff)
   % (weak corr for dc)
   corrcomp(basestrfrank(:),dc(:),...
            'basestrfrank','dc',[],repmat(catidx2,4,1));
   corrcomp(attstrfrank(:),dc(:),...
            'attstrfrank','dc',[],repmat(catidx2,4,1));
   
   % weak but ns corr
   corrcomp(basepredxc,dcpreddiff,...
            'basepredxc','dcpreddiff',[],catidx);
    
   
   corrcomp(outstrfrank(:),(attstrfrank(:)-outstrfrank(:)),...
            'dcbase','basestrfrankmean',[],repmat(catidx,attcount-1,1));
  
   
   
   
   figure(4);
   clf
   
   
    
   
   figure(1);
   clf
   
   subplot(2,2,1);
   bgidx=find(~isnan(dcdiff(:)+basestrfrankdiff(:)));
   corrcomp(basestrfrankdiff(bgidx),dcdiff(bgidx),...
            'pred base rank diff','dc diff',[],cidx(bgidx));
   
   subplot(2,2,2);
   bgidx=find(~isnan(basestrfrank(:)+dcdiffnorm(:)));
   corrcomp(basestrfrank(bgidx),dcdiff(bgidx),...
            'pred base rank','dc diff',[],cidx(bgidx));
  
   subplot(2,2,3);
   bgidx=find(~isnan(dcdiffnorm(:)+basestrfrankdiff(:)));
   corrcomp(basestrfrankdiff(bgidx),dcdiffnorm(bgidx),...
            'pred base rank diff','dc diff',[],cidx(bgidx));

   subplot(2,2,4);
   bgidx=find(~isnan(dcdiffnorm(:)+attstrfrankdiff(:)));
   corrcomp(attstrfrankdiff(bgidx),dcdiffnorm(bgidx),...
            'pred att rank diff','dc diff',[],cidx(bgidx));
  
end

keyboard

