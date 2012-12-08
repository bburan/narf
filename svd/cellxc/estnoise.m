

if ~exist('batchlist','var'),
   batchlist=[51 53]
end
if ~exist('forcereload','var'),
   forcereload=0;
end

fprintf('%s.m: batch=',mfilename);
fprintf('%d ',batchlist);
fprintf('\n');

resfile=sprintf('/auto/k5/david/tmp/estnoise.mat');

if ~forcereload & exist(resfile,'file');
   RELOAD=0;
   
   % don't need to re-run all the high-level analyses
   fprintf('loading saved results: %s\n',resfile);
   load(resfile);
else
   RELOAD=1;
   % find all sRunData entries for this batch
   batchcount=length(batchlist);
   
   clear res
   for batchidx=1:length(batchlist),
      rundata=dbgetrundata(batchlist(batchidx));
      
      for runidx=1:length(rundata);

         fprintf('loading %s%s.gz\n',rundata(runidx).respath,rundata(runidx).resfile);
         z=zload([rundata(runidx).respath,rundata(runidx).resfile,'.gz']);
         
         strf=z.strf;
         
         global BATQUEUEID
         BATQUEUEID =[];
         
         %
         % VALIDATION CEILING
         %
         
         % figure out pred files for the current batch
         predidx=1;
         [pcellfiledata,ptimes,pbatchdata]=...
             cellfiletimes(z.params.cellid,z.params.predbatch{predidx},1);
         
         % does this cell have data for batchid=predidx?
         if length(pcellfiledata)>0,
            
            predparams=z.params;
            predparams.stimfiles=pbatchdata.stimfiles;
            predparams.respfiles=pbatchdata.respfiles;
            predparams.stimcrfs=pbatchdata.stimcrfs;
            
            predparams.stimloadcmd=pbatchdata.stimloadcmd;
            predparams.stimloadparms=strsep(pbatchdata.stimloadparms,',');
            predparams.stimfiltercmd=pbatchdata.stimfiltercmd;
            predparams.stimfilterparms=strsep(pbatchdata.stimfilterparms,',');
            predparams.resploadcmd=pbatchdata.resploadcmd;
            predparams.resploadparms=strsep(pbatchdata.resploadparms,',');
            predparams.respfiltercmd=pbatchdata.respfiltercmd;
            predparams.respfilterparms=strsep(pbatchdata.respfilterparms,',');
            predparams.times=ptimes;
            
            tpredstartframe=ptimes(3).start;
            tpredstopframe=ptimes(3).stop;
            tpredfile=ptimes(3).fileidx;
            predparams.resploadparms{4}=0;
            
            [cdata.stim,cdata.resp]=...
                xcloadstimresp(tpredfile,tpredstartframe,...
                               tpredstopframe,predparams);
            tresp=cdata.resp;
            
            if size(tresp,2)>1,
               tresp=tresp(:,2:end);
               [tresp]=compact_raster_matrix3(tresp);
               
               nancounts=sum(isnan(tresp),1);
               
               firstwithnans=min([find(nancounts>(nancounts(1).*2))...
                                  size(tresp,2)+1]);
               tresp=tresp(:,1:firstwithnans-1);
            end
            repcount=size(tresp,2);
            fprintf('predidx=%d (bat %d) (%d reps): ceiling... ',...
                    predidx,z.params.predbatch{predidx},repcount);
            
            gg=find(~isnan(tresp(:)));
            [mu,alpha,beta]=reversepoisson(tresp(gg));
            rmax=zeros(repcount,1);
            for xcidx=1:repcount,
               gg=find(~isnan(tresp(:,xcidx)));
               %[mu,alpha,beta]=reversepoisson(tresp(gg,xcidx));
               rmax(xcidx)=singletrialceiling(tresp(gg,xcidx),alpha,beta);
            end
            
            xct=zeros(repcount,z.bootcount,z.bcount-1);
            xctceil=zeros(repcount,z.bootcount,z.bcount-1);
            xctnormal=zeros(z.bootcount,z.bcount-1);
            for bootidx=1:z.bootcount,
               for segidx=1:z.bcount-1,
                  r=xcpredict(strf(end,bootidx,segidx),cdata.stim);
                  
                  mresp=nanmean(tresp')';
                  gg=find(~isnan(r) & ~isnan(mresp));
                  xctnormal(bootidx,segidx)=xcov(r(gg),mresp(gg),0,'coeff');
                  
                  for xcidx=1:repcount,
                     gg=find(~isnan(r) & ~isnan(tresp(:,xcidx)));
                     if var(r(gg))>0 & var(tresp(gg,xcidx))>0
                        xct(xcidx,bootidx,segidx)=...
                            xcov(r(gg),tresp(gg,xcidx),0,'coeff');
                        xctceil(xcidx,bootidx,segidx)=...
                            xct(xcidx,bootidx,segidx)./rmax(xcidx);
                     end
                  end
               end
            end
         else
            % no val data for this batch
         end
         
         %
         % ESTIMATION CEILING
         %
         
         % figure out estimation noise ceiling:
         %('(v)predxc is bfrac X nl ( X boot)')
         %strf is nl X boot X bfrac
         
         % xct is single trial
         tpxc1=xct;
         
         % xctceil is val noise corrected
         tpxc=xctceil;
         
         tpxc1(tpxc1<0.02)=0.02;
         tpxc(tpxc<0.02)=0.02;
         
         x=z.params.bfracs(1:size(tpxc,3)).*sum(~isnan(z.resp(:,1)));
         x=repmat(x,[repcount 1]);
         y=tpxc.^2;
         
         tp=zeros(z.bootcount,2);
         xcestinf=zeros(z.bootcount,1);
         for ii=1:z.bootcount,
            
            ty=squeeze(y(:,ii,:));
            [xcestinf(ii),errestinf,tp(ii,:)]=fitceiling(x,ty);
            
            [xcestinf(ii),errestinf,tp(ii,:)]=...
                fitceiling(x,mean(y([1:ii-1 ii+1:end],:)));
            %[xcestinf(ii),errestinf,tp(ii,:)]=fitceiling(x(:)',ty(:)');
         end
         mmgood=find(~isnan(xcestinf));
         mxcestinf=nanmean(xcestinf);
         
         figure(1);
         clf
         plot(x(1,:),squeeze(mean(mean(tpxc1(:,mmgood,:),1),2)),'k--');
         hold on
         plot(x(1,:),squeeze(mean(mean(tpxc(:,mmgood,:),1),2)),'k-');
         plot(max(x(:)).*1.2,sqrt(mxcestinf),'o');
         hold off
         drawnow
         title(sprintf('%s',z.params.cellid));
         
         res(runidx,batchidx).cellid=z.params.cellid;
         res(runidx,batchidx).runid=rundata(runidx).id;
         res(runidx,batchidx).resfile=...
             [rundata(runidx).respath,rundata(runidx).resfile];
         res(runidx,batchidx).xct=xct;
         res(runidx,batchidx).xctceil=xctceil;
         res(runidx,batchidx).rmax=rmax;
         
         res(runidx,batchidx).onetrialmean=...
             squeeze(mean(mean(tpxc1(:,mmgood,:),1),2));
         res(runidx,batchidx).valmaxmean=...
             squeeze(mean(mean(tpxc(:,mmgood,:),1),2));
         
         res(runidx,batchidx).xcestinf=xcestinf;
         res(runidx,batchidx).mmgood=mmgood;
         res(runidx,batchidx).mxcestinf=mxcestinf;
         
      end
   end
   
   %done with stuff that gets saved to the processed data
   %file. rest of script is mostly display
   
   % save res to file
   save(resfile);
end

disp('finished loading');
keyboard

return
         



cellcount=size(pp,1);
for ii=1:size(pp,1),
   [xcestinf,errestinf,tp]=fitceiling(x,pp(ii,:));
   pp2(ii,end)=xcestinf;
   
   if ~isnan(xcestinf),
      bb=mean(1./(tp(:,2)*ones(size(aa))+tp(:,1)*(1./aa)),1);
      ppfit(ii,:)=bb;
   end
   
   fprintf('%s (%d): %.3f-->%.3f\n',cellids{ii},ii,pp(ii,end),xcestinf);
end

gidx=find(~isnan(pp2(:,end)));
length(gidx)
mean(pp2(gidx,:))
sqrt(mean(pp2(gidx,:)))

%figure(1);
%clf
%plot(aa,ppfit,'--');
%hold on
%plot(bfracs,pp2,'o-');
%hold off

%yy=mean(pp2(:,1:end-1));
%[xcestinf,errestinf,tp]=fitceiling(x,yy);

unix('rm /auto/k1/david/tmp/res.dat');
for ii=1:cellcount,
   xcsepfullres(cellids{ii},51);
   drawnow
end
for ii=1:cellcount,
   xcsepfullres(cellids{ii},53);
   drawnow
end


r=load('/auto/k1/david/tmp/res.dat');
figure(1);
clf

imidx0=find(~isnan(r(:,end-1)) & ~isnan(r(:,end)) & r(:,end)<1.5 & r(:,2)==53);
pfidx0=find(~isnan(r(:,end-1)) & ~isnan(r(:,end)) & r(:,end)<1.5 & r(:,2)==51);

masters=intersect(r(imidx0,1),r(pfidx0,1));
imidx=imidx0(ismember(r(imidx0,1),masters));
pfidx=pfidx0(ismember(r(pfidx0,1),masters));

pcount=length(imidx);
subplot(1,2,1);
errorbar(bfracs(1:4),mean(r(imidx,7:10)),std(r(imidx,7:10))./sqrt(pcount),'k+--');
hold on
adjrat=r(imidx,12)./r(imidx,11);
r2valmax=r(imidx,7:10).*repmat(adjrat,[1 4]);
errorbar(bfracs(1:4),mean(r2valmax),std(r2valmax)./sqrt(pcount),'k+-');
errorbar(1,mean(r(imidx,11)),std(r(imidx,11))./sqrt(pcount),'k+-');
errorbar(1,mean(r(imidx,12)),std(r(imidx,12))./sqrt(pcount),'ko-');
hold off
axis([0 1.2 0 1]);

title('image domain model');


pcount=length(pfidx);
subplot(1,2,2);
errorbar(bfracs(1:4),mean(r(pfidx,7:10)),std(r(pfidx,7:10))./sqrt(pcount),'k+--');
hold on
adjrat=r(pfidx,12)./r(pfidx,11);
r2valmax=r(pfidx,7:10).*repmat(adjrat,[1 4]);
errorbar(bfracs(1:4),mean(r2valmax),std(r2valmax)./sqrt(pcount),'k+-');
errorbar(1,mean(r(pfidx,11)),std(r(pfidx,11))./sqrt(pcount),'k+-');
errorbar(1,mean(r(pfidx,12)),std(r(pfidx,12))./sqrt(pcount),'ko-');
hold off
axis([0 1.2 0 1]);

title('fourier power model');

tier1=mean([r(imidx,10) r(imidx,12)./r(imidx,11).*r(imidx,10) r(imidx,12) ...
            r(pfidx,10) r(pfidx,12)./r(pfidx,11).*r(pfidx,10) r(pfidx,12)])
s1=std([r(imidx,10) r(imidx,12)./r(imidx,11).*r(imidx,10) r(imidx,12) ...
        r(pfidx,10) r(pfidx,12)./r(pfidx,11).*r(pfidx,10) r(pfidx,12)])./...
   sqrt(length(imidx));
tier2=[mean([r(imidx0,10) r(imidx0,12)./r(imidx0,11).*r(imidx0,10) ...
             r(imidx0,12)]) ...
       mean([r(pfidx0,10) r(pfidx0,12)./r(pfidx0,11).*r(pfidx0,10) ...
             r(pfidx0,12)])]
s2=[std([r(imidx0,10) r(imidx0,12)./r(imidx0,11).*r(imidx0,10) ...
         r(imidx0,12)]) ./ sqrt(length(imidx0)) ...
    std([r(pfidx0,10) r(pfidx0,12)./r(pfidx0,11).*r(pfidx0,10) ...
         r(pfidx0,12)]) ./ sqrt(length(pfidx0))]

figure(2);
subplot(1,2,1);
bar([tier1]);
hold on
for ii=1:length(s1),
   plot([ii ii],tier1(ii)+[0 s1(ii)],'k-');
end
hold off

title(sprintf('%d cells good in both',length(imidx)));
axis([0 7 0 1]);

colormap(gray);


subplot(1,2,2);
bar([tier2]);
hold on
for ii=1:length(s2),
   plot([ii ii],tier2(ii)+[0 s2(ii)],'k-');
end
hold off

title(sprintf('im%d/pf%d cells good in one',length(imidx0),length(pfidx0)));
axis([0 7 0 1]);

set(gcf,'PaperOrientation','portrait','PaperPosition',[1.5 2.5 6 2.5]);
print -depsc /auto/k1/david/docs/2004/dissertation/eps/3.noisesum.eps

