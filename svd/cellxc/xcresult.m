% function xcresult(outfile,nlidxsave,savefigs);
function xcresult(outfile,nlidxsave,savefigs);

fprintf('xcresult.m: %s\n',outfile);

outfile=nslmakelocal(outfile);

if strcmp(outfile(end-2:end),'.gz'),
   zload(outfile);
else
   load(outfile);
end

if exist('nlidxsave','var') && ~isempty(nlidxsave),
   nlidx=nlidxsave;
elseif isfield(params,'nlidxsave'),
   nlidx=params.nlidxsave;
end
if ~exist('savefigs','var'),
    savefigs=0;
end
if nlidx>length(strf),
   nlidx=1
end

if params.batch==10 | params.batch==11,
   nlidx=3;
end
if ismember(params.batch,[51 52 55 56 59 60 76 80]),
   xc(:,:,3)=xc(:,:,nlidx);
   nlidx=3;
end
if iscell(xc),
   xc=xc{end};
end
if size(xc,3)<nlidx,
   xc(:,:,nlidx)=xc(:,:,1);
end
fprintf('chose to show nlidx=%d (%s)\n',nlidx,strf(nlidx).nltype);

hf=strf(nlidx).h;

if length(iconside)==1,
   iconside=[round(sqrt(iconside)) round(sqrt(iconside))];
end

stephf=cumsum(hf,2);
steptime=sum(stephf,1);
maxidx=min(find(steptime==max(steptime)));
kern=stephf(:,maxidx);

% examples of kernel at different cutoffs
sampcount=6;
sampidx=[round(linspace(1,params.sfscount,sampcount))];

if sigrange(1)==0,
   siguse=sigrange(2);
else
   siguse=sigrange(1);
end

if size(mH,2)>=size(hf,2)-params.maxlag(1),
   smm=mH(:, (1:size(hf,2))-params.maxlag(1), sampidx,1);
   sms=eH(:, (1:size(hf,2))-params.maxlag(1), sampidx,1) .* siguse;
else
   smm=mH(:, :, sampidx,1);
   sms=eH(:, :, sampidx,1) .* siguse;
end

smd=abs(smm)./(sms+(sms==0));
if ~params.shrinkage, % old--shrinkage filter
   tsf=smm.*(smd>1);
elseif sum(smd(:)>1)>0,
   smd=(1-smd.^(-2));
   smd=smd.*(smd>0);
   smd(find(isnan(smd)))=0;
   tsf=smm.*smd;
else
   tsf=smm;
end

if ~exist('lambda','var'),
   lambda=1:params.sfscount;
   keyboard
end
for ii=1:sampcount,
   tl=lambda(sampidx(ii));
   titles{ii}=sprintf('%s: STRF sample %.2f NLxc=%0.3f',...
                      params.cellid,tl,xc(sampidx(ii),1));
end

%disp('forcing kernfmt to be pfftgr!!!!!');
%params.kernfmt='pfftgr';

if 0 & ((strcmp(params.kernfmt,'sca') | strcmp(params.kernfmt,'pixel'))),
   loadSCAmatrix
   for ii=1:sampcount,
      showSCAKernels(tsf(:,:,ii), S, 20);
   end
elseif isfield('params','batch') & ismember(params.batch,[115 118:121]),
   f1=figure;
   % to emphasize small coefficients, use this:
   showkern((abs(tsf)).^(1/3).*sign(tsf),params.kernfmt,iconside,titles,1,16);
else
   f1=figure;
   showkern(tsf,params.kernfmt,iconside,titles,1,16);
end

f2=figure;
subplot(2,1,1);
imagesc(xc(:,:,nlidx,1,1));
hold on
plot(strf(nlidx,end,1).parms.sigfit,strf(nlidx,end,1).parms.sfsfit,'x');
hold off
colormap(hot);
colorbar;
title('fit xc');

cnfidx=1;
if isfield(params,'predbatch') && isfield(params,'batch') && ...
      ~isempty(params.batch),
   cnfidx=find(cat(1,params.predbatch{:})==params.batch);
end
if isempty(cnfidx),
   cnfidx=1;
end
if length(predres(cnfidx).predxc)>0 & ~isnan(predres(cnfidx).predxc(nlidx)),
   subplot(2,1,2);
   pp=predres(cnfidx).mod_psth{1}(:,1);
   ppnl=predres(cnfidx).mod_psth{1}(:,nlidx);
   rr=predres(cnfidx).act_resp{1};
   gidx=find(~isnan(pp) & ~isnan(rr));
   pp=pp(gidx);
   ppnl=ppnl(gidx);
   rr=rr(gidx);
   
   bincount=10;
   [pp ppidx]=sort(pp);
   ppnl=ppnl(ppidx);
   rr=rr(ppidx);
   edgeidx=round(linspace(1,length(pp)+1,bincount+1));
   
   ppbinned=zeros(bincount,1);
   rrbinned=zeros(bincount,1);
   for bidx=1:bincount
      ppbinned(bidx)=mean(pp(edgeidx(bidx):(edgeidx(bidx+1)-1)));
      ppnlbinned(bidx)=mean(ppnl(edgeidx(bidx):(edgeidx(bidx+1)-1)));
      rrbinned(bidx)=mean(rr(edgeidx(bidx):(edgeidx(bidx+1)-1)));
   end
   
   if length(pp)>500,
      [tt1, tt2]=sort(rand(size(pp)));
      pp=pp(tt2(1:500));
      rr=rr(tt2(1:500));
   end
   
   plot(ppbinned,rrbinned,'rx-');
   hold on
   plot(ppnlbinned,rrbinned,'k--');
   scatter(pp,rr,'.');
   hold off
   xlabel('predicted resp');
   ylabel('actual resp');
   legend('linear','outputnl');
else
   plot(steptime,'k-');
   hold on
   plot(maxidx,steptime(maxidx),'ro');
   hold off
   title('step response (o peak)');
end

if isnumeric(params.predbatch),
   params.predbatch={params.predbatch};
end

if length(params.predbatch)>1,
   bstr=[];
   for ii=1:length(params.predbatch),
      if ~isnan(predxc(ii,nlidx,1,1)),
         bstr=[bstr, sprintf('%d: %.2f/',params.predbatch{ii},predxc(ii,nlidx,1,1))];
      end
   end
   bstr=bstr(1:end-1);
else
   bstr=sprintf('%.2f',predxc(nlidx));
end

sfsfit=strf(nlidx,1,1).parms.sfsfit;
sigfit=strf(nlidx,1,1).parms.sigfit;
if ~isfield(params,'batch') | isempty(params.batch),
   params.batch=0;
end

titles={sprintf('%s/%d Impulse response (fitxc=%.2f predxc=%s)',...
                params.cellid,params.batch,...
                xc(sfsfit,sigfit,nlidx,1,1),bstr),...
        sprintf('%s Step response',params.cellid)};

if strcmp(params.kernfmt,'spect'),
   hlen=size(hf,2);
else
   hlen=min([12 size(hf,2)]);
end
if exist('resp','var'),
   respcount=size(resp,2);
elseif params.respfmtcode==0
   respcount=size(predxc,4);
end

if strcmp(params.kernfmt,'sca'), % | strcmp(params.kernfmt,'pixel'),
   loadSCAmatrix
   showSCAKernels(hf(:,1:hlen), S, 20);
   
elseif respcount>1,
   
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
   
   f3=figure;
   showkern(cat(3,hf(:,1:hlen,:)),params.kernfmt,iconside,titles,1,16);
   
else
   f3=figure;
   showkern(cat(3,hf(:,1:hlen),stephf(:,1:hlen)),params.kernfmt,...
            iconside,titles,1,16);
end

if savefigs,
    disp('printing results to jpg');
    
    if onseil==1
        localoutpath=['/homes/svd/data/web/results/batch',num2str(params.batch),'/'];
        finaloutpath=['/auto/data/web/results/batch',num2str(params.batch),'/'];
    else
        localoutpath=['/auto/data/web/results/batch',num2str(params.batch),'/'];
    end
    
    sfigure(f1);
    fullpage portrait
    sfigure(f3);
    fullpage portrait
    drawnow
    
    %['mkdir -p ',localoutpath]
    unix(['mkdir -p ',localoutpath]);
    print(['-f',num2str(f1)],'-djpeg','-r150',...
          sprintf('%s%s.1.jpg',localoutpath,params.cellid));
    print(['-f',num2str(f2)],'-djpeg','-r150',...
          sprintf('%s%s.2.jpg',localoutpath,params.cellid));
    print(['-f',num2str(f3)],'-djpeg','-r150',...
          sprintf('%s%s.3.jpg',localoutpath,params.cellid));
    
    if onseil==1
        %['ssh bhangra.isr.umd.edu mkdir -p ',finaloutpath]
        unix(['ssh bhangra.isr.umd.edu mkdir -p ',finaloutpath]);
        %['scp ',localoutpath,params.cellid,...
        %      '*jpg svd@bhangra.isr.umd.edu:',finaloutpath]
        unix(['scp ',localoutpath,params.cellid,...
              '*jpg svd@bhangra.isr.umd.edu:',finaloutpath]);
    end
end

% only execute following code for stephen
if ~strcmp(getenv('USER'),'david'),
   return
end

predbatchcount=length(params.predbatch);
if predbatchcount==1 & params.predbatch{1}==0,
   predbatchcount=0;
end

if predbatchcount>0,
   for predidx=1:predbatchcount,
      
      % figure out pred files for the current batch
      [pcellfiledata,ptimes,pbatchdata]=...
          cellfiletimes(params.cellid,params.predbatch{predidx},1);
      
      % does this cell have data for batchid=predidx?
      if length(pcellfiledata)>0,
         predparams=params;
         predparams.stimfiles={};
         predparams.respfiles={};
         predparams.stimcrfs=[];
         for ii=1:length(pcellfiledata),
            predparams.respfiles{ii}=[pcellfiledata(ii).path,...
                    pcellfiledata(ii).respfile];
         end
         
         tpredstartframe=ptimes(3).start;
         tpredstopframe=ptimes(3).stop;
         tpredfile=ptimes(3).fileidx;
         
         tresp=respload(predparams.respfiles{tpredfile(1)},...
                        'r',1,1,0);
         if size(tresp,2)>1,
            tresp=tresp(tpredstartframe:tpredstopframe,2:end);
            tresp=compact_raster_matrix3(tresp);
            
            firstwithnans=min([find(sum(isnan(tresp))>0) size(tresp,2)+1]);
            
            if firstwithnans>2,
               tresp=tresp(:,1:firstwithnans-1);
            end
         else
            tresp=tresp(tpredstartframe:tpredstopframe);
         end
         
         %keyboard
         
         repcount=size(tresp,2);
         fprintf('predidx=%d (bat %d) (%d reps): ceiling... ',...
                 predidx,params.predbatch{predidx},repcount);
         
         gg=find(~isnan(tresp(:,1)));
         [mu,alpha,beta]=reversepoisson(tresp(gg,1));
         rmax=singletrialceiling(tresp(gg,1),alpha,beta);
         %figure(gcf);
         drawnow
         
         xct=zeros(repcount,1);
         fprintf(' a=%.3f, b=%.3f\n',alpha,beta);
         for nn=1:size(strf,1),
            tpred=predres(predidx).mod_psth{1}(:,1,nn);
            
            if size(tresp,1)<length(tpred),
               tpred=tpred(1:end+params.maxlag(1));
            end
            if size(tresp,1)<length(tpred),
               tpred=tpred(params.maxlag(2)+1:end);
            end
            drawnow
            if size(tresp,1)>length(tpred)
               tresp=tresp(1:end+params.maxlag(1),:);
            end
            
            for xcidx=1:repcount,
               if length(tpred)~=length(tresp(:,xcidx)),
                  keyboard
               end
               gg=find(~isnan(tpred) & ~isnan(tresp(:,xcidx)));
               xct(xcidx)=xcov(tpred(gg),tresp(gg,xcidx),0,'coeff');
            end
            
            projr=sqrt(mean(xct.^2)./rmax.^2);
            fprintf(' nl=%d: %.3f/%.3f ---> %.3f (vs. %.3f)\n',...
                    nn,sqrt(mean(xct.^2)),rmax,projr,...
                    predres(predidx).predxc(nn));
            
         end
      else
         % no val data for this batch
      end
   end
end
