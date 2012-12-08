function pred=kvapred(vstrf,stim,bresp,params);



% figure out number of valid fixations in each attentional condition
ncount=squeeze(sum(~isnan(bresp(:,end,:)),1));
cumncount=[0;cumsum(ncount(2:end))];
anyokidx=find(sum(~isnan(bresp(:,end,2:end)),3));

% choose a random set of fixations for each bootstrap. this should
% avoid bias in the non-randomized noiseidx=1
vcount=length(anyokidx);

ss={};
for attidx=2:attcount,
   aokidx=find(~isnan(bresp(:,end,attidx)));
   [vv,ss{attidx-1}]=sort(rand(size(aokidx)));
   ss{attidx-1}=aokidx(ss{attidx-1});
end

vidx=[];
for bootidx=1:params.bootcount,
   for attidx=2:attcount,
      acount=length(ss{attidx-1});
      ppidx=(round((bootidx-1)/params.bootcount*acount)+1):...
            round(bootidx/params.bootcount*acount);
      
      vidx=[vidx; ss{attidx-1}(ppidx)];
   end
end

% generate random attention state sets
sn=zeros(length(anyokidx),params.noisecount+1);
for noiseidx=1:params.noisecount+1;
   [tt,sn(:,noiseidx)]=sort(rand(size(anyokidx)));
end


%%
%% DO THE PREDICTIONS
%%

realattidx=zeros(size(bresp,1),1);
for attidx=2:attcount,
   tokidx=find(~isnan(bresp(:,1,attidx)));
   realattidx(tokidx)=attidx-1;
end

for noiseidx=1:params.noisecount+1,
   fprintf('noiseidx=%3d: ',noiseidx);
   
   resp=ones(size(bresp))*nan;
   resp(anyokidx,:,1)=bresp(anyokidx,:,1);
   if noiseidx==1 | ~params.randcnf,
      resp(:,:,2:end)=bresp(:,:,2:end);
   else
      % pick random atts for each fixation. effectively the same as 0?
      for attidx=2:attcount,
         nn=anyokidx(sn(cumncount(attidx-1)+1:cumncount(attidx),noiseidx));
         resp(nn,:,attidx)=bresp(nn,:,1);
      end
   end
   
   % assume there's only ONE response dimension!
   rbinidx=1;
   
   counter=zeros(blen,nlcount);
   for nlidx=1:nlcount,
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % JACKKNIFE PREDICTIONS
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      % set up vectors to hold preds from different models
      rpred=zeros(blen,1);
      rpred0=zeros(blen,1);
      
      % generate pred for each bootstrapped segment
      for bootidx=1:params.bootcount;
         a0idx=vidx((round((bootidx-1)/params.bootcount*vcount)+1):...
                    round(bootidx/params.bootcount*vcount));
         
         % pull out stimulus
         for attidx=2:attcount
            
            % figure out range to predict
            aidx=a0idx(find(~isnan(resp(a0idx,rbinidx,attidx))));
            counter(aidx,nlidx)=counter(aidx,nlidx)+1;
            
            % get kernel
            tstrf=vstrf(nlidx,attidx-1,rbinidx,bootidx);
            tH=tstrf.h;
            tmS=tstrf.mS';
            tnltype=tstrf.nltype;
            tnlparms=tstrf.nlparms;
            attname=tstrf.parms.attname;
            
            % do the linear prediction
            estim=bstim(aidx,:)-repmat(tmS,length(aidx),1);
            tpred=estim * tH;
            
            if ~isempty(tnltype) & ~strcmp(tnltype,'none'),
               tpred=feval(tnltype,tnlparms,tpred);
            end
            
            
            % measure vector angle between each stimulus and the kernel
            estim0=sqrt(sum(estim.^2,2));
            estim0=estim./repmat(estim0,[1 spacecount]);
            if norm(tH)>0,
               tH0=tH./norm(tH);
            else
               tH0=zeros(size(tH));
            end
            rpred0(aidx)=acos(estim0 * tH0)./pi*180;
            
            if 0,
               figure(3);
               plot(rpred(anyokidx),'g')
               figure(4);
               plot(tstrf.h);
               [bootidx attidx tstrf.nlparms']
               pause
            end
            
            rpred(aidx)=tpred;
         end
         
         % rectify... this seems like an ok thing to do. since
         % actual responses are always >=0 and this is unbiased
         %rpred(find(rpred<0))=0;
         
         
         ract=bresp(a0idx,rbinidx,1);
         if std(rpred(a0idx))>0 & std(ract)>0,
            xcboot(noiseidx,nlidx,bootidx)=xcov(rpred(a0idx),ract,0,'coeff');
         end
      end
      
      % evaluate against actual resp
      for attidx=1:attcount,
         tokidx=find(~isnan(bresp(:,1,attidx)));
         ract=bresp(tokidx,rbinidx,1);
         
         % within each attention state, cc will be biased toward
         % higher values because dc/gain changes won't matter
         if length(ract)>0 & std(ract)>0,
            xccross(noiseidx,nlidx,attidx)=xcov(rpred(tokidx),ract,0,'coeff');
         else
            xccross(noiseidx,nlidx,attidx)=0;
         end
      end
      
      for pidx=1:paircount,
         tokidx=find(~isnan(bresp(:,1,pairidx(pidx,1)+1)) | ...
                     ~isnan(bresp(:,1,pairidx(pidx,2)+1)));
         ract=bresp(tokidx,rbinidx,1);
         if length(ract)>0 & std(ract)>0,
            xcpair(noiseidx,nlidx,pidx)=xcov(rpred(tokidx),ract,0,'coeff');
         end
      end
      
      ract=bresp(anyokidx,rbinidx,1);
      if noiseidx==1 & nlidx==1,
         % get expected xc for randomly shuffled pred... to
         % determine whether this is a reasonable pred or not
         [xc(noiseidx,nlidx),randxc,tt,p]=...
             randxcov(rpred(anyokidx),ract,0,100);
      elseif length(ract)>0 & std(ract)>0 & std(rpred(anyokidx))>0,
         xc(noiseidx,nlidx)=xcov(rpred(anyokidx),ract,0,'coeff');
      end
      
      % not doing full xc.  average over booted xc
      xc(noiseidx,attidx)=mean(xcboot(noiseidx,attidx,:));
      
      fprintf(' %.3f',xc(noiseidx,nlidx));
      
      if noiseidx==1,
         rprec(:,nlidx)=rpred;
         rangle(:,nlidx)=rpred0;
      end
   end
   
   fprintf('\n');
   
   % update queue if active
   dbsetqueue;
end
