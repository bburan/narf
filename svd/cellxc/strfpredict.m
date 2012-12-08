% function [r]=strfpredict(strf,stim)
%
% strf      - strf structure
% stim      - time X spacecount :  pre-loaded movie matrix
%
function [p,p_lin]=strfpredict(strf,stim)

H=strf.h;
T1=-strf.zerobin+1;
T2=size(strf.h,2)+T1-1;

stimlen=size(stim,1);
spacecount=size(stim,2);

switch(strf.architecture),
 case 'linear',
  p_lin=zeros(stimlen,1);
  for t1=T1:T2,
     gidx=max(t1+1,1):min(t1+stimlen,stimlen);
     p_lin(gidx)=p_lin(gidx)+stim(gidx-t1,:)*H(1:spacecount,t1-T1+1);
  end
  
 case 'linear_sep',
  % same as linear, just don't sum over space
  p_lin=zeros(stimlen,spacecount);
  for t1=T1:T2,
     gidx=max(t1+1,1):min(t1+stimlen,stimlen);
     p_lin(gidx,:)=p_lin(gidx,:)+stim(gidx-t1,:)*H(1:spacecount,t1-T1+1,end);
  end
  
 case 'order2',
   
  % as of now, order 2 only works if spacecount==1
  if spacecount>1,
     error('sorry, order2 model not coded for spacecount>1');
  end
  
  p_lin=zeros(stimlen,1);
  
  shiftstim=zeros(stimlen,T2-T1+1);
  for os=0:(T2-T1),
     % t2 is always >= t1
     shiftstim((1+os):end,os+1)=stim(1:(end-os),:).*...
         stim((1+os):end,:);
  end
  
  for t1=T1:T2,
     gidx=max(t1+1,1):min(t1+stimlen,stimlen);
     p_lin(gidx)=p_lin(gidx)+stim(gidx-t1,:)*H(1,t1-T1+1,end);
     
     for t2=t1:T2,
        os=t2-t1;
        if os>0
           p_lin(gidx)=p_lin(gidx)+shiftstim(gidx-t1,os+1)* 2 *...
               H(t2-T1+2,t1-T1+1,end);
        else
           p_lin(gidx)=p_lin(gidx)+shiftstim(gidx-t1,os+1)*...
               H(t2-T1+2,t1-T1+1,end);
        end
     end
  end
  
 case 'stc',
   
  % as of now, order 2 only works if spacecount==1
  if spacecount>1,
     error('sorry, stc model not coded for spacecount>1');
  end
  if isfield(strf.parms,'h2'),
     H2=strf.parms.h2;
     maxs=size(H2,1)-1;
  else
     
     m2=H(2:end,:);
     [u,s,v]=svd(m2);
     ds=diag(s);
     maxs=max([0;find(ds>sum(ds)./10)]);
     
     % diagnostics
     %ds((maxs+1):end)=0;
     %imagesc(u(:,1)*ds(1)*v(:,1)' + u(:,2)*ds(2)*v(:,2)' +...
     %        u(:,3)*ds(3)*v(:,3)'+ u(:,4)*ds(4)*v(:,4)'-m2);
     
     H2=[H(1,:)./norm(H(1,:));u(:,1:maxs)'];
  end
  
  % same as linear, just don't sum over space
  p_lin=zeros(stimlen,1+maxs);
  for ss=1:(1+maxs),
     for t1=T1:T2,
        gidx=max(t1+1,1):min(t1+stimlen,stimlen);
        p_lin(gidx,ss)=p_lin(gidx,ss)+stim(gidx-t1,:)*H2(ss,t1-T1+1,end);
     end
  end
  
 case 'ML',

  r0=stim;
  cellcount=size(r0,2);
  resplen=size(r0,1);
  
  % smooth response to generate appropriate integration time
  smoothwin=strf.parms.smoothwin;
  fprintf('Response window=%d bins\n',smoothwin);
  smfilt=[zeros(smoothwin-1,1);ones(smoothwin,1)];
  r=conv2(r0,smfilt,'same');
  
  if isfield(strf.parms,'adddelay') && strf.parms.adddelay,
     % two responses, one offset by a smoothbin
     r=permute(r,[1 3 2]);
     rshift=shift(r,strf.parms.smoothwin);
     rshift(1:strf.parms.smoothwin,:,:)=nan;
     r=cat(2,r,rshift);
     r=reshape(r,size(r,1),size(r,2).*size(r,3));
     cellcount=size(r,2);
  end
  
  % bin the response
  disp('binning response');
  rbins=strf.parms.rbins;
  maxr=size(rbins,1)-1;
  ff=find(sum(isnan(r),2)==0);
  rbinned=zeros(size(r));
  for cc=1:cellcount,
     rbins(end,cc)=inf;
     for ii=1:maxr,
        rbinned(ff(find(r(ff,cc)>=rbins(ii,cc) & ...
                        r(ff,cc)<rbins(ii+1,cc))),cc)=ii;
     end
  end
  
  disp('calculating ML stim');
  plagcount=strf.parms.maxlag(2)-strf.parms.maxlag(1)+1;
  s=zeros(plagcount,length(r));
  if cellcount==2,
      for r1=1:maxr,
          for r2=1:maxr,
              rii=find(rbinned(:,1)==r1 & rbinned(:,2)==r2);
              s(:,rii)=repmat(strf.h(:,r1,r2),[1 length(rii)]);
          end
      end
  elseif cellcount==3,
      fprintf('cellcount==3\n');
      for r1=1:maxr,
          for r2=1:maxr,
              for r3=1:maxr,
                  rii=find(rbinned(:,1)==r1 & rbinned(:,2)==r2 &...
                           rbinned(:,3)==r3);
                  s(:,rii)=repmat(strf.h(:,r1,r2,r3),[1 length(rii)]);
              end
          end
      end
  end
  p_lin=s';
  
 case 'lnR',
  r0=stim;
  cellcount=size(r0,2);
  resplen=size(r0,1);
  
  % smooth response to generate appropriate integration time
  smoothwin=strf.parms.smoothwin;
  fprintf('Response window=%d bins\n',smoothwin);
  smfilt=[zeros(smoothwin-1,1);ones(smoothwin,1)];
  r=conv2(r0,smfilt,'same');
  
  if isfield(strf.parms,'adddelay') && strf.parms.adddelay,
     % two responses, one offset by a smoothbin
     r=permute(r,[1 3 2]);
     rshift=shift(r,strf.parms.smoothwin);
     rshift(1:strf.parms.smoothwin,:,:)=nan;
     r=cat(2,r,rshift);
     r=reshape(r,size(r,1),size(r,2).*size(r,3));
     cellcount=size(r,2);
  end
  
  p_lin=r*strf.h;
end

if ~isfield(strf,'nltype') || isempty(strf.nltype),
   p=p_lin;
   return;
end

switch strf.nltype

 case 'none',
  
  % linear, do nothing
  p=p_lin;
   
 case 'post-sum thresh',
   
  p=thresh(strf.nlparms,p_lin);
  
 case 'exp thresh',
   
  p=expthresh(strf.nlparms,p_lin);
   
 case 'pre-sum thresh',
  
  p=hinge2(strf.nlparms,p_lin);
 otherwise
  
  % assume nltype is a string
  p=feval(strf.nltype,strf.nlparms,p_lin);
end
