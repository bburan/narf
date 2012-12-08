% function rbinned=respvarbins(r,startbins,stopbins,meanrat[=1]);
%
% created SVD 4/23/02
%
% bin 2nd dimension of resp in a bunch of different bins defined by
% startbins and stopbins.  
% given:  rsize=size(r)=[resplen respcount attcount]
% rbinned(:,ii,:)=mean(resp(:,startbins(ii):stopbins(ii),:),2)
% of course the mean is appropriately normalized, given the
% possibility that some of the entries may be Nan.
%
function rbinned=respvarbins(r,startbins,stopbins,meanrat);

disp('respvarbins.m: rebinning response');

if ~exist('meanrat','var'),
   meanrat=1;
end
if ~exist('randcount','var'),
   randcount=0;
end

rsize=size(r);
resplen=rsize(1);
respcount=rsize(2);
if length(rsize)>2,
   attcount=rsize(3);
else
   attcount=1;
end

bincount=length(startbins);

if length(find(isnan(r)))>0,
   usenan=1;
else
   usenan=0;
   r(r==-1)=nan;
end

rbinned=zeros(resplen,bincount,attcount);
for ii=1:bincount
   if ~meanrat,
      % bin value = mean over all frames in that bin
      tr=r(:,startbins(ii):stopbins(ii),:);
      rnan=isnan(tr);
      rnz=sum(~rnan,2);
      tr(find(rnan))=0;
      tr=sum(tr,2)./(rnz+(rnz==0));
      tr(find(rnz==0))=nan;
      rbinned(:,ii,:)=tr;
   else
      % bin value = projection onto mean resp per fixation
      tr=r(:,startbins(ii):stopbins(ii),:);
      
      % take out pfth's with fewer than 3 valid bins
      ss=sum(~isnan(tr(:,:)),2);
      tr(ss<3,:)=nan;
      rnan=isnan(tr);
      
      % find mean response per fixation
      mbin=nanmean(tr(:,:,1));
      
      % tile over entire response matrix, normalizing according to
      % the number of valid bins in each flash
      mscale=repmat(mbin,[size(tr,1) 1 size(tr,3)]);
      mscale(find(rnan))=0;
      msum=sum(mscale,2);
      msum(find(msum==0))=nan;
      mscale=mscale./repmat(msum,[1 stopbins(ii)-startbins(ii)+1 1]);
      tr(rnan)=0;
      
      % project each response onto the normalized mean
      % response vector
      rbinned(:,ii,:)=sum(tr.*mscale,2);
   end
end

% take sqrt to distribute noise more normally (? JLG suggestion 12/02)
% has no effect
%rbinned(~isnan(rbinned))=sqrt(rbinned(~isnan(rbinned)));

if randcount>0,
   % add some random "attentional states" on to end as a control
   if usenan,
      rrand=ones([resplen,bincount,randcount])*nan;
   else
      rrand=-ones([resplen,bincount,randcount]);
   end
   if attcount>1,
      mcount=floor(sum(sum(~isnan(rbinned(:,1,2:end))))./(attcount-1));
      %attonidx=find(sum(~isnan(rbinned(:,1,2:end)),3));
   else
      mcount=sum(~isnan(rbinned(:,1)));
      %attonidx=find(~isnan(rbinned(:,1)));
   end
   
   % all attentional states are fair game
   attonidx=find(~isnan(rbinned(:,1)));
   
   rridx=0;   
   for ii=1:ceil(randcount/(attcount-1)),
      [rr,tgoodidx]=sort(rand(length(attonidx),1));
      tgoodidx=attonidx(tgoodidx);
      for jj=1:(attcount-1),
         rridx=rridx+1;
         if rridx<=randcount,
            ttgood=tgoodidx(((jj-1)*mcount+1):(jj*mcount));
            rrand(ttgood,:,rridx)=rbinned(ttgood,:,1);
         end
      end
   end
   rbinned=cat(3,rbinned,rrand);
end

