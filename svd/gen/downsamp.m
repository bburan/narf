function x1=downsamp(x0,r,dim)

if ~exist('dim','var'),
    dim=1;
end

outsize=ceil(size(x0,dim)/r);

if dim==1,
   x0(outsize*r,1,1)=0;
   x1=zeros(outsize,size(x0,2),size(x0,3));
   for aa=1:size(x0,2),
       for bb=1:size(x0,3),
           t=reshape(x0(:,aa,bb),r,outsize);
           x1(:,aa,bb)=nanmean(t)';
       end
   end
elseif dim==2
   x0(1,outsize*r,1)=0;
   x1=zeros(size(x0,1),outsize,size(x0,3));
   for aa=1:size(x0,1),
       for bb=1:size(x0,3),
           t=reshape(x0(aa,:,bb),r,outsize);
           x1(aa,:,bb)=nanmean(t);
       end
   end
elseif dim==3
   x0(1,1,outsize*r)=0;
   x1=zeros(size(x0,1),size(x0,2),outsize);
   for aa=1:size(x0,1),
       for bb=1:size(x0,2),
           t=reshape(x0(aa,bb,:),r,outsize);
           x1(aa,bb,:)=nanmean(t);
       end
   end
end
