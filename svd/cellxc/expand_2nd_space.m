function x2=expand_2nd_space(x);

n=size(x,1);

if n<2,
   disp('fewer than 2 channels, no 2nd order expansion possible');
   x2=x;
   return
end

pairs=nchoosek(1:n,2);
paircount=size(pairs,1);
fprintf('second order expansion: %d to %d channels\n',...
        n,paircount+n);

x2=zeros(n+paircount,size(x,2));
x2(1:n,:)=x;
mm=nanmean(x,2);
for pp=1:paircount,
   % subtract mean?  problematic because mean changes
   %x2(pp+n,:)=(x(pairs(pp,1),:)-mm(pairs(pp,1))).*...
   %    (x(pairs(pp,2),:)-mm(pairs(pp,2)));
   % for now, don't subtract mean. this shouldn't really hurt,
   % should it? (since there's meaning to both cells firing at
   % once?)
   % actually, now mean subtraction is taken care of in loadsiteraster
   x2(pp+n,:)=x(pairs(pp,1),:).*x(pairs(pp,2),:);
end

