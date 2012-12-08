function [u,s,v,cc]=songdbspace(tevalmtx,bootcount)

if ~exist('bootcount','var') || bootcount<1,
   bootcount=1;
end

u=[];
s=[];
v=[];
ecount=size(tevalmtx,1);
stepsize=ecount./bootcount;
ucount=size(tevalmtx,2);
testshuffle=shuffle(1:ecount);
for bb=1:bootcount,
   keeprange=ones(ecount,1);
   if bootcount>1,
      keeprange(testshuffle(round((bb-1).*stepsize+1):...
                            round(bb.*stepsize)))=nan;
   end
   
   % compute correlation matrix for big users
   cc=zeros(ucount);
   n=zeros(ucount);
   for u1=1:ucount,
      for u2=u1:ucount,
         gii=find(~isnan(tevalmtx(:,u1)) & ~isnan(tevalmtx(:,u2)) &...
                  ~isnan(keeprange));
         if length(gii)>0,
            n(u1,u2)=length(gii);
            n(u2,u1)=n(u1,u2);
            cc(u1,u2)=mean(tevalmtx(gii,u1).*tevalmtx(gii,u2));
            cc(u2,u1)=cc(u1,u2);
         end
      end
   end
   cc(isnan(cc))=0;
   
   % regularize
   lown=find(n<10);
   cc(lown)=cc(lown).*n(lown)./10;
   
   % project songs onto 1...svdcount eigenvectors. this info
   % will be saved in dEigProj
   
   [u(:,:,bb),s(:,:,bb),v(:,:,bb)]=svd(cc);
end

svdcount=20;

xc=zeros(svdcount,ucount);
ee=zeros(svdcount,ucount);
for svduse=1:svdcount,
   rr=1:svduse;
   
   predeval=zeros(ecount,ucount);
   
   for bb=1:bootcount,
      if bootcount>1,
         fitrange=testshuffle([1:round((bb-1).*stepsize) ...
                    round(bb.*stepsize+1):end]);
         valrange=testshuffle(round((bb-1).*stepsize+1):round(bb.*stepsize));
      end
      
      ei=v(:,rr,bb); %*(s(rr,rr,bb)^-1);
      
      for uu=1:ucount,
         a=tevalmtx(fitrange,:);
         a=a(find(~isnan(a(:,uu))),:);
         y=a(:,uu);
         a(isnan(a))=0;
         x=a(:,[1:(uu-1) (uu+1):end])*...
           u([1:(uu-1) (uu+1):end],rr,bb);
         
         b=regress(y,x);
         
         a=tevalmtx(valrange,:);
         a(find(isnan(a)))=0;
         evaleigs=a(:,[1:(uu-1) (uu+1):end])*...
           u([1:(uu-1) (uu+1):end],rr,bb);
         
         predeval(valrange,uu)=evaleigs(:,rr)*b;
      end
   end
   
   for uu=1:ucount,
      ii=find(~isnan(tevalmtx(:,uu)));
      xc(svduse,uu)=xcov(tevalmtx(ii,uu),predeval(ii,uu),0,'coeff');
      ee(svduse,uu)=sqrt(mean((tevalmtx(ii,uu)-predeval(ii,uu)).^2));
      % ee already normalized since nanstd(tevalmtx,uu)==1
   end
   
   [svduse sqrt(nanmean(nanmean((predeval-tevalmtx).^2)))]
end

ucc=sum(~isnan(tevalmtx));
[ucc,ui]=sort(-ucc);
ucc=-ucc;
subplot(3,1,1);
plot(ucc./10000);
hold on;
plot(xc([1 5 10],ui)');
hold off
subplot(3,1,2);
imagesc(xc(:,ui));
colorbar
subplot(3,1,3);
plot([mean(xc(:,ui(1:15)),2) mean(xc(:,ui(1:25)),2) mean(xc(:,ui(1:end)),2) ]);

keyboard