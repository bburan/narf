% function [cc,n,mm,ss]=songdbcorrmtx(evalmtx);
%
% corrmtx of sparsely fillted evalmtx. each column normalized to
% have mean 0 and variance 1 before corr calculation (returned in
% mm and ss).  n tells the number of values that went into each cc
% entry.
%
function [cc,n,mm,ss]=songdbcorrmtx(evalmtx);

mm=nanmean(evalmtx);
ss=nanstd(evalmtx);
songcount=size(evalmtx,1);
usercount=size(evalmtx,2);

evalmtx=(evalmtx-repmat(mm,[songcount,1])) ./ ...
        repmat(ss,[songcount,1]);

cc=zeros(usercount);
n=zeros(usercount);
for u1=1:usercount,
   for u2=u1:usercount,
      n(u1,u2)=sum(~isnan(evalmtx(:,u1).*evalmtx(:,u2)));
      n(u2,u2)=n(u1,u1);
      cc(u1,u2)=nanmean(evalmtx(:,u1).*evalmtx(:,u2));
      cc(u2,u1)=cc(u1,u2);
   end
end

cc(isnan(cc))=0;


