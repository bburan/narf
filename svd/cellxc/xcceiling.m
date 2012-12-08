% function [ccexp,ccp,ccr]=xcceiling(pp,r2);
%
function [ccexp,ccp,ccr]=xcceiling(pp,r2);

M=size(r2,2);

ccp2=xcov(pp,mean(r2,2),0,'coeff').^2;

ccp0=zeros(M,1);
ccr0=zeros(M,1);
for mm=1:M
   
   ra=mean(r2(:,1:round(M/2)),2);
   rb=mean(r2(:,round(M/2)+1:end),2);
   r2=shift(r2,[0 1]);
   
   if sum(abs(ra))>0 & sum(abs(rb))>0,
      ccr0(mm)=xcov(ra,rb,0,'coeff').^2;
   else
      ccr0(mm)=0;
   end
   
   if sum(abs(pp))>0 & sum(abs(r2(:,1)))>0,
      ccp0(mm)=xcov(pp,r2(:,1),0,'coeff').^2;
   else
      ccp0(mm)=0;
   end

end
ccp=mean(ccp0);
ccr2=mean(ccr0);

ccrat= (1+sqrt(1./ccr2)) ./ (2 - M + M * sqrt(1./ccr2));
ccr=   2                 ./ (2 - M + M * sqrt(1./ccr2));
ccexp=ccp./ccr;

%keyboard

