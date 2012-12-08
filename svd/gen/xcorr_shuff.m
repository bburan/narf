% function [xc,xc_err,xc_shuff,xc_shuff_err]=xcorr_shuff(d1,d2,maxbin);
%
% xc(tau) is cross correlation between d1 and d2: d1(t)*d2(t-tau)
% ie, for positive tau, d1 follows d2
%     for negative tau, d2 follows d1
%
% created SVD 2007-08
%
function [xc,xc_err,xc_shuff,xc_shuff_err]=xcorr_shuff(d1,d2,maxbin);

repcount=size(d1,2);

% not used currently, subtract mean of all trials
md1=mean(d1(:));
md2=mean(d2(:));

shuffcount=5;
txc=zeros(maxbin*2+1,repcount);
sxc=zeros(maxbin*2+1,shuffcount.*repcount);

for ii=1:repcount,
   % subtract mean for each trial (xcov)
   md1=mean(d1(:,ii));
   md2=mean(d2(:,ii));
   txc(:,ii)=xcorr(d1(:,ii)-md1,d2(:,ii)-md2,maxbin,'unbiased');
end

for ss=1:shuffcount,
   for ii1=1:repcount,
      ii2=mod(ii1+round(repcount.*(ss./(shuffcount+1))),repcount)+1;
      %ii2=ceil(rand.*(repcount-1));
      %if ii2>=ii1,
      %   ii2=ii2+1;
      %end
      % subtract mean for each trial (xcov)
      md1=mean(d1(:,ii1));
      md2=mean(d2(:,ii2));
      ii=(ss-1).*repcount+ii1;
      sxc(:,ii)=xcorr(d1(:,ii1)-md1,d2(:,ii2)-md2,maxbin,'unbiased');
   end
end

xc=mean(txc,2);
xc_err=std(txc,0,2)./sqrt(repcount);
xc_shuff=mean(sxc,2);
xc_shuff_err=std(sxc,0,2)./sqrt(repcount);

return
sfigure(2);
clf
errorshade((-maxbin:maxbin)',xc_shuff,xc_shuff_err,[1 0 0]);
hold on
errorshade((-maxbin:maxbin)'+maxbin.*2+1,xc,xc_err,[0 0 1]);
hold off
axis tight
