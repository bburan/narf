% function r=raw_nl2(beta,pred);
%
function r=raw_nl2(beta,pred);

pp=beta{1};
rr=beta{2};
if length(beta)>2,
    v=beta{3};
    pred=pred*v;
end

if isempty(pp),
    warning('empty beta, returning linear prediction');
    r=pred;
    return
end

bincount=size(pp);
tcount=size(pred,1);
r=zeros(tcount,1);

if std(pp(:))>0 && std(rr(:))>0,
   xx=pp;
   xx=[xx(1,:)-(xx(2,:)-xx(1,:))*5;
      xx;
      xx(end,:)+(xx(end,:)-xx(end-1,:))*5];
   yy=rr;
   yy=[yy(1,:)-(yy(2,:)-yy(1,:))*5;
      yy;
      yy(end,:)+(yy(end,:)-yy(end-1,:))*5];
   yy=[yy(:,1)-(yy(:,2)-yy(:,1))*5 ...
      yy ...
      yy(:,end)+(yy(:,end)-yy(:,end-1))*5];
   [x1,x2] = meshgrid(xx(:,1),xx(:,2));

   inrange=find(pred(:,1)>=xx(1,1) & pred(:,1)<xx(end,1) &...
      pred(:,2)>=xx(1,2) & pred(:,2)<xx(end,2));
      
   %I = INTERP1(X,Y,XI,'method','extrap')
   %if sum(diff(xx)==0)>0,
   %    xx=xx+eps.*(1:length(xx))';
   %end
   ff=find(diff(round(xx./max(xx(:)).*1000000))==0);
   if ~isempty(ff),
      keepidx=setdiff(1:length(xx),ff+1);
      xx=xx(keepidx);
      yy=yy(keepidx);
   end
   r=interp2(x1,x2,yy,pred(:,1),pred(:,2),'linear');
      
end

