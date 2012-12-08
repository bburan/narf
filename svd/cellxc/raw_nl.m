% function r=raw_nl(beta,pred);
%
function r=raw_nl(beta,pred);

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

bincount=size(pp,1);
dimcount=size(pp,2);

tcount=size(pred,1);
r=zeros(tcount,1);

for dd=1:dimcount,
   
   xx=pp(:,dd);
   xx=[xx(1)-(xx(2)-xx(1))*5;
       xx;
       xx(end)+(xx(end)-xx(end-1))*10];
   yy=rr(:,dd);
   yy=[yy(1)-(yy(2)-yy(1))*5;
       yy;
       yy(end)+(yy(end)-yy(end-1))*10];
   
   if std(xx)>0 && std(yy)>0,
      
      inrange=find(pred(:,dd)>=xx(1) & pred(:,dd)<xx(end));
      
      %I = INTERP1(X,Y,XI,'method','extrap')
      if sum(diff(xx)==0)>0,
          xx=xx+eps.*(1:length(xx))';
      end
      tp=interp1(xx,yy,pred(inrange,dd),'linear');
      
      r(inrange)=r(inrange)+tp;
      
      orange=find(pred(:,dd)>xx(end));
      r(orange)=r(orange)+yy(end);
      urange=find(pred(:,dd)<xx(1));
      r(urange)=r(urange)+yy(1);
      
   end
end
