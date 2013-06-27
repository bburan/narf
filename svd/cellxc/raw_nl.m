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
    r=nanmean(pred,2);
    return
end

bincount=size(pp,1);
dimcount=size(pp,2);

tcount=size(pred,1);
r=zeros(tcount,1);

if dimcount==1,
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
else
    % assume dimcount==2,
    m1=4.5;m2=m1-1;
    pp=cat(1,pp(1,:)-(pp(2,:)-pp(1,:)).*5,...
           pp,pp(end,:)+(pp(end,:)-pp(end-1,:)).*5);
    rr=[rr(1,1)*m1-rr(2,2).*m2 rr(1,:)*m1-rr(2,:).*m2 rr(1,end)*m1-rr(2,end-1).*m2;
        rr(:,1)*m1-rr(:,2).*m2 rr rr(:,end)*m1-rr(:,end-1).*m2;
        rr(end,1)*m1-rr(end-1,2).*m2 rr(end,:)*m1-rr(end-1,:).*m2 rr(end,end)*m1-rr(end-1,end-1).*m2];
        
    if sum(abs(pp(:,1)),1)==0,
       pp(:,1)=linspace(-0.00001,0.00001,length(pp));
    end
    if sum(abs(pp(:,2)),1)==0,
       pp(:,2)=linspace(-0.00001,0.00001,length(pp));
    end
    r=interp2(pp(:,1),pp(:,2),rr,pred(:,1),pred(:,2),'linear');
    
    if 0,
        sfigure(1);
        clf
        surf(pp(2:end-1,1),pp(2:end-1,2),rr(2:end-1,2:end-1));
        surf(pp(:,1),pp(:,2),rr);
        drawnow
    end
end

    
