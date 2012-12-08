function t0=find_raw_nl(pred,resp,verbose);

if ~exist('verbose','var');
   verbose=0;
end
if nansum(resp)==0,
    warning('response is zero vector, no raw nl');
    t0={[],[]};
    return
end

dimcount=size(pred,2);
bincount=25;

pp=zeros(bincount,dimcount);
rr=zeros(bincount,dimcount);
rre=zeros(bincount,dimcount);

if 0 && dimcount>2,
    pred0=pred;
    [u,s,v]=svd(pred(:,2:end));
    v=[1 zeros(1,dimcount-1);
       zeros(dimcount-1,1) v];
    pred=pred*v;
else
    v=[];
end

residual_resp=resp;
for dd=1:dimcount,
   [ss,si1]=sort(pred(:,dd));
   tb=bincount;
   b=[];
   
   while length(b)<bincount && tb<length(si1),
       tb=tb+1;
       edges1=round(linspace(1,length(si1)+1,tb));
       [b,ui,uj]=unique(ss(edges1(1:(end-1)))');
   end
   edgeend=edges1(end);
   edges1=edges1(ui);
   
   if length(b)<bincount,
      b=[b repmat(b(end),[1 bincount-length(b)])];
      edges1=[edges1 repmat(edges1(end),[1 bincount-length(edges1)])];
      
      %t0={[],[]};
      %return
   end
   edges1=[edges1 edgeend];
   
   for bb=1:bincount,
      pp(bb,dd)=mean(pred(si1(edges1(bb):(edges1(bb+1)-1)),dd));
      rr(bb,dd)=mean(residual_resp(si1(edges1(bb):(edges1(bb+1)-1))));
      nn=sqrt(edges1(bb+1)-edges1(bb));
      if edges1(bb+1)>edges1(bb),
         rre(bb,dd)=std(residual_resp(si1(edges1(bb):(edges1(bb+1)-1))))./...
             (nn+(nn==1));
      end
   end
   pp(isnan(pp))=0;
   rr(isnan(rr))=0;
   
   rr(:,dd)=gsmooth(rr(:,dd),1);
   %[pp(:,dd),rr(:,dd)]
   intpred=raw_nl({pp(:,dd),rr(:,dd)},pred(:,dd));
   residual_resp=residual_resp-intpred;
   if dd==1,
       res_resp1=residual_resp;
   end
end




if verbose,
    % joint analysis
    bincount=10;
    d1=2;
    d2=3;
    rr2=zeros(bincount,bincount);
    [ss,si1]=sort(pred(:,d1));
    edges1=round(linspace(1,length(si1)+1,bincount+1));
    [ss,si2]=sort(pred(:,d2));
    edges2=round(linspace(1,length(si2)+1,bincount+1));
    for bb1=1:bincount,
        for bb2=1:bincount,
            tt=intersect(si1(edges1(bb1):(edges1(bb1+1)-1)),si2(edges2(bb2):(edges2(bb2+1)-1)));
            if length(tt)>0,
                rr2(bb1,bb2)=mean(res_resp1((tt)));
            end
        end
    end

    clf
    subplot(2,1,1);
    errorbar(pp,rr,rre);
    
    subplot(2,1,2);
    imagesc(rr2);
end

if isempty(v),
    t0={pp,rr};
else
    t0={pp,rr,v};
end

