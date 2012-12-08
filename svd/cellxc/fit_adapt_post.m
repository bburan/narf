function beta=fit_adapt_post(pred,resp);

p1=pred;
p2=conv2(pred(:),[0 0 0 0 1/3 1/3 1/3]','same');
ggidx=find(~isnan(resp(:)+p1(:)+p2(:)));
resp=resp(ggidx);
pred=pred(ggidx);

pv=std(pred);
if pv==0,
   pv=1;
end


beta0=[0 1./pv 1];
lb=[max(0,min(pred))    0.1./pv  0.5];
ub=[max(pred).*0.9+0.1  3./pv     5  ];

fitopt=optimset('Display','off');
beta=lsqcurvefit('adapt_post',beta0,pred,resp,lb,ub,fitopt);
