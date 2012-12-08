function strfout=strfdimremap(strf,stim,resp);

H=strf.h;
T1=-strf.zerobin+1;
T2=size(strf.h,2)+T1-1;

stimlen=size(stim,1);
spacecount=size(stim,2);
  
if ~strcmp(strf.architecture,'order2'),
   error('dimremap only works for order2 strfs');
end

if spacecount>1,
   error('sorry, dim remap not coded for spacecount>1');
end

m2=H(2:end,:);
[u,s,v]=svd(m2);
ds=diag(s);
maxs=max([0;find(ds>sum(ds)./10)]);
if maxs>2,
   maxs=2;
   disp('hard limit on maxs at 2');
end

p_lin=zeros(stimlen,1+maxs);
H2=[H(1,:)./(norm(H(1,:))+(norm(H(1,:))==0));u(:,1:maxs)'];

goodidx=find(~isnan(resp));
bincount=20;

REMAP=0;
if maxs>0 && REMAP,
   
   or=linspace(0,pi,21)';
   or=or(1:(end-1));
   xc=zeros(length(or),maxs);
   
   newH=H2;
   for dd=1:(maxs+1),
      if dd==1,
         d1=1;
         d2=2;
      else
         d1=1;
         d2=dd;
      end
      
      % same as linear, just don't sum over space
      p_lin=zeros(stimlen,1+maxs);
      for ss=1:(1+maxs),
         for t1=T1:T2,
            gidx=max(t1+1,1):min(t1+stimlen,stimlen);
            p_lin(gidx,ss)=p_lin(gidx,ss)+stim(gidx-t1,:)*newH(ss,t1-T1+1,end);
         end
      end
      if dd>1,
         thisfit=find_raw_nl(p_lin(goodidx,1:(dd-1)),resp(goodidx));
         thispred=raw_nl(thisfit,p_lin(:,1:(dd-1)));
         respres=resp-thispred;
         xcov(thispred(goodidx),resp(goodidx),0,'coeff')
      else
         respres=resp;
      end
      
      thisfit={};
      for ii=1:length(or),
         tp=cos(or(ii)).*p_lin(:,d1) + sin(or(ii)).*p_lin(:,d2);
         
         thisfit{ii}=find_raw_nl(tp(goodidx),respres(goodidx));
         thispred=raw_nl(thisfit{ii},tp);
         
         xc(ii,dd)=xcov(thispred(goodidx),respres(goodidx),0,'coeff');
      end
      
      mm=min(find(xc(:,dd)==max(xc(:,dd))));
      fprintf('best %d(%d,%d) fit for %.1f deg\n',dd,d1,d2,or(mm)/pi*180);
      tpmax=cos(or(mm)).*p_lin(:,d1) + sin(or(mm)).*p_lin(:,d2);
      
      newH(dd,:)=cos(or(mm)).*newH(d1,:) + sin(or(mm)).*newH(d2,:);
      newH(dd,:)=newH(dd,:)./std(newH(dd,:));
      %plot(newH');
      %drawnow
   end

   p_lin=zeros(stimlen,1+maxs);
   for ss=1:(1+maxs),
      for t1=T1:T2,
         gidx=max(t1+1,1):min(t1+stimlen,stimlen);
         p_lin(gidx,ss)=p_lin(gidx,ss)+stim(gidx-t1,:)*newH(ss,t1-T1+1,end);
      end
   end
   thisfit=find_raw_nl(p_lin(goodidx,:),resp(goodidx));
   thispred=raw_nl(thisfit,p_lin);
   respres=resp-thispred;
   xcov(thispred(goodidx),resp(goodidx),0,'coeff')
   
   % flip sign for display purposes
   for dd=1:size(p_lin,2),
      xx=xcov(p_lin(goodidx,dd),resp(goodidx),0,'coeff');
      if xx<0,
         newH(dd,:)=-newH(dd,:);
         thisfit{1}(:,dd)=flipud(-thisfit{1}(:,dd));
         thisfit{2}(:,dd)=flipud(thisfit{2}(:,dd));
      end
   end
else
   % force orthogonal
   newH=H2;
   
   % option: orthogonalize dimensions of strf
   if 0
   for dd=2:(1+maxs),
      t=newH(1:(dd-1),:)*newH(dd,:)';
      t=repmat(t,[1,size(newH,2)]).*newH(1:(dd-1),:);
      t=newH(dd,:)-sum(t,1);
      newH(dd,:)=t./norm(t)
   end
   end
   
   p_lin=zeros(stimlen,1+maxs);
   for ss=1:(1+maxs),
      for t1=T1:T2,
         gidx=max(t1+1,1):min(t1+stimlen,stimlen);
         p_lin(gidx,ss)=p_lin(gidx,ss)+stim(gidx-t1,:)*newH(ss,t1-T1+1,end);
      end
   end
   thisfit=find_raw_nl(p_lin(goodidx,:),resp(goodidx));
end

strfout=strf;
strfout.architecture='stc';
strfout.nlparms=thisfit;
strfout.nltype='raw_nl';
strfout.parms.h2=newH;

return


%resp=resp(gidx);
%stim=stim(gidx,:);
%p_lin=p_lin(gidx,:);

res_resp1=resp;

bincount=20;
linpred=p_lin;
for d1=2:(maxs-1),
   if d1>1,
      thisfit=find_raw_nl(linpred(goodidx,1:(d1-1)),resp(goodidx));
      thispred=raw_nl(thisfit,linpred(:,1:(d1-1)));
      res_resp1=resp-thispred;
   end
   
   d2=d1+1;
   
   rr2=zeros(bincount,bincount);
   [ss,si1]=sort(linpred(goodidx,d1));
   edges1=round(linspace(1,length(si1)+1,bincount+1));
   [ss,si2]=sort(linpred(goodidx,d2));
   edges2=round(linspace(1,length(si2)+1,bincount+1));
   for bb1=1:bincount,
      for bb2=1:bincount,
         tt=intersect(si1(edges1(bb1):(edges1(bb1+1)-1)),si2(edges2(bb2):(edges2(bb2+1)-1)));
         if length(tt)>0,
            rr2(bb1,bb2)=mean(res_resp1(goodidx(tt)));
         else
            rr2(bb1,bb2)=nan;
         end
      end
   end
   
   ff=find(isnan(rr2));
   gg=find(~isnan(rr2));
   [yy,xx]=meshgrid(1:bincount,1:bincount);
   xi=xx(gg);yi=yy(gg);
   rr3=rr2;
   for ii=1:length(ff),
      [jj,kk]=ind2sub(size(rr2),ff(ii));
      dd=abs((jj-xi).^2+(kk-yi).^2);
      nn=find(dd==min(dd));
      pp=sub2ind(size(rr2),xi(nn),yi(nn));
      rr3(jj,kk)=mean(rr2(pp));
   end
   rr3=gsmooth(rr3,[1 1]);
   
   or=(0:2.5:178.5)';
   testvar=zeros(size(or));
   for ii=1:length(or),
      tim=imrotate(rr3,or(ii),'bilinear');
      tim(tim==0)=nan;
      testvar(ii)=nansum(nanstd(tim));
   end
   testvar=cconv2(testvar,ones(5,1));
   plot(or,testvar);
   
   minordiff=25;
   mm1=min(find(testvar==max(testvar)));
   mm1dist=or-or(mm1);
   mm1dist(mm1dist>90)=mm1dist(mm1dist>90)-180;
   testvar2=testvar;
   testvar2(abs(mm1dist)<minordiff)=0;
   mm2=min(find(testvar2==max(testvar2)));
   fprintf('or(%d)=%.1f or(%d)=%.1f\n',d1,or(mm1),d2,or(mm2));
   
   newH=H2;
   newH(d1,:)=cos(-or(mm1)./180*pi)*H2(d1,:) + sin(-or(mm1)./180*pi)*H2(d2,:);
   %newH(d2,:)=cos(-or(mm2)./180*pi)*H2(d1,:) + sin(-or(mm2)./180*pi)*H2(d2,:);
   newH(d2,:)=cos((-or(mm2)+90)./180*pi)*H2(d1,:) + sin((-or(mm2)+90)./180*pi)*H2(d2,:);
   H2=newH;
   
   linpred=zeros(stimlen,1+maxs);
   for ss=1:(1+maxs),
      for t1=T1:T2,
         gidx=max(t1+1,1):min(t1+stimlen,stimlen);
         linpred(gidx,ss)=linpred(gidx,ss)+stim(gidx-t1,:)*H2(ss,t1-T1+1,end);
      end
   end
   
   keyboard
   
end

