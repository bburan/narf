% x is time X space
% y is time X 1

function h=blasso(x,y,minlag,maxlag,stepsize,tolerance);

if ~exist('minlag','var'),
   minlag=0;
end
if ~exist('maxlag','var'),
   maxlag=10;
end

if ~exist('stepsize','var'),
   stepsize=nanstd(y(:))./std(x(:)) ./ 100;
end
if ~exist('tolerance','var'),
   tolerance=0.001;
end

timecount=size(x,1);
spacecount=size(x,2);
tbincount=maxlag-minlag+1;

mS=nanmean(x);
sS=nanstd(x);


% compute cross correlation between stimulus at various time
% lags, stim(:-tt,:), and observed response, resp(:)
ttRR=mean(y.^2);
ttSS=repmat((sS(:)+(sS(:)==0)).^2,[1 tbincount]);
h=zeros(spacecount,tbincount);
tresp=y;
rgoodidx=find(~isnan(y));
rpred=ones(size(y)).*nan;
cumpred=zeros(size(y));

pstep=x.*stepsize; % this is the incremental prediction for each step

for stepidx=1:500,
   ttSR=zeros(spacecount,tbincount);
   tr=tresp;
   
   if stepidx>1,
      ttEE=ones(spacecount,tbincount).*inf;
      for tt=minlag:maxlag,
         trg=rgoodidx(find(rgoodidx-tt>0 & rgoodidx-tt<=timecount));
         tt0=tt-minlag+1;
         for xidx=find(h(:,tt0)~=0)',
            pp=cumpred(trg)-sign(h(xidx,tt0)).*pstep(trg-tt,xidx);
            ttEE(xidx,tt0)=sqrt(sum((y(trg)-pp).^2));
         end
      end
      
      maxh=min(find(abs(ttEE(:))==min(abs(ttEE(:)))));
      [xx,tt0]=ind2sub([spacecount tbincount],maxh);
      tt=tt0+minlag-1;
      
      hnew=h;
      coeff=-sign(hnew(xx,tt0)).*stepsize;
      hnew(xx,tt0)=hnew(xx,tt0)+coeff;
      
      gammanew=ttEE(xx,tt0) + lambda*sum(abs(hnew));
      
      trg=rgoodidx(find(rgoodidx-tt>0 & rgoodidx-tt<=timecount));
      lasterr=sqrt(sum((y(trg)-cumpred(trg)).^2));
      gammaold=lasterr + lambda*sum(abs(h));
      
      [gammanew gammaold gammanew-gammaold lambda sum(abs(h))./ ...
       stepsize]
      %if length(find(h(:,tt0)~=0)')>1
       %  keyboard
     % end
   else
      gammanew=1; gammaold=0;
   end
   
   if gammanew-gammaold<-tolerance,
      forward=0;
      %keyboard
   else
      forward=1;
      
      ttEE=zeros(spacecount,tbincount,2);
      for tt=minlag:maxlag,
         trg=rgoodidx;
         trg=trg(find(trg-tt>0 & trg-tt<=timecount));
         
         for xidx=1:spacecount,
            ttEE(xidx,tt-minlag+1,1)=...
                sqrt(sum((pstep(trg-tt,xidx)-tresp(trg)).^2));
            ttEE(xidx,tt-minlag+1,2)=...
                sqrt(sum((-pstep(trg-tt,xidx)-tresp(trg)).^2));
         end
      end
      maxh=min(find(abs(ttEE(:))==min(abs(ttEE(:)))));
      [xx,tt0,stepsign]=ind2sub([spacecount tbincount 2],maxh);
      tt=tt0+minlag-1;
      
      if stepsign==1,
         coeff=stepsize;
      else
         coeff=-stepsize;
      end
      
      % if that's bigger than regression result, take regression
      % result
      %maxcoeff=ttSR(xx,tt0)./sS(xx).^2./2;
      %if abs(coeff)>abs(maxcoeff),
      %   coeff=maxcoeff.*sign(coeff);
      %end
   end
   
   % predict cross-val for subtraction
   tridx=[1:length(tresp)]';
   tridx=tridx(find(tridx-tt>0 & tridx-tt<=length(tresp) & ...
                    ~isnan(tresp(tridx))));
   
   rpred(tridx)=coeff.*x(tridx-tt,xx);
   cumpred=cumpred+rpred;
   
   ff=find(~isnan(rpred) & ~isnan(tresp));
   otmag=nanstd(tresp(ff));
   ntmag=nanstd(tresp(ff)-rpred(ff));
   fprintf('step=%d(%d): xx=%d tt=%d ss=%d EE=%0.2f step=%0.3f rvar=%.4f->%.4f\n',...
           stepidx,forward*2-1,xx,tt,stepsign,ttEE(xx,tt0,stepsign),coeff,otmag,ntmag);
   if otmag<ntmag,
      disp('backwards??');
      %keyboard;
   end
   
   tresp(ff)=tresp(ff)-rpred(ff);
   
   % update parameters
   h(xx,tt0)=h(xx,tt0)+coeff;
   
   if stepidx==1,
      % initialize tolerance
      lambda=(sqrt(sum(y(ff).^2))-sqrt(sum((y(ff)-cumpred(ff)).^2)))./stepsize;
   elseif forward,
      lambdanew=(sqrt(sum((tresp(ff)+rpred(ff)).^2))-...
                 sqrt(sum(tresp(ff).^2)))./stepsize;
      lambda=min(lambda,lambdanew);
   end
   if lambda<0,
      break
   end
end
