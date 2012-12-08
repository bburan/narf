function t0=find_raw_nl2(linpred,resp,stim,maxlag);

if ~exist('stim','var');
   verbose=0;
elseif length(stim)==1,
    verbose=1;
elseif size(stim,1)<size(linpred,1),
    verbose=2;
    eigs=stim;
else
    verbose=3;
end
if ~exist('maxlag','var'),
    maxlag=12;
end

dimcount=size(linpred,2);
if dimcount>2,
   warning('find_raw_nl2: only using first two dimensions');
end
d1=1;
d2=2;

gidx=find(~isnan(resp) & ~isnan(linpred(:,1)));
bincount=15;
if verbose>2,
    mstim=zeros(maxlag,bincount,bincount);
end
pp=zeros(bincount,2);
rr2=zeros(bincount,bincount);
[ss,si1]=sort(linpred(gidx,d1));
edges1=round(linspace(1,length(si1)+1,bincount+1));
[ss,si2]=sort(linpred(gidx,d2));
edges2=round(linspace(1,length(si2)+1,bincount+1));
for bb1=1:bincount,
   pp(bb1,1)=mean(linpred(gidx(si1(edges1(bb1):(edges1(bb1+1)-1))),d1));
   pp(bb1,2)=mean(linpred(gidx(si2(edges2(bb1):(edges2(bb1+1)-1))),d2));
   
   for bb2=1:bincount,
      tt=intersect(si1(edges1(bb1):(edges1(bb1+1)-1)),si2(edges2(bb2):(edges2(bb2+1)-1)));
      if length(tt)>0,
         rr2(bb1,bb2)=mean(resp(gidx(tt)));
         rn(bb1,bb2)=length(tt);
         if verbose>2,
             for jj=1:maxlag,
                 tt=tt(tt>=maxlag);
                 mstim(jj,bb1,bb2)=mean(stim(tt-jj+1));
             end
         end
         
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
   tpp=sub2ind(size(rr2),xi(nn),yi(nn));
   rr3(jj,kk)=mean(rr2(tpp));
end
rr3=gsmooth(rr3,[1 1]);

if verbose==1,
   [xi,yi]=meshgrid(linspace(pp(1,2),pp(end,2),bincount*2),...
                    linspace(pp(1,1),pp(end,1),bincount*2));
   z=interp2(pp(:,2),pp(:,1),rr3,xi,yi,'bilinear');
   zn=interp2(pp(:,2),pp(:,1),rn,xi,yi,'bilinear');
   imagesc(pp(:,2),pp(:,1),z*100);
   hold on
   contour(xi,yi,zn,[5 5],'k','LineWidth',1);
   hold off
   xlabel('p2');
   ylabel('p1');
   axis xy
   colorbar
    
elseif verbose==2,
    figure;
    
    
    [m1,m2]=ind2sub(size(rr3),find(rr3==max(rr3(:))));
    optstim=fliplr(pp(m1,1).*eigs(1,:)+pp(m2,2).*eigs(2,:));
    hs=[];
    mx=0;mn=0;
    for ii=1:3,
        for jj=1:3,
            if ii==2 && jj==2,
               subplot(3,3,ii+(jj-1)*3);
                [xi,yi]=meshgrid(linspace(pp(1,2),pp(end,2),bincount*2),...
                                 linspace(pp(1,1),pp(end,1),bincount*2));
                z=interp2(pp(:,2),pp(:,1),rr3,xi,yi,'bilinear');
                zn=interp2(pp(:,2),pp(:,1),rn,xi,yi,'bilinear');
                imagesc(pp(:,2),pp(:,1),z*100);
                hold on
                contour(xi,yi,zn,[5 5],'k','LineWidth',1);
                hold off
                axis xy
                colorbar
            else
               hs=cat(1,hs,subplot(3,3,ii+(jj-1)*3));
               if jj==2,
                  bb1=find(abs(pp(:,1))==min(abs(pp(:,1))));
               else
                  bb1=round((4-jj)/3*bincount);
               end
               if ii==2,
                  bb2=find(abs(pp(:,2))==min(abs(pp(:,2))));
               else
                  bb2=round(ii/3*bincount);
               end
               tstim=pp(bb1,1).*eigs(1,:)+pp(bb2,2).*eigs(2,:);
               tstim=fliplr(tstim);
               
               plot(optstim,'g--','LineWidth',2);
               hold on
               plot([1 size(eigs,2)],[0 0],'k--');
               plot(tstim,'k','LineWidth',2);
               hold off
               aa=axis;
               mx=max(mx,aa(4));
               mn=min(mn,aa(3));
               axis off;
               
                %axis([0 maxlag+1 ]);
            end
        end
    end
    
    for ii=1:length(hs),
       axes(hs(ii));
       axis([0 size(eigs,2)+1 mn mx]); 
    end
    
    fullpage landscape
elseif verbose==3,
    figure;
    
    for ii=1:3,
        for jj=1:3,
            subplot(3,3,ii+(jj-1)*3);
            if ii==2 && jj==2,
                [xi,yi]=meshgrid(linspace(pp(1,2),pp(end,2),bincount*2),...
                                 linspace(pp(1,1),pp(end,1),bincount*2));
                z=interp2(pp(:,2),pp(:,1),rr3,xi,yi);
                imagesc(pp(:,2),pp(:,1),z);
                xlabel('p2');
                ylabel('p1');
            else
                if jj==1,
                    bb1=1:3;
                elseif jj==2,
                    bb1=round(bincount/2)+(-1:1);
                else
                    bb1=(bincount-2):bincount;
                end
                if ii==1,
                    bb2=1:3;
                elseif ii==2,
                    bb2=round(bincount/2)+(-1:1);
                else
                    bb2=(bincount-2):bincount;
                end
                plot(mean(mean(mstim(:,bb1,bb2),2),3));
                axis([0 maxlag+1 min(mstim(:)) max(mstim(:))]);
            end
        end
    end
    fullpage landscape
    
end

t0={pp,rr3};

