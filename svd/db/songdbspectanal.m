
dbopen;

minevals=8;
rtcount=6;
spcount=5;
minuserevals=50;

[evalmtx,users,songid]=songevalmtx(minevals);
fullsongid=songid;
fullevalmtx=evalmtx;
ecount=sum(~isnan(evalmtx));
keepuseridx=find(ecount>minuserevals);
fullevalmtx=fullevalmtx(:,keepuseridx);
tevalmtx=evalmtx(:,keepuseridx);
users={users{keepuseridx}};
ucount=size(tevalmtx,2);
fprintf('minuserevals %d leaves %d users with ratings for stats\n',...
        minuserevals,ucount);

sql=['SELECT * FROM dEigProj order by songid'];
eigdata=mysql(sql);
songcount=length(eigdata);

ecell=struct2cell(eigdata);

rtmtx=zeros(songcount,rtcount);
spmtx=zeros(songcount,spcount);

esongid=cat(1,ecell{1,:});
ratedsongs=find(ismember(esongid,songid));
for ii=1:rtcount,
   rtmtx(:,ii)=cat(1,ecell{ii+1,:});
end
for ii=1:spcount,
   spmtx(:,ii)=cat(1,ecell{ii+11,:});
end
rspmtx=spmtx(ratedsongs,:);


d=[tevalmtx rspmtx];
mm=nanmean(d);
ss=nanstd(d);
d=(d-repmat(mm,[length(d),1])) ./ repmat(ss,[length(d),1]);

cc=zeros(ucount+spcount);
n=zeros(ucount+spcount);
for u1=1:(ucount+spcount),
   for u2=u1:(ucount+spcount),
      n(u1,u2)=sum(~isnan(d(:,u1).*d(:,u2)));
      n(u2,u1)=n(u1,u2);
      cc(u1,u2)=nanmean(d(:,u1).*d(:,u2));
      cc(u2,u1)=cc(u1,u2);
   end
end
cc(isnan(cc))=0;
ccns=cc(1:ucount,1:ucount);

[u,s,v]=svd(cc);
[uns,sns,vns]=svd(ccns);

user='david';
useridx=find(strcmp(user,users));
eigmax=20;

disp('finding optimal number of eigenvectors...');
ee=zeros(eigmax,2,ucount);
for useridx=1:ucount,
   testids=find(~isnan(d(:,useridx)));
   for svdidx=1:eigmax,
      zmtx=d;
      zmtx(find(isnan(zmtx)))=0;
      zmtx(:,useridx)=0;
      emtx=zmtx(testids,:)*u(:,1:svdidx);
      tmtx=emtx*u(:,1:svdidx)';
      
      ee(svdidx,2,useridx)=xcov(tmtx(:,useridx),...
                        d(testids,useridx),0,'coeff');
      
      % exclude spectral analysis dimensions
      emtx=zmtx(testids,1:ucount)*uns(:,1:svdidx);
      tmtx=emtx*uns(:,1:svdidx)';
      
      ee(svdidx,1,useridx)=xcov(tmtx(:,useridx),...
                        d(testids,useridx),0,'coeff');
      
   end
end










evalcount=size(d,1);
jackcount=10;
jstep=evalcount./jackcount;
E=zeros(spcount,jackcount);

for jj=1:jackcount,
   validx=round(((jj-1).*jstep+1):jj*jstep);
   estidx=setdiff(1:evalcount,validx);
   
   estin=d(estidx,ucount+(1:spcount));
   estout=d(estidx,1:ucount);
   valin=d(validx,ucount+(1:spcount));
   valout=d(validx,1:ucount);
   
   sta=zeros(ucount,spcount);
   for u1=1:(ucount),
      ff=find(~isnan(estout(:,u1)));
      sta(u1,:)=estout(ff,u1)'*estin(ff,:)./length(ff).*length(estin);
   end
   
   for ss=1:spcount,
      
      C=estin(:,1:ss)'*estin(:,1:ss);
      h=sta(:,1:ss) * C^-1;
      
      predout=valin(:,1:ss)*h';
      
      ff=find(~isnan(valout));
      E(ss,jj)=mean((predout(ff)-valout(ff)).^2);
   end
   
end
