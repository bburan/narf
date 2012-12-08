% function r=xcsupptest(cellid,batch)
%
% load results from resfile in sRunData test suppressive output NL
%
% r=0 if no entries found in db, = [cclin, ccsupp] otherwise
%
function r=xcsupptest(batchid)

rundata=dbgetrundata(batchid);
r=zeros(length(rundata),2);

for ii=1:length(rundata);
   r(ii,:)=xcsupptest0(rundata(ii).cellid,batchid);
end

keyboard

function r=xcsupptest0(runidx,batch);

dbopen;
if isnumeric(runidx),
   sql=['SELECT * from sRunData WHERE id=',num2str(runidx)];
   rundata=mysql(sql);
   cellid=rundata.cellid;
else
   cellid=runidx;
   goodbatch=zeros(1,length(batch));
   sql=['SELECT * from sRunData WHERE cellid="',cellid,'"',...
        ' AND batch=',num2str(batch)];
   rundata=mysql(sql);
end

if length(rundata)==0,
   disp('no entries found in db!');
   if nargout>0,
      r=0;
   end
   return
end

resfile=[rundata(1).respath,rundata(1).resfile,'.gz'];

fprintf('%s.m: %s\n',mfilename,resfile);
zload(resfile);

nlidx=params.nlidxsave;
if nlidx>length(strf),
   nlidx=1;
end
fprintf('chose to show nlidx=%d\n',nlidx);

global BATQUEUEID
clear BATQUEUEID

xcloadfiles;

hf=strf(nlidx).h;

tstim=stim'-repmat(strf(nlidx).mS,[1 movlen]);
seplinpred=kernpredict(strf(nlidx).h,tstim,spacecount,0);
hp=find(strf(nlidx).hspace>0);
hn=find(strf(nlidx).hspace<0);
tr=[sum(seplinpred(:,hp),2) -sum(seplinpred(:,hn),2)];

rlen=min([3000 length(resp)]);
gidx=find(~isnan(resp(1:rlen)));
nlparms=fitsuppnorm(tr(gidx,:),resp(gidx));

[cstim,cresp]=xcloadstimresp(times(3).fileidx,times(3).start,...
                             times(3).stop,params);

pr1=xcpredict(strf(nlidx),cstim);

tstrf=strf(nlidx);
tstrf.nltype='suppnorm';
tstrf.nlparms=nlparms;
pr2=xcpredict(tstrf,cstim);

gidx=find(~isnan(cresp));

r=[xcov(pr1(gidx),cresp(gidx),0,'coeff'),xcov(pr2(gidx),cresp(gidx),0,'coeff')]



