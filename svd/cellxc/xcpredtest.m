function [predxc,prederr]=xcpredtest(cellid,batch,predbatch);

dbopen;
goodbatch=zeros(1,length(batch));
sql=['SELECT * from sRunData WHERE cellid="',cellid,'"',...
     ' AND batch=',num2str(batch)];
rundata=mysql(sql);

if length(rundata)==0,
   disp('no entry found in db!');
   if nargout>0,
      r=0;
   end
   return
end

resfile=[rundata.respath,rundata.resfile,'.gz'];
fprintf('loading: %s\n',resfile);
zload(resfile);

global BATQUEUEID
BATQUEUEID=[];

batchcount=size(strf,1);

if ~exist('predbatch','var'),
   predbatch=1;
end

[cellfiledata,times]=cellfiletimes(cellid,sparams.estbatch(predbatch));

predparams=params;
predparams.stimfiles={};
predparams.respfiles={};
predparams.stimcrfs=[];

fileidx=[];
starttimes=[]; stoptimes=[];
for ii=1:length(cellfiledata),
   
   predparams.stimfiles{ii}=[cellfiledata(ii).stimpath,...
                    cellfiledata(ii).stimfile];
   predparams.respfiles{ii}=[cellfiledata(ii).path,...
                    cellfiledata(ii).respfile];
   if cellfiledata(ii).stimfilecrf>0,
      predparams.stimcrfs(ii)=cellfiledata(ii).stimfilecrf;
   else
      tstimpix=strsep(cellfiledata(ii).stimiconside,',');
      if length(tstimpix)>0,
         tstimpix=tstimpix{1};
      end
      sql=['SELECT * FROM gCellMaster WHERE cellid="',...
           params.cellid,'"'];
      celldata=mysql(sql);
      predparams.stimcrfs(ii)=tstimpix./celldata.rfsize;
   end
   
   if cellfiledata(ii).repcount>=6,
      fileidx=[fileidx;ii];
      starttimes=[starttimes;1];
      stoptimes=[stoptimes;cellfiledata(ii).resplen];
   end
end

if (length(fileidx)==0 | sum(stoptimes)<1000) & ...
       (length(fileidx)==0 |fileidx(1)~=1) ...
   disp('not much multi-trial data. also using first file');
   
   fileidx=[1;fileidx];
   starttimes=[1; starttimes];
   stoptimes=[cellfiledata(1).resplen;stoptimes];
end

clear cdata

[cdata.stim,cdata.resp]=xcloadstimresp(fileidx,starttimes,stoptimes,predparams);
tpredres=xcval(strf(:),predparams,cdata);
[tpredres.predxc tpredres.predfix]
predxc=tpredres.predxc;
prederr=tpredres.prederr;

