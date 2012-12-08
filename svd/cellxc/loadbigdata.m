% function [res,bigdata]=loadbigdata(batches,nloutidx,cnfidx);
function [res,bigdata]=loadbigdata(batches,nloutidx,cnfidx);

% load pred results for each batch
expcount=size(batches,1);
domcount=size(batches,2);
cellmaster=zeros(expcount,domcount,200);
totlen=zeros(expcount,domcount,200);

ssidx=[];
for expidx=1:expcount,
   for domidx=1:domcount,
      
      sql=['SELECT sResults.*,singleid,cellid',...
           ' FROM sResults INNER JOIN sRunData',...
           ' ON sResults.runid=sRunData.id',...
           ' WHERE sResults.batch=',num2str(batches(expidx,domidx)),...
           ' ORDER BY sRunData.cellid'];
%           ' AND singleid>0 ORDER BY sRunData.cellid'];
%           ' AND not(sRunData.cellid like "model%")',...
      resdata=mysql(sql);

      sid=cat(2,resdata.singleid);
      ssidx=cat(2,ssidx,setdiff(sid,[0 ssidx]));
      
      for runidx=find(sid~=0),
         singlematch=find(sid(runidx)==ssidx);
         
         cellmaster(expidx,domidx,singlematch)=1;
         
         % record prediction results for current run
         eval(char(resdata(runidx).matstr));
         
         runid=resdata(runidx).runid;
         preddata(runid).runid=runid;
         
         %figure out interesting parameters of run (eg, total file length)
         cellid=resdata(runidx).cellid;
         cellfiledata=cellfiletimes(cellid,batches(expidx,domidx));
         
         preddata(runid).totlen=sum(cat(1,cellfiledata.resplen));
         preddata(runid).spikes=sum(cat(1,cellfiledata.spikes));
         if preddata(runid).totlen>0,
            preddata(runid).meanrate=...
                preddata(runid).spikes./preddata(runid).totlen .* 1000/14;
         else
            preddata(runid).meanrate=0;
         end
         
         if ~isfield(preddata(runid),'predp'),
            preddata(runid).predp=preddata(runid).predxc.*nan;
         end
         if ~isfield(preddata(runid),'prederr'),
            preddata(runid).prederr=preddata(runid).predxc.*nan;
         end
         if ~isfield(preddata(runid),'predinf'),
            preddata(runid).predinf=preddata(runid).predxc.*nan;
         end
         if ~isfield(preddata(runid),'predone'),
            preddata(runid).predone=preddata(runid).predxc.*nan;
         end
         if ~isfield(preddata(runid),'predmse'),
            preddata(runid).predmse=preddata(runid).predxc.*nan;
         end
         if ~isfield(preddata(runid),'predfix'),
            preddata(runid).predfix=preddata(runid).predxc.*nan;
         end
         if ~isfield(preddata(runid),'predfixerr'),
            preddata(runid).predfixerr=preddata(runid).predxc.*nan;
         end
         
         % this could be made slicker!
         bigdata(expidx,domidx,singlematch)=preddata(runid);
         bigdata(expidx,domidx,singlematch).cellid=cellid;
      end
   end
end

% don't need this any more
clear preddata

fprintf('bigdata matrix loaded expcount=%d domcount=%d\n',expcount,domcount);

% preddata(ttt(1)).predxc: batchcount x strfcount x attcount x latcount

goodcells=find(squeeze(sum(sum(cellmaster,1),2)));
cellcount=length(goodcells);
cellmaster=cellmaster(:,:,goodcells);
bigdata=bigdata(:,:,goodcells);

ii=min(find(cellmaster(1,1,:)));

if ~exist('cnfidx','var'),
   cnfidx=repmat(1:size(bigdata(1,1,ii).predxc,1),[length(batches) 1]);
end
cnfcount=size(cnfidx,2);

if ~exist('nloutidx','var'),
   nloutidx=repmat(1:size(bigdata(1,1,ii).predxc,2),[length(batches) 1]);
end
nlcount=size(nloutidx,2);
respcount=size(bigdata(1,1,ii).predxc,4);

% predxc: cell X expclass X in-model X cnfclass X out-model
predxc=ones(cellcount,expcount,domcount,cnfcount,nlcount,respcount).*nan;
predinf=ones(cellcount,expcount,domcount,cnfcount,nlcount,respcount).*nan;
predp=ones(cellcount,expcount,domcount,cnfcount,nlcount,respcount);
prederr=ones(cellcount,expcount,domcount,cnfcount,nlcount,respcount).*nan;
predfix=ones(cellcount,expcount,domcount,cnfcount,nlcount,respcount).*nan;
predfixerr=ones(cellcount,expcount,domcount,cnfcount,nlcount,respcount).*nan;
predmse=ones(cellcount,expcount,domcount,cnfcount,nlcount,respcount).*nan;
expxc=zeros(cellcount,expcount,domcount,nlcount,respcount);
sigfit=ones(cellcount,expcount,domcount,nlcount,respcount);
sfsfit=ones(cellcount,expcount,domcount,nlcount,respcount);
meanrate=zeros(cellcount,expcount).*nan;
totlen=ones(cellcount,expcount).*nan;
spikes=ones(cellcount,expcount).*nan;
rundataid=ones(cellcount,expcount).*nan;
celllist=cell(cellcount,1);

for expidx=1:expcount,
   for domidx=1:domcount,
      tgoodidx=find(cellmaster(expidx,domidx,:));
      tgoodidx=tgoodidx(:)';
      
      for cellidx=tgoodidx,
         celllist{cellidx}=bigdata(expidx,domidx,cellidx).cellid;
         predxc(cellidx,expidx,domidx,:,:,:)=...
             reshape(bigdata(expidx,domidx,cellidx).predxc(...
                cnfidx(expidx,:),nloutidx(expidx,:),1,:),...
                     [1 1 1 cnfcount nlcount respcount]);
         if prod(size(bigdata(expidx,domidx,cellidx).predinf))>0,
            predinf(cellidx,expidx,domidx,:,:,:)=...
                reshape(bigdata(expidx,domidx,cellidx).predinf(...
                   cnfidx(expidx,:),nloutidx(expidx,:),1,:),...
                        [1 1 1 cnfcount nlcount respcount]);
         end
         if prod(size(bigdata(expidx,domidx,cellidx).predp))>0,
            predp(cellidx,expidx,domidx,:,:,:)=...
                reshape(bigdata(expidx,domidx,cellidx).predp(...
                   cnfidx(expidx,:),nloutidx(expidx,:),1,:),...
                        [1 1 1 cnfcount nlcount respcount]);
         end
         if prod(size(bigdata(expidx,domidx,cellidx).prederr))>0,
            prederr(cellidx,expidx,domidx,:,:,:)=...
                reshape(bigdata(expidx,domidx,cellidx).prederr(...
                   cnfidx(expidx,:),nloutidx(expidx,:),1,:),...
                        [1 1 1 cnfcount nlcount respcount]);
         end
         if prod(size(bigdata(expidx,domidx,cellidx).expxc))>0,
            expxc(cellidx,expidx,domidx,:,:)=...
                reshape(bigdata(expidx,domidx,cellidx).expxc(...
                   :,nloutidx(expidx,:),1)',...
                        [1 1 1 nlcount respcount]);
         end
         if isfield(bigdata(expidx,domidx,cellidx),'predfix') & ...
               length(bigdata(expidx,domidx,cellidx).predfix) > 0,
            predfix(cellidx,expidx,domidx,:,:,:)=...
                reshape(bigdata(expidx,domidx,cellidx).predfix(...
                   cnfidx(expidx,:),nloutidx(expidx,:),1,:),...
                        [1 1 1 cnfcount nlcount respcount]);
            predfixerr(cellidx,expidx,domidx,:,:,:)=...
                reshape(bigdata(expidx,domidx,cellidx).predfixerr(...
                   cnfidx(expidx,:),nloutidx(expidx,:),1,:),...
                        [1 1 1 cnfcount nlcount respcount]);
         end
         
         %predmse(cellidx,expidx,domidx,:,:,:)=...
         %    reshape(bigdata(expidx,domidx,cellidx).predmse(...
         %       cnfidx(expidx,:),nloutidx(expidx,:),1,:),...
         %            [1 1 1 cnfcount nlcount respcount]);
         
         if domidx==1,
            totlen(cellidx,expidx)=bigdata(expidx,domidx,cellidx).totlen;
            spikes(cellidx,expidx)=bigdata(expidx,domidx,cellidx).spikes;
            meanrate(cellidx,expidx)=bigdata(expidx,domidx,cellidx).meanrate;
         end
      end
   end
end

[celllist,si]=sort(celllist);

res.celllist=celllist;
res.predxc=predxc(si,:,:,:,:);
res.predinf=predinf(si,:,:,:,:);
res.predp=predp(si,:,:,:,:);
res.prederr=prederr(si,:,:,:,:);
res.predmse=predmse(si,:,:,:,:);
res.predfix=predfix(si,:,:,:,:);
res.predfixerr=predfixerr(si,:,:,:,:);
res.expxc=expxc(si,:,:,:,:);
res.totlen=totlen(si,:,:,:,:);
res.spikes=spikes(si,:,:,:,:);
res.meanrate=meanrate(si,:,:,:,:);

