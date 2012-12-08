% function [xc,pred,resp]=xcgetpreds(cellid,batch, 
%                                nlidx[=nlidxsave],cnfidx[=withinclass])
%
% load results from resfile in sRunData and return pred results
%
% xc correlation coefficient between non-nan time bins
% pred - predicted response
% resp - observed response
%
function [xc,pred,resp]=xcgetpreds(runidx,batch,nlidx0,cnfidx0)

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
   xc=nan;
   pred=[];
   resp=[];
   return
end


resfile=[rundata(1).respath,rundata(1).resfile,'.gz'];
if strcmp(resfile(end-2:end),'.gz'),
   zload(resfile);
else
   load(resfile);
end

if ~exist('nlidx0','var'),
   nlidx0=params.nlidxsave;
end

if ~exist('cnfidx0','var'),
   cnfidx0=1;
   if isfield(params,'predbatch') & isfield(params,'batch') & ...
         ~isempty(params.batch),
      cnfidx0=find(cat(1,params.predbatch{:})==params.batch);
   end
end

fprintf('cell %s batch %d nlidx=%d cnfidx=%d\n',...
        cellid,batch,nlidx0,cnfidx0);
fprintf('resfile=%s\n',resfile);

resp=predres(cnfidx0).act_resp{1};
pred=predres(cnfidx0).mod_psth{1}(:,nlidx0);

xc=predxc(cnfidx0,nlidx0);
keyboard

% figure out pred files for the current batch
[pcellfiledata,ptimes,pbatchdata]=...
    cellfiletimes(params.cellid,params.predbatch{cnfidx0},1);

% does this cell have data for batchid=predidx?
if length(pcellfiledata)>0,
   predparams=params;
   predparams.stimfiles={};
   predparams.respfiles={};
   predparams.stimcrfs=[];
   for ii=1:length(pcellfiledata),
      predparams.respfiles{ii}=[pcellfiledata(ii).path,...
                    pcellfiledata(ii).respfile];
   end
   
   tpredstartframe=ptimes(3).start;
   tpredstopframe=ptimes(3).stop;
   tpredfile=ptimes(3).fileidx;
   
   tresp=respload(predparams.respfiles{tpredfile(1)},...
                  'r',1,1,0);
   if size(tresp,2)>1,
      tresp=tresp(tpredstartframe:tpredstopframe,2:end);
      tresp=compact_raster_matrix3(tresp);
      
      firstwithnans=min([find(sum(isnan(tresp))>0) size(tresp,2)+1]);
      
      if firstwithnans>2,
         tresp=tresp(:,1:firstwithnans-1);
      end
   else
      tresp=tresp(tpredstartframe:tpredstopframe);
   end
   
   repcount=size(tresp,2);
   fprintf('predidx=%d (bat %d) (%d reps): ceiling... ',...
           cnfidx0,params.predbatch{cnfidx0},repcount);
   
   gg=find(~isnan(tresp(:,1)));
   [mu,alpha,beta]=reversepoisson(tresp(gg,1));
   rmax=singletrialceiling(tresp(gg,1),alpha,beta);
   %figure(gcf);
   drawnow
   
   xct=zeros(repcount,1);
   fprintf(' a=%.3f, b=%.3f\n',alpha,beta);
   for nn=1:size(strf,1),
      tpred=predres(cnfidx0).mod_psth{1}(:,1,nn);
      
      % deal with extra tails on tpred added by xcloadstimresp 
      if size(tresp,1)<length(tpred) & ptimes(3).start>params.maxlag(2),
         tpred=tpred(params.maxlag(2)+1:end);
      end
      drawnow
      if size(tresp,1)<length(tpred),
         tpred=tpred(1:size(tresp,1));
      end
      drawnow
      if size(tresp,1)>length(tpred)
         tresp=tresp(1:end+params.maxlag(1),:);
      end

      for xcidx=1:repcount,
         if length(tpred)~=length(tresp(:,xcidx)),
            keyboard
         end
         gg=find(~isnan(tpred) & ~isnan(tresp(:,xcidx)));
         xct(xcidx)=xcov(tpred(gg),tresp(gg,xcidx),0,'coeff');
      end
      
      projr=sqrt(mean(xct.^2)./rmax.^2);
      fprintf(' nl=%d: %.3f/%.3f ---> %.3f (vs. %.3f)\n',...
              nn,sqrt(mean(xct.^2)),rmax,projr,...
              predres(cnfidx0).predxc(nn));
      
      pred=tpred;
      resp=tresp;
      
   end
else
   % no val data for this batch
end



%            cd /auto/k2/share/matlab/wafting/
%            addmikiepath
%            pCCmetric=pCC(tpred(gg),sqrt(mean(xct.^2)));
%            fprintf(' pCCmetric=%.3f\n',pCCmetric);
