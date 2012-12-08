% function [ccr,ccrM,M,vpredxc,xcvstrials]=xcfredfrac(runidx,batch);
%
% figure out asymptotic prediction scores based on fits to single/finite
% trial predictions. use the hsu and theunissen algorithm
%
% the idea is that noise goes as 1/N * sigma(neuron)
% that is, it shrinks linearly as one over the number of
% trials. when 1/N goes to 0 (ie, infinite trials), the noise
% should be zero. any residual error is the quality of the model
% doing the prediction
%
% returns:
% ccr - ratio of single trial cc^2 to expected noise-free cc^2
% ccrM - ratio of single trial cc^2 to M trial cc^2
%
function [ccr,ccrM,M,vpredxc,xcvstrials]=xcfredfrac(runidx,batch);

disp('xcfredfrac.m:');

global BATQUEUEID

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

outfile=[rundata(1).respath,rundata(1).resfile,'.gz'];
zload(outfile);

BATQUEUEID=[];
ccr=ones(1,length(params.predbatch)).*nan;
ccrM=ones(1,length(params.predbatch)).*nan;
M=zeros(1,length(params.predbatch));

nlidx=params.nlidxsave;
if params.fitboot==1,
   bootcount=1;
end


for predidx=1:length(params.predbatch),
   xcvstrials(predidx).m=nan;
   xcvstrials(predidx).s=nan;
   xcvstrials(predidx).cxy=nan;
   xcvstrials(predidx).n=nan;
   pdata(predidx).alive=1;
   
   fprintf('Computing fredfrac for cell %s batch %d.\n',...
           params.cellid,params.predbatch{predidx});
   
   % figure out pred files for the current batch
   [pcellfiledata,ptimes,pbatchdata]=...
       cellfiletimes(params.cellid,params.predbatch{predidx});
   
   % does this cell have data for batchid=predidx?
   if length(pcellfiledata)>0,
      predparams=params;
      predparams.stimfiles={};
      predparams.respfiles={};
      predparams.stimcrfs=[];
      for ii=1:length(pcellfiledata),
         predparams.stimfiles{ii}=[pcellfiledata(ii).stimpath,...
                    pcellfiledata(ii).stimfile];
         predparams.respfiles{ii}=[pcellfiledata(ii).path,...
                    pcellfiledata(ii).respfile];
         if pcellfiledata(ii).stimfilecrf>0,
            predparams.stimcrfs(ii)=pcellfiledata(ii).stimfilecrf;
         else
            tstimpix=strsep(pcellfiledata(ii).stimiconside,',');
            if length(tstimpix)>0,
               tstimpix=tstimpix{1};
            end
            sql=['SELECT * FROM gCellMaster WHERE cellid="',...
                 params.cellid,'"'];
            celldata=mysql(sql);
            predparams.stimcrfs(ii)=tstimpix./celldata.rfsize;
         end
      end
      
      tpredstartframe=ptimes(3).start;
      tpredstopframe=ptimes(3).stop;
      tpredfile=ptimes(3).fileidx;
      
      predparams.stimloadcmd=pbatchdata.stimloadcmd;
      predparams.stimloadparms=strsep(pbatchdata.stimloadparms,',');
      
      predparams.stimfiltercmd=pbatchdata.stimfiltercmd;
      predparams.stimfilterparms=strsep(pbatchdata.stimfilterparms,',');
      predparams.resploadcmd=pbatchdata.resploadcmd;
      predparams.resploadparms=strsep(pbatchdata.resploadparms,',');
      predparams.respfiltercmd=pbatchdata.respfiltercmd;
      predparams.respfilterparms=strsep(pbatchdata.respfilterparms,',');
      
      % do asymptotic thingy
      ro=respload(predparams.respfiles{tpredfile});
      ro=ro(tpredstartframe:tpredstopframe,:);
      r=compact_raster_matrix3(ro(:,2:end));      
      
      if size(r,1)>2 & size(r,2)>2,
         
         r2=r;
         
         okidx=find(~isnan(r2(:,1)));
         
         r2=r2(okidx,:);
         r2=r2(:,find(sum(isnan(r2),1)==0));
         fprintf('%d repeats of validation data\n',size(r2,2));
         
         if params.predsmoothsigma > 0 & params.respfmtcode==0,
            %fprintf('Smoothing actual response with optimal filter...\n');
            %pfilt=[1/9 2/9 1/3 2/9 1/9]';
            
            fprintf('(rr: smsig=%.1f)',params.predsmoothsigma);
            tt=(-10:10)';
            pfilt=exp(-tt.^2/(2*params.predsmoothsigma.^2))./...
                  (sqrt(2*pi)*params.predsmoothsigma);
            pfilt=pfilt(find(pfilt > 0.1*max(pfilt)));
            fr=r2;
            fr(isnan(fr))=0;
            fr=conv2(fr,pfilt,'same');
            r2(~isnan(r2))=fr(~isnan(r2));
         end
         fprintf('\n');
         
         [ccr(predidx),ccrM(predidx)]=frednoise(r2);
         M(predidx)=size(r2,2);
         
         ncount=round((M(predidx)+1)/2);
         
         curveres(predidx).m=zeros(ncount,bootcount);
         curveres(predidx).s=zeros(ncount,bootcount);
         curveres(predidx).n=zeros(ncount,bootcount);
         curveres(predidx).cxy=zeros(ncount,M(predidx),bootcount);
         
         for bootidx=1:bootcount,
            pp=predres(predidx).mod_psth{1}(okidx,:,nlidx+ ...
                                            (bootidx-1)*params.nloutparm);
            minpp=min(find(~isnan(pp)));
            pp=pp(minpp:end);
            r2=r2(1:length(pp),:);
            %ff=find(~isnan(pp));
            %xcov(pp(ff),mean(r2(ff,:),2),0,'coeff')
            [curveres(predidx).m(:,bootidx),...
             curveres(predidx).s(:,bootidx),...
             curveres(predidx).n(:,bootidx),...
             curveres(predidx).cxy(:,:,bootidx)]=predtrialcurve(pp,r2);
         end
         
         xcvstrials(predidx).m=mean(curveres(predidx).m,2);
         % kloodgy
         %xcvstrials(predidx).s=mean(curveres(predidx).s,2);
         % realistic, error between trials
         xcvstrials(predidx).s=std(curveres(predidx).m,0,2) ...
             .* sqrt(bootcount);
         xcvstrials(predidx).cxy=curveres(predidx).cxy;
         xcvstrials(predidx).n=curveres(predidx).n(:,1);
         
      else
         fprintf('only one trial in cnf data!\n');
         M(predidx)=1;
      end
      
   end
end

% save important results

