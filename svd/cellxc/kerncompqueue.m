% kerncompqueue.m:  wrapper for kerncomp4.m attentional modulation testing
%
% basic flow: 
% 1. read next un-run "kerncomp" entry in tQueue
% 2. load up parameters and call kerncomp.m
%
% created SVD 11/10/02  branched off of cellxcqueue.m
%


dbopen;

disp('Clearing all variables from memory...');
if debugcheck,
   clear all
   debugon;
else
   clear all
   debugoff;
end

for ii=1:0,
   figure(ii);
end
drawnow;

% REDOXC: (1) do the xc and decorr; (0) LOAD RESULTS FROM LAST RUN (0)
REDOXC=1;
rand('state',sum(100*clock));

global BATQUEUEID
BATQUEUEID=0;
while BATQUEUEID>=0,
   
   BATQUEUEID=dbgetnextqueue('kerncomp',BATQUEUEID);
   
   if BATQUEUEID<0,
      % no more entries in queue
      return
   end
   
   sql=['SELECT * from tQueue WHERE id=',num2str(BATQUEUEID)];
   queuedata=mysql(sql);
   
   sql=['SELECT * from sRunData WHERE id=',num2str(queuedata.rundataid)];
   rundata=mysql(sql);
   cellid=rundata.cellid;
   
   [cellfiledata,times,batchdata]=cellfiletimes(cellid,rundata.batch);
   
   starttimes=times(1).start;
   stoptimes=times(1).stop;
   fitfile=times(2).fileidx;
   fitstartframe=times(2).start;
   fitstopframe=times(2).stop;
   predfile=times(3).fileidx;
   predstartframe=times(3).start;
   predstopframe=times(3).stop;
   
   stimfiles={};
   respfiles={};
   resplens=zeros(length(cellfiledata),1);
   stimpix=zeros(length(cellfiledata),1);
   stimcrfs=zeros(length(cellfiledata),1);
   totlen=0;
   for ii=1:length(cellfiledata),
      stimfiles{ii}=[cellfiledata(ii).stimpath,cellfiledata(ii).stimfile];
      respfiles{ii}=[cellfiledata(ii).path,cellfiledata(ii).respfile];
      resplens(ii)=cellfiledata(ii).resplen;
      tstimpix=strsep(cellfiledata(ii).stimiconside,',');
      stimpix(ii)=tstimpix{1};
      stimcrfs(ii)=cellfiledata(ii).stimfilecrf;
   end
   
   % predbatch--batches containing other stim classes
   if isempty(batchdata.predbatch),
      predbatch={batchdata.id};
   else
      predbatch=strsep(batchdata.predbatch,',');
   end
   
   stimfmtcode=batchdata.stimfmtcode;
   respfmtcode=batchdata.respfmtcode;
   
   % resp load cmd. format: resploadcmd(respfile,p1,p2,...)
   resploadcmd=batchdata.resploadcmd;
   resploadparms=strsep(batchdata.resploadparms,',');
   
   % resp filter - eg, resample or smooth 
   respfiltercmd=batchdata.respfiltercmd;
   respfilterparms=strsep(batchdata.respfilterparms,',');
   
   % stim load cmd. should be of format
   %   stimloadcmd(stimfile,startframe,endframe,p1,p2,...)
   stimloadcmd=batchdata.stimloadcmd;
   stimloadparms=strsep(batchdata.stimloadparms,',');
   
   % stim filter - eg, convert to phase-sep fourier domain.
   stimfiltercmd=batchdata.stimfiltercmd;
   stimfilterparms=strsep(batchdata.stimfilterparms,',');
   stimwindowcrf=batchdata.stimwindowcrf;
   
   kernfmt=batchdata.kernfmt;
   
   maxlag=[batchdata.minlag batchdata.maxlag];
   resampfmt=batchdata.resampfmt;     % ie bootstrapping or
   resampcount=batchdata.resampcount; % number of resampled kernels
                                      % something else
   dotSA=batchdata.decorrtime;
   dosSA=batchdata.decorrspace;
   
   neigs=0:batchdata.sfsstep:(batchdata.sfsstep*batchdata.sfscount-1);
   sfscount=batchdata.sfscount;
   sffiltsigma=batchdata.sffiltsigma;
   nloutparm=batchdata.nloutparm;
   
   % preload old results for saving if not recalculated
   sql=['SELECT * FROM sResults WHERE runid=',num2str(rundata.id)];
   resdata=mysql(sql);
   if length(resdata)>0,
      eval(char(resdata.matstr));
      predxc=preddata(rundata.id).predxc;
      if isfield(preddata(rundata.id),'perf'),
         %perf=preddata(rundata.id).perf;
         %perfatt=preddata(rundata.id).perfatt;
         %perfatt0=preddata(rundata.id).perfatt0;
         %pstrf=preddata(rundata.id).pstrf;
         %pstrfatt=preddata(rundata.id).pstrfatt;
         %nondampedcount=preddata(rundata.id).nondampedcount;
      else
         perf=[];
         perfatt=[];
         perfatt0=[];
         pstrf=[];
         pstrfatt=[];
         nondampedcount=[];
      end
   else
      predxc=nan;
   end

   postanal='kerncomp';
   
   % this is for the old paired attention comparison. may want
   % to go back 
   %resampcount=1;
   %noisecount=300;
   %kerncomp3;  % kerncomp or kerncomp2 or kerncomp3 or kerncomp4?
   
   noisecount=500;
   kerncomp4;  % kerncomp or kerncomp2 or kerncomp3 or kerncomp4?
   
   % save the results
   daterun=date;
   skernfile=rundata.kernfile;
   sfullfile=[rundata.respath,skernfile];
   disp(sprintf('Saving kerncomp results to %s...',skernfile));
   
   % then save everything else (inclusive to avoid future screw-ups
   % involving accidentally not saving new important stuff)
   tpath=getenv('TMP');
   save([tpath,'/',skernfile]);
   
   % compress and copy over network if necessary
   if checkbic,
      unix(['gzip ',tpath,'/',skernfile]);
      if checkbic==1,
         unix(['jsend ',tpath,'/',skernfile,'.gz ',sfullfile,'.gz']);
      elseif checkbic==2,
         unix(['bsend ',tpath,'/',skernfile,'.gz ',sfullfile,'.gz']);
      end
      delete([tpath,'/',skernfile,'.gz']);
   else
      eval(['!gzip -cf ',tpath,'/',skernfile,'>',sfullfile,'.gz']);
      delete([tpath,'/',skernfile]);
   end
   
   % construct a string to dump to the results table for later summarizing
   sres=sprintf('preddata(%d).cellid=''%s'';',rundata.id,cellid);
   sres=sprintf('%s preddata(%d).predxc=%s;',...
                sres,rundata.id,mat2string(predxc));
   sres=sprintf('%s preddata(%d).pxc=%s;',...
                sres,rundata.id,mat2string(pxc));
   sres=sprintf('%s preddata(%d).pxcp=%s;',...
                sres,rundata.id,mat2string(pxcp));
   sres=sprintf('%s preddata(%d).nondampedcount=%s;',...
                sres,rundata.id,mat2string(nondampedcount));
   sres=sprintf('%s preddata(%d).ncount=%s;',...
                sres,rundata.id,mat2string(ncount));
   
   % determine whether there's already a results record and delete
   % it so that it can be replaced.  maybe this should be
   % superceded by an archiving scheme someday?
   sql=['SELECT * FROM sResults',...
        ' WHERE runid=',num2str(rundata.id)];
   resdata=mysql(sql);
   if length(resdata)>0,
      mysql(['DELETE FROM sResults WHERE id=',num2str(resdata(1).id)]);
   end
   sqlinsert('sResults',...
             'runid',rundata.id,...
             'batch',rundata.batch,...
             'matstr',sres);
   
   % record that we're done with this queue entry
   dbsetqueue(BATQUEUEID,1000,1);
end






