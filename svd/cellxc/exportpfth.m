% function z=exportpfth(cellid,batch)
%
% load results from kernfile in sRunData and display attentional
% modulation info
%
% z=data from kernfile
%
function z=exportpfth(runidx,batch)

outpath='/auto/k1/david/tmp/';

dbopen;
if isnumeric(runidx),
   sql=['SELECT * from sRunData WHERE id=',num2str(runidx)];
   rundata=mysql(sql);
   batchcount=length(rundata);
   cellid=rundata.cellid;
else
   cellid=runidx;
   goodbatch=zeros(1,length(batch));
   batchcount=0;
   for ii=1:length(batch),
      sql=['SELECT * from sRunData WHERE cellid="',cellid,'"',...
           ' AND batch=',num2str(batch(ii))];
      trd=mysql(sql);
      if ~isempty(trd),
         goodbatch(ii)=1;
         batchcount=batchcount+1;
         rundata(batchcount)=trd;
      end
   end
   batch=batch(find(goodbatch));
end

if batchcount==0,
   disp('no entries found in db!');
   if nargout>0,
      r=0;
   end
   return
end

rcsetstrings;

for runidx=1:batchcount,
   [cellfiledata,times,batchdata]=cellfiletimes(cellid,rundata(1).batch);
   
   stimfiles={};
   respfiles={};
   for ii=1:length(cellfiledata),
      stimfiles{ii}=[cellfiledata(ii).stimpath,cellfiledata(ii).stimfile];
      respfiles{ii}=[cellfiledata(ii).path,cellfiledata(ii).respfile];
   end
   
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
   
   fidx=1;
   
   bstim=feval(stimloadcmd,stimfiles{fidx},1,0,stimloadparms{:});
   
   % filter stimulus segment if selected
   if ~isempty(stimfiltercmd),
      bstim=feval(stimfiltercmd,bstim,stimfilterparms{:});
   end
   % reshape to space X time if necessary
   iconside=size(bstim);
   iconside=iconside(1:(end-1));
   if length(iconside)>=2,
      % reshape all spatial dims into one
      bstim=reshape(bstim,prod(iconside),size(bstim,length(iconside)+1));
   end
   bstim=bstim'; % take transpose to put time in rows
   spacecount=size(bstim,2);
   
   disp(sprintf('Loading response: %s...',respfiles{fidx}));
   % resp=resploadatt(respfile,respvarname);
   % resp is time X space(phase) X attcode
   bresp=feval(resploadcmd,respfiles{fidx},resploadparms{:});   
   
   % filter response (eg resample, pick attentional state, etc)
   if ~isempty(respfiltercmd),
      bresp=feval(respfiltercmd,bresp,respfilterparms{:});
   end
   
   rlen=size(bresp,1);
   attcount=size(bresp,3);
   
   targ=zeros(rlen,1);
   for ii=2:attcount,
      targ(find(~isnan(bresp(:,1,ii))))=ii-1;
   end
   
   %respscale=respfilterparms{2}(1)-respfilterparms{1}(1);
   
   A=[targ bstim bresp(:,1,1).*1000];
   fprintf('output matrix A is %d X %d.\n',size(A,1),size(A,2));
   
   outfile=[outpath,cellid,'.',runclassstr{batchdata.runclassid+1},...
            '.',batchdata.kernfmt, ...
            '.',respfilefmtstr{batchdata.respfmtcode+1},'.dat'];
   fprintf('Saving to %s...\n',outfile);
   save(outfile,'-ascii','A');
end

return
fprintf('Loading %s...\n',[rundata(1).respath,rundata(1).kernfile,'.gz']);
z=zload([rundata(1).respath,rundata(1).kernfile,'.gz']);





