dbopen(1);

sql=['SELECT * FROM sCellFile where respfmtcode=0 and stimfmtcode=0'];
cellfiledata=mysql(sql);

for ii=1:length(cellfiledata),
   respfile=[cellfiledata(ii).path,cellfiledata(ii).respfile];
   stimfile=[cellfiledata(ii).stimpath,cellfiledata(ii).stimfile];
   
   fprintf('%s: ',cellfiledata(ii).respfile);
   
   if exist(respfile,'file'),
      
      stimfilecrf=cellfiledata(ii).stimfilecrf;
      spikecount=cellfiledata(ii).spikes;
      nreps=cellfiledata(ii).repcount;
      resplen=cellfiledata(ii).resplen;
      respfiletype=cellfiledata(ii).respfiletype;
      
      if ~(stimfilecrf>0 & resplen>0 & nreps>0 & spikecount>0),
         
         if strcmp(cellfiledata(ii).respfile(end-3:end),'.mat'),
            respfiletype=0;
            r=load(respfile);
            if isfield(r,'r'),
               r=[r.psth,r.r];
            else
               r=getfield(r,cellfiledata(ii).respvarname);
            end
         else
            respfiletype=1;
            r=respload(respfile,cellfiledata(ii).respvarname);
         end
         
         resplen=sum(~isnan(r(:,1)),1);
         nreps=size(r,2);
         r(find(r==-1))=nan;
         if nreps > 1,
            nshow=sum(isnan(r(:,2:nreps)),2);
            nshow=nshow(find(nshow<nreps-1));
            nshow=max(nshow);
            if isempty(nshow),
               nshow=1;
            end
            nreps=nreps-nshow-1;
            if nreps==0,  % maybe one frame always got dropped?
               nreps=1;
            end
         end
         
         r(isnan(r))=0;
         if nreps==0,
            nreps=0;
            spikecount=0;
         elseif size(r,2)==1,
            nreps=1;
            spikecount=sum(r);
         else
            %spikecount=sum(sum(r(:,2:size(r,2))));
            spikecount=nansum(r(:,1));
         end
         spikecount=round(spikecount);
      end
      
      [framecount,stimwindowsize,ss,pp,filetype]=imfileinfo(stimfile);
      
      sampleframe=loadimfile(stimfile,1,1);
      s=size(sampleframe);
      stimiconside='';
      for ss=1:length(s);
         stimiconside=[stimiconside,',',num2str(s(ss))];
      end
      stimiconside=stimiconside(2:end);
      if cellfiledata(ii).stimfmtcode==0, % ie, pixel
         stimwindowsize=s(1);
      else
         % need to look at raw file for size
         keyboard
      end
      
      sql=['SELECT * FROM gCellMaster WHERE id=',...
           num2str(cellfiledata(ii).masterid)];
      celldata=mysql(sql);
      
      if length(celldata)>0 & ~isempty(celldata.rfsize) & celldata.rfsize>0,
         stimfilecrf=stimwindowsize./celldata.rfsize;
      else
         stimfilecrf=cellfiledata(ii).stimfilecrf;
      end
      
      if cellfiledata(ii).runclassid==1 & ...
             strcmp(cellfiledata(ii).info(1:6),'dbcrap'),
         stimwindowsize=cellfiledata(ii).stimwindowsize;
         stimfilecrf=cellfiledata(ii).stimfilecrf;
      end
      
      if cellfiledata(ii).respfiletype~=respfiletype | ...
             cellfiledata(ii).resplen~=resplen | ...
             cellfiledata(ii).spikes~=spikecount| ...
             cellfiledata(ii).repcount~=nreps| ...
             cellfiledata(ii).stimwindowsize~=stimwindowsize| ...
             round(cellfiledata(ii).stimfilecrf*10)~=round(stimfilecrf*10),
             
         
         fprintf('field: OLD --> NEW\n');
         fprintf('type:  %d  --> %d\n',cellfiledata(ii).respfiletype,...
                 respfiletype);
         fprintf('resplen:  %d  --> %d\n',cellfiledata(ii).resplen,resplen);
         fprintf('spikes:  %d  --> %d\n',cellfiledata(ii).spikes,spikecount);
         fprintf('reps:  %d  --> %d\n',cellfiledata(ii).repcount,nreps);
         fprintf('winsize:  %d  --> %d\n',cellfiledata(ii).stimwindowsize,...
                 stimwindowsize);
         fprintf('icon:  %s  --> %s\n',cellfiledata(ii).stimiconside,...
                 stimiconside);
         fprintf('crfs:  %.1f --> %.1f\n',cellfiledata(ii).stimfilecrf,...
                 stimfilecrf);
         
         yn=input('change? [n]','s'),
      
         if length(yn)>0 & yn(1)=='y',
            sql=['UPDATE sCellFile SET',...
                 ' spikes=',num2str(spikecount),',',...
                 ' respfiletype=',num2str(respfiletype),',',...
                 ' repcount=',num2str(nreps),',',...
                 ' resplen=',num2str(resplen),',',...
                 ' stimwindowsize=',num2str(stimwindowsize),',',...
                 ' stimfilecrf=',num2str(stimfilecrf),',',...
                 ' stimiconside="',stimiconside,'"',...
                 ' WHERE id=',num2str(cellfiledata(ii).id)];
            %keyboard
            mysql(sql);
         end
      end
      
   else
      fprintf('not found\n');
   end
end


