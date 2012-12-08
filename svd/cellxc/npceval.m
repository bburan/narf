% npceval
%
% 1. connect to npc database
% 2. query for next new predset entry
% 3. check if user hasn't overrun quota
% 4. if not, load predfile and val resp file for each cell
% 5. eval xc, save to predfile for each cell
% 6. other pred scores?
% 7. save avg stats to predset.
% 8. loop to 2
%
%
%


dbopen('localhost','david','nine1997','npc');
mysql('use npc');


% figure out the name of this machine:
[ww,hh]=unix('hostname');

while 1,
   % loop forever
   sql='SELECT * FROM nPredSet WHERE isnull(evaltime)';
   preddata=mysql(sql);
   
   
   for ii=1:length(preddata),
      errorflag=0;
      predsetid=preddata(ii).id;
      dsid=preddata(ii).datasetid;
      infile=[preddata(ii).path preddata(ii).predfile];
      [pp,bb,ee]=fileparts(infile);
      logfile=fullfile(pp,[bb,'.log'])
      
      fid=fopen(logfile,'w');
      npclog(fid,'Opening file %s (user %s, algorithm %s v %.2f)...',...
              preddata(ii).predfile,preddata(ii).userid,...
              preddata(ii).algorithm,preddata(ii).version);
      
      
      try
         pp=load([preddata(ii).path preddata(ii).predfile]);
         npclog(fid,'Loaded successfully.');
      catch 
         pp=[];
         errorflag=1;
         errormsg='ERROR 1: Could not open file.';
      end
      
      if ~errorflag,
         % ok, loaded the file. 
         % now check to make sure that all the variables are defined
         if ~isfield(pp,'prediction'),
            errorflag=2;
            errormsg='ERROR 2: Missing structure array ''prediction''.';
         elseif ~isfield(pp.prediction,'cellid'),
            errorflag=3;
            errormsg='ERROR 3: Missing field ''cellid''.';
         elseif ~isfield(pp.prediction,'response'),
            errorflag=4;
            errormsg='ERROR 4: Missing field ''response''.';
         end
      end
      
      if ~errorflag,
         sql=['SELECT * FROM nDataSet WHERE id=',num2str(dsid)];
         setdata=mysql(sql);
         
         npclog(fid,'Evaluating each prediction vector...\n');
         
         predxc=zeros(length(pp.prediction),1);
         for jj=1:length(pp.prediction),
            sql=['SELECT * FROM nDataFile WHERE datasetid=',num2str(dsid),...
                 ' AND cellid="',pp.prediction(jj).cellid,'"'];
            filedata=mysql(sql);
            
            rr=load([filedata.path filedata.valrespfile]);
            
            rridx=find(~isnan(rr.resp));
            
            tpred=pp.prediction(jj).response;
            if length(tpred)<max(rridx),
               errorflag=5;
               terrorflag=5;
               errormsg=['ERROR 5: Prediction vector size error for' ...
                         ' some or all neurons.'];
               results='Prediction vector size error';
               predxc(jj)=0;
               npclog(fid,...
                'Cell %s. ERROR 5: Expected %d x 1 prediction vector. cc=0',...
                      pp.prediction(jj).cellid,length(tpred));
            else
               terrorflag=0;
               results='Prediction evaluated successfully';
               if std(pp.prediction(jj).response(rridx))>0 & ...
                  std(rr.resp(rridx))>0,
std(pp.prediction(jj).response(rridx))
std(rr.resp(rridx))
                  predxc(jj)=xcov(pp.prediction(jj).response(rridx),...
                                  rr.resp(rridx),0,'coeff');
               end
               npclog(fid,'Cell %s (data %d x 1): cc=%.3f',...
                      pp.prediction(jj).cellid,length(tpred),...
                      predxc(jj));
            end
            
            % delete existing prediction entry (shouldn't be one?)
            sql=['DELETE FROM nPred WHERE predsetid=',num2str(preddata(ii).id),...
                 ' AND datafileid=',num2str(filedata.id)];
            mysql(sql);
            
            % record outcome for this prediction set for this cell
            sqlinsert('nPred',...
                      'datasetid',dsid,...
                      'datafileid',filedata.id,...
                      'predsetid',predsetid,...
                      'userid',preddata(ii).userid,...
                      'errorflag',terrorflag,...
                      'results',results,...
                      'predxc',round(predxc(jj)*1000)/1000,...
                      'addedby',preddata(ii).addedby,...
                      'addeddate',datestr(now,31),...
                      'info','npceval.m');
            
         end
         
      end
      
      % record outcome to prediction set placeholder
      if errorflag
         sql=['UPDATE nPredSet set evaltime=now(),'...
              ' logfile="', basename(logfile),'",',...
              ' errorflag=',num2str(errorflag),',',...
              ' results="',errormsg,'"',...
              ' WHERE id=',num2str(predsetid)];
         mysql(sql);
         
         npclog(fid,errormsg);
         npclog(fid,'Dataset validation NOT completed successfully.');
      else
         sql=['UPDATE nPredSet set evaltime=now(),'...
              ' logfile="', basename(logfile),'",',...
              ' errorflag=0,results="Prediction evaluated successfully",',...
              ' predxc=',sprintf('%.3f',mean(predxc)),...
              ' WHERE id=',num2str(predsetid)];
         mysql(sql);
         npclog(fid,'Dataset validation completed successfully.');
      end
   end
   
   % tell the db that daemon is still alive
   sql=['UPDATE nDaemon set lasttick=now(),host="' hh '"'];
   mysql(sql);
   
   % wait before infinite loop
   pause(10);
   
end
