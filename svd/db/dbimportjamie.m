%  read in all of jamie's raw data

dbopen;
rootpath='/auto/data/archive/shadowfax/Monks/2001/'

filelist=jls([rootpath,'*/*.*.[0-9][0-9][0-9].gz']);
lastcellid='';
for fileidx=1:length(filelist),
   [respfile,resppath]=basename(filelist{fileidx});
   % remove .gz from tail
   respfile=respfile(1:end-3);
   filedetails=strsep(respfile,'.',1);
   cellid=filedetails{1};
   task=filedetails{2};
   tasknum=filedetails{3};
   
   pathdetails=strsep(resppath(1:end-1),'/',1);
   penname=pathdetails{end};
   
   % add/update information in gDataRaw
   skipfile=0;
   if ~sum([findstr(cellid,'m0') findstr(cellid,'z0')]) | ...
         strcmp(cellid,'m0000') | strcmp(cellid,'z0000'),
      skipfile=1;
   elseif findstr(task,'train'),
      runclass='FS';
      runclassid=3;
   elseif findstr(task,'zvs'),
      runclass='FS';
      runclassid=3;
   elseif findstr(task,'imview'),
      runclass='IMV';
      runclassid=9;
   elseif findstr(task,'eyecal'),
      runclass='EYEC';
      runclassid=8;
   else
      skipfile=1;
   end
   fprintf('%s: %s\n',task,filelist{fileidx});
   if ~skipfile
      keyboard;
   end
   
   if ~skipfile & ~strcmp(lastcellid,cellid)
   
      mdbinfo=frinfo(cellid);
      
      if isempty(mdbinfo),
         if strcmp(cellid(1:2),'mo'),
            animal='model';
         elseif ismember(cellid(1),['m','M']),
            animal='mac';
         elseif ismember(cellid(1),['z','Z']),
            animal='ziggy';
         elseif ismember(cellid(1),['r','R']),
            animal='reese';
         end
         
         % need to guess well
         sql=['SELECT well,cellid FROM gCellMaster',...
              ' WHERE penname like "',penname(1:end-2),'%" ORDER BY cellid'];
         welldata=mysql(sql);
         if length(welldata)==0,
            sql=['SELECT well,cellid FROM gCellMaster',...
                 ' WHERE cellid<"',cellid,'" ORDER BY cellid'];
            welldata=mysql(sql);
         end
         if length(welldata)>0,
            well=welldata(end).well;
            fprintf('guessing well %d\n',well);
         else
            sql=['SELECT max(well) as well FROM gCellMaster',...
                 ' WHERE animal="',animal,'"'];
            welldata=mysql(sql);
            well=welldata.well;
         end
         
         xoffset=0;
         yoffset=0;
         rfsize=0;
         eyecal='';
         ppd=0;
         area='?';
         owner='jm';
      else
         animal=mdbinfo.id{3};
         well=mdbinfo.id{4};
         
         xoffset=mdbinfo.rf{1};
         yoffset=mdbinfo.rf{2};
         rfsize=mdbinfo.rf{3}*2;
         if isnan(rfsize),
            xoffset=mdbinfo.fid{1};
            yoffset=mdbinfo.fid{2};
            rfsize=mdbinfo.fid{3}*2;
         end
         if isnan(rfsize),
            xoffset=0;
            yoffset=0;
            rfsize=0;
         end
         eyecal=mdbinfo.smap3{1};
         ppd=mdbinfo.ppd;
         if isfield(mdbinfo,'area'),
            area=mdbinfo.area;
         else
            area='v?';
         end
         if isfield(mdbinfo,'owner'),
            owner=mdbinfo.owner;
         else
            owner='jm?';
         end
      end
      
      sql=['SELECT * FROM gCellMaster where cellid="',cellid,'"'];
      celldata=mysql(sql);
      sql=['SELECT * FROM gPenetration where penname="',penname,'"'];
      matchedpendata=mysql(sql);
      
      if length(celldata)>1,
         disp('error: two cell entries!');
         keyboard
         
      elseif length(celldata)>0,
         % cell already exists. set masterid and penid appropriately
         % for updates
         
         masterid=celldata.id;
         
         sql=['SELECT * FROM gPenetration where id=',num2str(celldata.penid)];
         pendata=mysql(sql);
         
         % check to see if penetration is correct.
         if length(pendata)==0,
            disp('error: no pen for existing cell!');
            keyboard;
         elseif length(matchedpendata)==0 & ~strcmp(pendata.penname,penname),
            fprintf('cell %s, setting penname to %s\n',cellid,penname);
            penid=pendata.id;
            pendate=pendata.pendate;
         elseif ~strcmp(pendata.penname,penname),
            fprintf('cell %s, switching pen to existing penname %s\n',...
                    cellid,penname);
            penid=matchedpendata.id;  % swtiched to new penid!
            pendate=matchedpendata.pendate;
         else
            % cellmaster and penid are ok
            fprintf('cell %s matches correct penname: %s\n',cellid,penname);
            penid=pendata.id;
            pendate=pendata.pendate;
         end
      else
         % cell doesn't exist, add it
         
         pendate='';
         if length(matchedpendata)==0,
            % pen doesn't exist, create it
            fprintf('creating pentration penname=%s\n',penname);
            
            [aff,penid]=sqlinsert('gPenetration',...
                                  'penname',penname,...
                                  'well',well,...
                                  'who',owner,...
                                  'animal',animal,...
                                  'addedby','david',...
                                  'info','dbscratch.m');
         else
            penid=matchedpendata.id;
         end
         
         [aff,masterid]=sqlinsert('gCellMaster',...
                                  'cellid',cellid,...
                                  'animal',animal,...
                                  'well',well,...
                                  'penname',penname,...
                                  'penid',penid,...
                                  'addedby','david',...
                                  'info','dbscratch.m');
         
      end
      
      % update gCellMaster and gPenetration to reflect info in jamie's
      % mdb (via fileinfo.m);
      if length(pendate)==0,
         pendate=penname;
      end
      sql=['UPDATE gPenetration SET ',...
           'animal="',animal,'",',...
           'well=',num2str(well),',',...
           'etudeg=',num2str(ppd),',',...
           'penname="',penname,'",',...
           'pendate="',pendate,'",',...
           'who="',owner,'",',...
           'info="dbscratch.m"',...
           ' WHERE id=',num2str(penid)];
      mysql(sql);
      
      sql=['UPDATE gCellMaster SET ',...
           'animal="',animal,'",',...
           'well=',num2str(well),',',...
           'ppd=',num2str(ppd),',',...
           'rfsize=',num2str(rfsize),',',...
           'xoffset=',num2str(xoffset),',',...
           'yoffset=',num2str(yoffset),',',...
           'eyecal="',eyecal,'",',...
           'penname="',penname,'",',...
           'penid=',num2str(penid),',',...
           'area="',area,'",',...
           'info="dbscratch.m"',...
           ' WHERE id=',num2str(masterid)];
      mysql(sql);
      
      lastcellid=cellid;
   end
   
   if ~skipfile,
      sql=['SELECT * FROM gDataRaw',...
           ' WHERE masterid=',num2str(masterid),...
           ' AND respfile like "%',respfile,'%"'];
      rawdata=mysql(sql);
      if length(rawdata)==0,
         sql=['SELECT * FROM gDataRaw',...
              ' WHERE masterid=',num2str(masterid),...
              ' AND matlabfile like "%',respfile,'%"'];
         rawdata=mysql(sql);
      end
      
      if length(rawdata)==0,
         [aff,rawid]=sqlinsert('gDataRaw','masterid',masterid);
      elseif length(rawdata)>1,
         disp('more than one gRawData entry???');
         keyboard;
      else
         rawid=rawdata.id;
      end
      
      sql=['UPDATE gDataRaw SET ',...
           'cellid="',cellid,'",',...
           'runclassid=',num2str(runclassid),',',...
           'runclass="',runclass,'",',...
           'resppath="',resppath,'",',...
           'respfile="',respfile,'",',...
           'task="',task,'",',...
           'addedby="david",',...
           'info="dbscratch.m"',...
           ' WHERE id=',num2str(rawid)];
      mysql(sql);
      
%          '=',num2str(),',',...
%          '="',,'",',...
%| stimfile    | varchar(255)     | YES  |     | NULL    |                |
%| matlabfile  | varchar(255)     | YES  |     | NULL    |                |
%| stimpath    | varchar(255)     | YES  |     | NULL    |                |
%| stimspeedid | int(11)          | YES  |     | NULL    |                |
         
   end
   
end
