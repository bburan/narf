
dbopen;
fid=fopen('/home/svd/nsl/users/svd/Weights-all ferrets.txt','r');
if fid==-1,
   error('%s: Error. File not found.',mfilename);
end
buffer=char(fread(fid,inf,'char'))';
fclose(fid);

eolidx=find(buffer==10);
lasteolidx=0;

bdata=strsep(buffer,char(10),1);
for ii=1:length(bdata),
   
   tdata=strsep(bdata{ii},sprintf('\t'),1);
   if ii==1,
      data=cell(length(eolidx),length(tdata));
   end
   for jj=1:length(tdata),
      data{ii,jj}=tdata{jj};
   end
end

dates={data{:,1}};
animals={data{531,:}};
animalidx=[];
for ii=1:33,
   animals{ii}=strrep(animals{ii},char(0),'');
   if ~isempty(animals{ii}) & ~isempty(strtrim(animals{ii})),
      animalidx=[animalidx;ii];
   end
end

for ii=1320:length(dates),
   if isempty(dates{ii}) | findstr('%',dates{ii}),
      dn=[];
   else
      try,
         dn=datenum(dates{ii});
      catch
         dn=[];
      end
   end
   if ~isempty(dn),
      dd=datestr(dn,29);
      
      fprintf('%s\t',dd);
      for jj=1:length(animalidx),
         wstr=data{ii,animalidx(jj)};
         if isempty(wstr),
            ww=[];
            fprintf(' \t');
         else
            ww=str2num(data{ii,animalidx(jj)});
            
         end
         if ~isempty(ww),
            fprintf('%.0f\t',ww);
            astr=strtrim(lower(animals{animalidx(jj)}));
            sql=['SELECT * FROM gPenetration'...
                 ' WHERE animal like "',astr,'"',...
                 ' AND pendate="',dd,'"'];
            wdata=mysql(sql);
            
            if length(wdata)>0,
               if wdata(1).training<2,
                  %keyboard
               end
               penid=wdata(1).id;
               sql=['UPDATE gPenetration SET',...
                    ' weight=',num2str(ww),...
                    ' WHERE id=',num2str(penid)];
               mysql(sql);
               %fprintf('updated gPenetration entry %d\n',penid);
            else
               %keyboard
               [aff,penid]=sqlinsert('gPenetration',...
                                     'penname',dd,...
                             'animal',astr,...
                             'pendate',dd,...
                             'weight',ww,...
                             'training',2,...
                             'addedby','david',...
                             'info','dbscanweights.m');
               %fprintf('added gPenetration entry %d\n',penid);
               
            end
         end
      end
      fprintf('\n');
      
   end
end
