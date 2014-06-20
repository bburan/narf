function dbcopybatch(batchin)
    
dbopen;

%batchin=251;
celllist=[];

% list for dgrat/dnat
%celllist=['("R148A","R156A","R158A","R162B","R164C","R166C","R168B",',...
%          '"R170A","R206B","R208D","R210A","R211A","R212B","R215B",',...
%          '"R217B","R219B","R220A","R221A","R223A","R225C")'];

% list for nat (long reviews)
%%celllist=['("R110A","R112A","R143C","R206B","R208D","R210A","R211A",',...
%          '"R212B","R213C","R214A","R215A","R215B","R217B","R219B",',...
%          '"R220A","R221A","R222A","R223A","R225A","R225C")'];

sql=['SELECT * from sBatch WHERE id=',num2str(batchin)];
batchdata=mysql(sql);

if isempty(celllist),
   sql=['SELECT * FROM sRunData WHERE batch=',num2str(batchin)];
else
   sql=['SELECT * FROM sRunData WHERE batch=',num2str(batchin),...
        ' AND cellid in ',celllist];
end

rundata=mysql(sql);

fprintf('Copying %d runs for batch %d (%s)...\n',length(rundata), ...
        batchin,batchdata.name);

batchname=input('Enter new name ['''']: ','s');

%rcsetstrings;
sql=['SELECT * FROM gRunClass WHERE id=',...
     num2str(batchdata.runclassid)];
rcdata=mysql(sql);
sbase=['.',rcdata.name,'.',batchdata.kernfmt];
ss=input(sprintf('resfile identifier [%s]: ',sbase),'s');
if ~isempty(ss),
   sbase=ss;
end

resbase=[sbase,'.res.mat'];
kernbase=[sbase,'.kern.mat'];


f=fieldnames(batchdata);
snames='INSERT INTO sBatch (';
svalues=') VALUES (';
bfirst=1;
for ii=1:length(f)
   if not(strcmp(f{ii},'id')),
      t=batchdata.(f{ii});
      if strcmp(f{ii},'name'),
         sv=['"',batchname,'"'];
      elseif strcmp(f{ii},'parmstring'),
         sv=['"',strrep(char(t),'"','\"'),'"'];
      elseif isempty(t),
            sv='null';
      elseif isnumeric(t),
         sv=num2str(t);
      else
         sv=['"',t,'"'];
      end
      
      if bfirst,
         snames=[snames,f{ii}];
         svalues=[svalues,sv];
         bfirst=0;
      else
         snames=[snames,',',f{ii}];
         svalues=[svalues,',',sv];
      end
   end
end

sql=[snames,svalues,');']
[res,aff,batchidx]=mysql(sql);

fprintf(['added batch %d. Use dbbatchfill to add cells and dbbatchcells ' ...
         'to set selection criteria\n'],batchidx);

return

for runidx=1:length(rundata),
   f=fieldnames(rundata);
   snames='INSERT INTO sRunData (';
   svalues=') VALUES (';
   bfirst=1;
   
   for ii=1:length(f)
      if not(strcmp(f{ii},'id')),
         t=getfield(rundata(runidx),f{ii});
         if strcmp(f{ii},'resfile'),
            sv=['"',rundata(runidx).cellid,resbase,'"'];
         elseif strcmp(f{ii},'kernfile'),
            sv=['"',rundata(runidx).cellid,kernbase,'"'];
         elseif strcmp(f{ii},'batch'),
            sv=num2str(batchidx);
         elseif strcmp(f{ii},'respath'),
            tpos=findstr(t,['/batch',num2str(rundata(runidx).batch)]);
            if isempty(tpos),
               sv=['"',t,'"'];
            else
               sv=['"',t(1:tpos),'batch',num2str(batchidx),'/"'];
               if ~exist([t(1:tpos),'batch',num2str(batchidx)],'dir'),
                  fprintf('creating dir %s\n',...
                          [t(1:tpos),'batch',num2str(batchidx)]);
                  unix(['mkdir ',t(1:tpos),'batch',num2str(batchidx)]);
               end
            end
         elseif isempty(t),
            sv='null';
         elseif isnumeric(t),
            sv=num2str(t);
         else
            sv=['"',t,'"'];
         end
         
         if bfirst,
            snames=[snames,f{ii}];
            svalues=[svalues,sv];
            bfirst=0;
         else
            snames=[snames,',',f{ii}];
            svalues=[svalues,',',sv];
         end
      end
   end
   
   sql=[snames,svalues,');']
   [res,aff,runidx]=mysql(sql);
end

