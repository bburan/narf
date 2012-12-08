% function siteid=dbgetlastsite(Ferret, doingphysiology);
function siteid=dbgetlastsite(Ferret, doingphysiology);

if ~dbopen,
   siteid='';
   return;
end

sql=sprintf('SELECT * FROM gAnimal WHERE animal like "%s"',Ferret);
adata=mysql(sql);
if length(adata)==0,
    warning('%s not in celldb\n',Ferret);
    siteid='';
    return;
end
animal=adata(1).animal;

sql=['SELECT max(id) as maxid,to_days(now())-to_days(max(pendate)) as daylag',...
     ' FROM gPenetration' ...
     ' WHERE training=',num2str(1-doingphysiology),...
     ' AND animal="',animal,'"'];
lastpendata=mysql(sql);
daylag=lastpendata.daylag;

if isempty(lastpendata.maxid) | strcmp(lastpendata.maxid,'NULL'),
   
   warning(['No penetrations exist for this animal. Guessing SiteID' ...
            ' from scratch.']);
   if doingphysiology,
      siteid=[adata(1).cellprefix,'001a'];
   else
      siteid=[adata(1).cellprefix,'001Ta'];
   end
   return
   
else
   
   sql=['SELECT penid,max(siteid) as siteid FROM gCellMaster'...
       ' WHERE penid=',num2str(lastpendata.maxid),' GROUP BY penid'];
   tpendata=mysql(sql);
   
   if length(tpendata)==0,
       siteid='';
   else
       siteid=tpendata.siteid
   end
   if isempty(siteid),
       % no sites for this penetration yet
       sql=['SELECT * FROM gPenetration'...
           ' WHERE id=',num2str(lastpendata.maxid)];
       lastpendata=mysql(sql);
       siteid=[lastpendata.penname,'a'];
   end


   if daylag>=1,
       numidx=find(siteid>='0' & siteid<='9');
       pennum=str2num(siteid(numidx));
        
       siteid=[siteid(1:numidx(1)-1) sprintf('%03d',pennum+1) ...
           siteid(numidx(end)+1:end)];
       fprintf('guessing new penetration/site: %s\n',siteid);
           
   end
   
end
