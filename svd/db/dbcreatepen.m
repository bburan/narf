% function penid=dbcreatepen(params)
%
% created SVD 2005-11-21
%
function penid=dbcreatepen(params)

dbopen;

sql=sprintf('SELECT * FROM gAnimal WHERE animal like "%s"',...
            params.Ferret);
adata=mysql(sql);
if length(adata)==0,
   error('Ferret not found in cellDB!');
end
params.animal=adata(1).animal;

sql=sprintf('SELECT * FROM gUserPrefs WHERE realname like "%%%s%%"',...
            params.Tester);
udata=mysql(sql);
if length(udata)==0,
   error('Tester not found in cellDB!');
end
params.who=udata(1).userid;


sql=['SELECT max(id) as maxid FROM gPenetration' ...
     ' WHERE training in (0,1) AND animal="',params.animal,'"'];
lastpendata=mysql(sql);

if isempty(lastpendata.maxid),
   
   warning('No penetrations exist for this animal. Guessing info from scratch.');
   
   if ~isfield(params,'well'),
      params.well=1;
   end
   params.cellprefix=adata(1).cellprefix;
   params.ear='b';
   params.eye='';
   params.mondist=0;
   params.etudeg=0;
   params.racknotes='';
   params.speakernotes='';
   params.probenotes='';
   params.electrodenotes='';
else
   
   sql=['SELECT gPenetration.*, gAnimal.cellprefix FROM gPenetration'...
        ' INNER JOIN gAnimal ON gAnimal.animal=gPenetration.animal'...
        ' WHERE gPenetration.id=',num2str(lastpendata.maxid)];
   lastpendata=mysql(sql);
   
   
   fprintf('Guessing info from pen %s\n', lastpendata.penname);

   if ~isfield(params,'well'),
      params.well=lastpendata.well;
   end
   params.cellprefix=lastpendata.cellprefix;
   params.ear=lastpendata.ear;
   params.eye=lastpendata.eye;
   params.mondist=lastpendata.mondist * 1.0;
   params.etudeg=lastpendata.etudeg * 1.0;
   params.racknotes=lastpendata.racknotes;
   params.speakernotes=lastpendata.speakernotes;
   params.probenotes=lastpendata.probenotes;
   params.electrodenotes=lastpendata.electrodenotes;
end

params.numchans=params.NumberOfElectrodes;
if strcmp(params.Physiology,'No'),
    params.training=1;
else
    params.training=0;
end


sql=['SELECT * FROM gPenetration'...
     ' WHERE animal="',params.animal,'"',...
     ' AND pendate="',params.date,'" AND training=2'];
wdata=mysql(sql);

if length(wdata)>0,
   penid=wdata(1).id;
   sql=['UPDATE gPenetration SET',...
        ' penname="',params.penname,'",',...
        'animal="',params.animal,'",',...
        'well=',num2str(params.well),',',...
        'pendate="',params.date,'",',...
        'who="',params.who,'",',...
        'fixtime="',datestr(now,'HH:MM'),'",',...
        'ear="',params.ear,'",',...
        'numchans=',num2str(params.numchans),',',...
        'racknotes="',params.racknotes,'",',...
        'speakernotes="',params.speakernotes,'",',...
        'probenotes="',params.probenotes,'",',...
        'electrodenotes="',params.electrodenotes,'",',...
        'training=',num2str(params.training),',',...
        'addedby="',params.who,'",',...
        'info="dbcreatepen.m"',...
        ' WHERE id=',num2str(penid)];
   mysql(sql);
   fprintf('updated gPenetration entry %d\n',penid);
else
   %keyboard
   [aff,penid]=sqlinsert('gPenetration',...
                         'penname',params.penname,...
                         'animal',params.animal,...
                         'well',params.well,...
                         'pendate',params.date,...
                         'who',params.who,...
                         'fixtime',datestr(now,'HH:MM'),...
                         'ear',params.ear,...
                         'numchans',params.numchans,...
                         'racknotes',char(params.racknotes),...
                         'speakernotes',char(params.speakernotes),...
                         'probenotes',char(params.probenotes),...
                         'electrodenotes',char(params.electrodenotes),...
                         'training',params.training,...
                         'addedby',params.who,...
                         'info','dbcreatepen.m');
   fprintf('added gPenetration entry %d\n',penid);
end
