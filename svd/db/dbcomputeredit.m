% function dbcomputeredit(name,ext,owner[=''],allowothers[=1],maxproc[=1]);
function dbcomputeredit(name,ext,owner,allowothers,maxproc);

if strcmp(ext,'bic'),
   location=1;
elseif strcmp(ext,'Millennium'),
   location=2;
elseif strcmp(ext,'fet'),
   location=3;
else % ext is jlg
   location=0;
end

allowqueuemaster=1;

% check to see if computer is already in db
dbopen;
sql=['SELECT * FROM tComputer WHERE name="',name,'" and ext="',ext,'"'];
compdata=mysql(sql);

if length(compdata)>0,
   
   % yes, it's in the db, update it.
   
   fprintf('computer %s.%s exists. updating db.\n',name,ext);
   
   if ~exist('owner','var') | isempty(owner),
      owner=compdata.owner;
   end
   if ~exist('allowothers','var') | isempty(allowothers),
      allowothers=compdata.allowothers;
   end
   if ~exist('maxproc','var') | isempty(maxproc),
      maxproc=compdata.maxproc;
   end
   
   sql=['UPDATE tComputer',...
        ' SET owner="',owner,'",',...
        ' allowothers=',num2str(allowothers),',',...
        ' maxproc=',num2str(maxproc),...
        ' WHERE id=',num2str(compdata.id)]
   mysql(sql);
   
else
   
   % no, it's not in db, add it.
   
   fprintf('computer %s.%s does not exist. adding to db.\n',name,ext);
   
   if ~exist('owner','var'),
      owner='';
   end
   
   if ~exist('allowothers','var'),
      if isempty(owner),
         allowothers=1;
      else
         allowothers=0;
      end
   end
   
   if ~exist('maxproc','var'),
      maxproc=1;
   end
   
   sqlinsert('tComputer',...
             'name',name,...
             'location',location,...
             'maxproc',maxproc,...
             'allowqueuemaster',allowqueuemaster,...
             'ext',ext,...
             'owner',owner,...
             'allowothers',allowothers);
end

