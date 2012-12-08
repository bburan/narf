% matchcell2file.m
%
% after saving spike file, record information in celldb
%   figure out cell identities and establish appropriate links
%   between gDataRaw --> gSingleCell via gSingleRaw and sCellFile
%
% created SVD 2005-10-01
%
dbopen;

sql=['SELECT * FROM gDataRaw WHERE parmfile like "%',fname,'%"'];
rawdata1=mysql(sql);

if ONEFILE,
   rawcount=1;
   rawdata2=[];
else
   rawcount=2;
   sql=['SELECT * FROM gDataRaw WHERE parmfile like "%',f2name,'%"'];
   rawdata2=mysql(sql);
end

if length(rawdata1)==1,
   
   % figure out how many cells were identified
   spikecount=zeros(size(spk));
   chanstr={'a','b','c','d','e'};
   cellids={};
   for ii=1:size(spk,1),
      for jj=1:size(spk,2),
         spikecount(ii,jj)=length(spk{ii,jj});
      end
      if extras.numChannels>1,
         cellids{ii}=sprintf('%s-%s%d',siteid,chanstr{str2num(chanNum)},ii);
      else
         cellids{ii}=sprintf('%s-%d',siteid,ii);
      end
   end
   unitcount=max(find(sum(spikecount,2)>0));
   
   % are there already more cells in db than measured units?
   sql=['SELECT id,cellid FROM gSingleCell WHERE masterid=',...
        num2str(rawdata1.masterid),...
        ' AND channum=',chanNum,' ORDER BY id'];
   singledata=mysql(sql);
   if length(singledata) > unitcount,
      disp('too many pre-existing units in db!');
      yn=input('delete extras ([y]/n)? ','s');
      if isempty(yn) | yn(1)~='n',
         for ii=unitcount+1:length(singledata),
            sql=['DELETE FROM sCellFile WHERE singleid=', ...
                 num2str(singledata(ii).id)];
            mysql(sql);
            sql=['DELETE FROM gSingleRaw WHERE singleid=', ...
                 num2str(singledata(ii).id)];
            mysql(sql);
            sql=['DELETE FROM gSingleCell WHERE id=',num2str(singledata(ii).id)];
            mysql(sql);
         end
      end
   end
   
   % cool, found a match. get info for relevant raw data files
   siteid=rawdata1.cellid;
   [rawdata,site,celldata,rcunit,rcchan]=dbgetsite(siteid);
   rawids=cat(1,rawdata.id);
   singleids=cat(1,celldata.id);
   if length(rawdata2)>0,
      r1=[find(rawdata1.id==rawids);
          find(rawdata2.id==rawids)];
   else
      r1=find(rawdata1.id==rawids);
   end
   if extras.numChannels>1,
      cc=nanmin(rcchan);
      cc=find(cc==str2num(chanNum));
   else
      cc=1:size(rcchan,2);
   end
   rawdata=rawdata(r1);
   rawids=rawids(r1);
   unitmaster=rcunit(1,:);
   rcunit=rcunit(r1,:);
   rcchan=rcchan(r1,:);
   
   % two possibilities: cells either exist or don't exist yet in db
   while length(cc)<unitcount,
      unitmaster=[unitmaster 0];
      rcunit=[rcunit zeros(length(r1),1)];
      rcchan=[rcchan zeros(length(r1),1)];
      singleids=[singleids;-1];
      cc=[cc length(singleids)];
      fprintf('adding cell %s\n',cellids{length(cc)});
   end
   
   rcchan(:,cc)=str2num(chanNum);
   
   for ii=1:length(cc),
      rcunit(:,cc(ii))=ii;
   end
   
   % save mapping between channels and cells
   dbsetsite(siteid,rcunit,rcchan,singleids,rawids);
   
   sql=['SELECT id,cellid FROM gSingleCell WHERE masterid=',...
        num2str(rawdata1.masterid),' ORDER BY cellid'];
   singledata=mysql(sql);
   singleids=cat(1,singledata.id);
   
   % save sorted cell info in sCellFile
   % respfilefmt="meska"; respfmtcode=2; addedby=sorter;
   % info=mfilename; cellid=cellids{}
   % spikes?
   % from gDataRaw etc: runclassid, masterid, singleid, singlerawid
   respfile={};resppath={};
   [respfile{1},resppath{1}] = basename([destin '.spk.mat']);
   [respfile{2},resppath{2}] = basename([destin2 '.spk.mat']);
   stimfile={};stimpath={};
   [stimfile{1},stimpath{1}] = basename(source);
   [stimfile{2},stimpath{2}] = basename(source2);
   if ~isempty(sorter),
      addedby=sorter;
   else
      addedby='david';
   end
   if rawcount==2,
      extras(2)=extras2;
   end
   
   for rr=1:rawcount
      
      delay = get(extras(rr).torcList.tag,'Onset');
      for ii=cc,
         % figure out singlerawid for this cell/raw combo
         sql=['SELECT * FROM gSingleRaw WHERE rawid=',num2str(rawids(rr)),...
              ' AND singleid=',num2str(singleids(ii))];
         singlerawdata=mysql(sql);
         if length(singlerawdata)>0,
            singlerawid=singlerawdata.id;
         else
            singlerawid=-1;
         end
         
         if spikecount(find(cc==ii),rr) ==0,
            sql=['DELETE FROM sCellFile',...
                 ' WHERE singleid=',num2str(singleids(ii)),...
                 ' AND rawid=',num2str(rawids(rr)),...
                 ' AND respfmtcode=0'];
            mysql(sql);
         else
            % check to see if sorted entry already exists for this cell/rawid
            sql=['SELECT * FROM sCellFile',...
                 ' WHERE singleid=',num2str(singleids(ii)),...
                 ' AND rawid=',num2str(rawids(rr)),...
                 ' AND respfmtcode=0'];
            cellfiledata=mysql(sql);
            
            if length(cellfiledata)==0,
               sqlinsert('sCellFile',...
                         'cellid',cellids{rcunit(rr,ii)},...
                         'masterid',rawdata(rr).masterid,...
                         'rawid',rawids(rr),...
                         'runclassid',rawdata(rr).runclassid,...
                         'path',resppath{rr},...
                         'resplen',extras(rr).npoint,...
                         'repcount',extras(rr).sweeps,...
                         'respfile',respfile{rr},...
                         'respvarname','sortinfo',...
                         'respfiletype',rcunit(rr,ii),...
                         'respfmtcode',0,...
                         'respfilefmt','meska',...
                         'stimfile',stimfile{rr},...
                         'stimfilecrf',delay,...
                         'stimfilefmt','par',...
                         'stimfmtcode',1,...
                         'addedby',addedby,...
                         'info',mfilename,...
                         'stimpath',stimpath{rr},...
                         'spikes',length(spksav1{rcunit(rr,ii)}),...
                         'singleid',singleids(ii),...
                         'channum',rcchan(rr,ii),...
                         'unit',rcunit(rr,ii),...
                         'singlerawid',singlerawid);
            else
               sql=['UPDATE sCellFile SET',...
                    ' masterid=',num2str(rawdata(rr).masterid),',',...
                    ' runclassid=',num2str(rawdata(rr).runclassid),',',...
                    ' path="',resppath{rr},'",',...
                    ' resplen=',num2str(extras(rr).npoint),',',...
                    ' repcount=',num2str(extras(rr).sweeps),',',...
                    ' respfile="',respfile{rr},'",',...
                    ' respvarname="sortinfo",',...
                    ' respfiletype=',num2str(rcunit(rr,ii)),',',...
                    ' respfmtcode=0,',...
                    ' respfilefmt="meska",',...
                    ' stimfile="',stimfile{rr},'",',...
                    ' stimfilecrf=',num2str(delay),',',...
                    ' stimfilefmt="par",',...
                    ' stimfmtcode=1,',...
                    ' addedby="',addedby,'",',...
                    ' info="',mfilename,'",',...
                    ' stimpath="',stimpath{rr},'",',...
                    ' spikes=',num2str(length(spksav1{rcunit(rr,ii)})),',',...
                    ' singleid=',num2str(singleids(ii)),',',...
                    ' singlerawid=',num2str(singlerawid),',',...
                    ' channum=',num2str(rcchan(rr,ii)),',',...
                    ' unit=',num2str(rcunit(rr,ii)),...
                    ' WHERE id=',num2str(cellfiledata.id)];
               mysql(sql);
            end
         end
      end
      
      sql=['UPDATE gDataRaw SET matlabfile="',resppath{rr},respfile{rr},'"',...
           ' WHERE id=',num2str(rawids(rr))];
      mysql(sql);
   end
   
elseif length(rawdata1)>1,
   disp('uh oh. multiple matches in celldb');
   keyboard
   
else
   disp('oh no! no matches found in celldb');
end


