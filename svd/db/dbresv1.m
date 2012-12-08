% function [predxc,cellid,prederr]=dbresv1(batchid,predmetric[=0],verbose[=1])
%
% predxc - valset X nlidx X respdim X cellid matrix
% cellid - cell array with cellids corresponding to the fourth
%          dimension of predxc
% prederr - valset X nlidx X respdim X cellid matrix with 95%
%           confidence intervals for each predxc entry
%
function [predxc,cellid,prederr]=dbresv1(batchid,predmetric,verbose)

dbopen;

if ~exist('batchid','var'),
   batchid=26
end
if ~exist('predmetric','var'),
   predmetric=0;
end
if ~exist('verbose','var'),
   verbose=1;
end

% concatenate vector of batch ids into a string for SQL "in" command
batstr='(';
for ii=1:length(batchid),
   batstr=[batstr,num2str(batchid(ii)),','];
end
batstr(end)=')';

sql=['SELECT sResults.*,masterid,cellid',...
     ' FROM sResults INNER JOIN sRunData ON sResults.runid=sRunData.id',...
     ' WHERE sResults.batch in ',batstr,...
     ' AND (sRunData.batch in (127,128) OR not(cellid like "%model%"))',...
     ' ORDER BY cellid'];
resdata=mysql(sql);

ttt=cat(1,resdata.runid);

for ii=1:length(resdata),
   eval(char(resdata(ii).matstr));
end

cellid={preddata(ttt).cellid};
cellcount=length(cellid);

% preddata(ttt(1)).predxc: batchcount x strfcount x attcount x latcount
s=[size(preddata(ttt(1)).predxc) 1 1 1];

if predmetric==0,
   predxc=cat(5,preddata(ttt).predxc);
   prederr=cat(5,preddata(ttt).prederr);
elseif predmetric==1,
   predxc=cat(5,preddata(ttt).predfix);
   prederr=cat(5,preddata(ttt).predfixerr);
elseif predmetric==2,
   predxc=cat(5,preddata(ttt).predinf);
   prederr=cat(5,preddata(ttt).prederr);
end

%predxc=cat(5,preddata(ttt).predfix);
predxc=reshape(predxc(:,:,1,:,:),s(1),s(2),s(4),cellcount);
prederr=reshape(prederr(:,:,1,:,:),s(1),s(2),s(4),cellcount);
if ismember(batchid,[24 26 37 51 80 83 86]);
   badcellidx=find(strcmp(cellid,'r0150B') | strcmp(cellid,'93G83A'));
   predxc(:,:,:,badcellidx)=nan;
   prederr(:,:,:,badcellidx)=nan;
end
predxc(isinf(predxc))=nan;


predxc=cat(4,predxc,...
           permute(nanmean(permute(predxc,[4 1 2 3])),[2 3 4 1]),...
           permute(nanmean(permute(predxc.*abs(predxc),[4 1 2 3])),[2 3 4 1]));
prederr=cat(4,prederr,...
           permute(nanmean(permute(prederr,[4 1 2 3])),[2 3 4 1]),...
           permute(nanmean(permute(prederr.*abs(prederr),[4 1 2 3])),[2 3 4 1]));

if prod(size(preddata(ttt(1)).predp))>0,
   predp=cat(5,preddata(ttt).predp);
   predp=reshape(predp(:,:,1,:,:),s(1),s(2),s(4),cellcount);
   predp=cat(4,predp,...
             permute(nanmean(permute(predp,[4 1 2 3])),[2 3 4 1]),...
             permute(nanmean(permute(predp,[4 1 2 3])),[2 3 4 1]));
else
   predp=zeros(size(predxc));
end

cellcount=cellcount+1;
cellid{cellcount}='mean';
cellcount=cellcount+1;
cellid{cellcount}='mean2';

PMIN=0.01;
s=[size(predxc) 1 1 1];

if ~verbose
   return;
end

for ii=1:cellcount,
   for respidx=1:s(3),
      if respidx==1,
         fprintf('%-8s',cellid{ii});
      else
         fprintf('%-8s','');
      end
      
      for batchidx=1:s(1),
         for nlidx=1:s(2),
            if ~isnan(predxc(batchidx,nlidx,respidx,ii)),
               fprintf(' %6.3f',predxc(batchidx,nlidx,respidx,ii));
               if predp(batchidx,nlidx,respidx,ii)==-1,
                  fprintf('x');
               elseif predp(batchidx,nlidx,respidx,ii)<PMIN,
                  fprintf('*');
               else
                  fprintf(' ');
               end
            else
               fprintf(' %-6s ',' x.xxx');
            end
            
         end
      end
      
      fprintf('\n');
   end
end

if nargout<1,
   clear cellid predxc
end

return


predcount=s(1);
latcount=s(4);
nlcount=s(2);
predxc=zeros(s(1),s(2)*(cellcount-1),latcount,nlcount);
predp=zeros(s(1),s(2)*(cellcount-1),latcount,nlcount);
for ii=1:cellcount-1,
   nltemp1=min([size(preddata(ttt(ii)).predxc,1) s(1)]);
   nltemp4=min([size(preddata(ttt(ii)).predxc,4) nlcount]);
   
   predxc(1:nltemp1,ii,:,1:nltemp4)=...
       preddata(ttt(ii)).predxc(1:nltemp1,1,:,1:nltemp4);
   if ~isempty(preddata(ttt(ii)).predp),
      predp(1:nltemp1,ii,:,1:nltemp4)=...
          preddata(ttt(ii)).predp(1:nltemp1,1,:,1:nltemp4);
   else
      predp(1:nltemp1,ii,:,1:nltemp4)=-1;
   end
end

PMIN=0.01;
predxc=cat(2,predxc,permute(nanmean(permute(predxc,[2 1 3 4])),[2 1 3 4]));
predp=cat(2,predp,permute(nanmean(permute(predp,[2 1 3 4])),[2 1 3 4]));
%for predidx=1:predcount,
%   for nlidx=1:nlcount,
%      goodidx=find(predp(predidx,1:cellcount-1,1,nlidx)<PMIN);
%      predxc(predidx,cellcount,1,nlidx)=...
%          nanmean(predxc(predidx,goodidx,1,nlidx));
%   end
%end

for ii=1:cellcount,
   if 1,
      predrange=1:size(predxc,1);
   elseif ismember(batchid,[30 31 32]),
      predrange=3;
   elseif ismember(batchid,[27 28 29]),
      predrange=2;
   else
      predrange=1;
   end
   for predidx=predrange,
      if predidx==predrange(1),
         fprintf('%-6s',cellid{ii});
      else
         fprintf('%-2s','');
      end
      %nlrange=1:nlcount; % [1 2 3 4 7 9];
      if nlcount>=9,
         nlrange=[2 5 8];
         %nlrange=[1 5 9];
      else
         nlrange=1:nlcount;
      end
      for nlidx=nlrange,
         if ~isnan(predxc(predidx,ii,:,nlidx)),
            fprintf(' %6.3f',predxc(predidx,ii,:,nlidx));
	    if predp(predidx,ii,:,nlidx)==-1,
	       fprintf('x');
	    elseif predp(predidx,ii,:,nlidx)<PMIN,
	       fprintf('*');
	    else
	       fprintf(' ');
	    end
         else
            fprintf(' %-6s ',' x.xxx');
         end
      end
   end
   fprintf('\n');
end



