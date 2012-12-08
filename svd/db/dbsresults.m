% function res=dbsreults(batchid)
%
% created SVD 2011-05-13 ripped off dbresv1.m
%
function res=dbsresults(batchid)

dbopen;

batstr='(';
for ii=1:length(batchid),
   batstr=[batstr,num2str(batchid(ii)),','];
end
batstr(end)=')';
sql=['SELECT distinct cellid',...
     ' FROM sRunData WHERE sRunData.batch in ',batstr,...
     ' AND (sRunData.batch in (127,128) OR not(cellid like "%model%"))',...
     ' ORDER BY cellid'];
celldata=mysql(sql);

cellids={celldata.cellid};
cellcount=length(cellids);
batcount=length(batchid);

res=[];
for batidx=1:length(batchid),
   
   sql=['SELECT sResults.*,masterid,cellid',...
        ' FROM sResults INNER JOIN sRunData ON sResults.runid=sRunData.id',...
        ' WHERE sResults.batch=',num2str(batchid(batidx)),...
        ' AND (sRunData.batch in (127,128) OR not(cellid like "%model%"))',...
        ' ORDER BY cellid'];
   resdata=mysql(sql);
   
   ttt=cat(1,resdata.runid);
   
   for ii=1:length(resdata),
      eval(char(resdata(ii).matstr));
      tp=preddata(ttt(ii));
      
      ff=fields(tp);
      cidx=find(strcmp(tp.cellid,cellids));
      
      for jj=1:length(ff),
         if strcmp(ff{jj},'cellid'),
            % skip
         else
            ss=[size(tp.(ff{jj}),1) size(tp.(ff{jj}),2) ...
                size(tp.(ff{jj}),4) size(tp.(ff{jj}),4)];
            if ii==1 && batidx==1,

               res.(ff{jj})=zeros([ss cellcount batcount]).*nan;
            end
            if size(tp.(ff{jj}),1)>size(res.(ff{jj}),1),
               tss=ss;
               tss(1)=tss(1)-size(res.(ff{jj}),1);
               res.(ff{jj})=cat(1,res.(ff{jj}),...
                                zeros([tss cellcount batcount]).*nan);
            end
            if ~isempty(tp.(ff{jj})),
               res.(ff{jj})(1:ss(1),:,:,:,cidx,batidx)=tp.(ff{jj});
            end
         end
      end
      
   end
   
end

