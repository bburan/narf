

sql=['SELECT gDataRaw.*,gData.svalue FROM gDataRaw LEFT JOIN gData ON' ...
     ' (gDataRaw.id=gData.rawid AND gData.name="Ref_NoiseType")',...
     ' WHERE gDataRaw.runclass in ("SP1","SPC","SP2","SP3","SP4")',...
     ' AND not(gDataRaw.bad)',...
     ' AND isnull(gData.rawid)'];
rawdata=mysql(sql);

[defparm,defperf]=dbReadData(48246);
fieldset=fields(defparm);

for ii=1:length(rawdata),
   [parm,perf]=dbReadData(rawdata(ii).id)

   if ~isempty(parm) && ~isfield(parm,'Ref_SNR'),
      keyboard
   end
   for ff=1:length(fieldset),
      if ~isfield(parm,fieldset{ff}),
         parm=setfield(parm,fieldset{ff},getfield(defparm,fieldset{ff}));
      end
   end
   
   dbWriteData(rawdata(ii).id,parm,0);
end







files={'amethyst.mat','bom.mat','diamond.mat','jill.mat',...
       'rai.mat','salsa.mat','zim.mat'};

cd /home/svd/data/STRFlists

fid=fopen('strfbiglist.txt','w');

for ii=1:length(files),
   a=load(files{ii});
   for jj=1:length(a.FNAME),
      fprintf(fid,'%4d %s\n',a.Lfreq(jj),a.FNAME(jj).name);
      
   end
end
fclose(fid);

   
   


load /home2/vdata/david/tmp/kvaparms/kvatunesum.batch82.mat

a=cat(1,res.cartvscc);
b=cat(1,res.cartvsccmin);

a=a(:,1:3);
b=b(:,1:3);

bestcatmax=zeros(size(a,1),1);
bestcatmin=zeros(size(a,1),1);
bestcatdiff=zeros(size(a,1),1);

for ii=1:size(a,1),
   bestcatmax(ii)=min(find(a(ii,:)==max(a(ii,:))));
   bestcatmin(ii)=min(find(b(ii,:)==min(b(ii,:))));
   bestcatdiff(ii)=min(find(a(ii,:)-b(ii,:)==max(a(ii,:)-b(ii,:))));
end


