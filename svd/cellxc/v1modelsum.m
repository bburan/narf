% v1modelsum.m

% callmaster : exp X in-model X cell
% predxc: exp X in-model X cell X out-model X cnf

dbopen;

batches=[24 37 26 50];   % pix fft power fft+power
clear bigdata preddata
%close all
for ii=1:2,
   figure(ii);
end
drawnow;

% load pred results for each batch
expcount=size(batches,1);
domcount=size(batches,2);
cellmaster=zeros(expcount,domcount,200);
totlen=zeros(expcount,domcount,200);
for expidx=1:expcount,
   for domidx=1:domcount,
      
      % load parms for current batch
      sql=['SELECT * FROM sBatch WHERE id=',num2str(batches(expidx,domidx))];
      batchdata=mysql(sql);
      stimfmtcode=batchdata.stimfmtcode;
      respfmtcode=batchdata.respfmtcode;
      
      sql=['SELECT sResults.*,masterid,cellid',...
           ' FROM sResults INNER JOIN sRunData',...
           ' ON sResults.runid=sRunData.id',...
           ' WHERE sResults.batch=',num2str(batches(expidx,domidx)),...
           ' AND masterid>0'];
      resdata=mysql(sql);
      
      for runidx=1:length(resdata),
         
         % record prediction results for current run
         eval(char(resdata(runidx).matstr));
         
         % this could be made slicker!
         bigdata(expidx,domidx,resdata(runidx).masterid).predxc=...
             preddata(resdata(runidx).runid).predxc;
         
         bigdata(expidx,domidx,resdata(runidx).masterid).predp=...
             preddata(resdata(runidx).runid).predp;
         
         bigdata(expidx,domidx,resdata(runidx).masterid).cellid=...
             preddata(resdata(runidx).runid).cellid;
         
         bigdata(expidx,domidx,resdata(runidx).masterid).runid=...
             resdata(runidx).runid;
         
         cellmaster(expidx,domidx,resdata(runidx).masterid)=1;
         
         cellid=resdata(runidx).cellid;
         %figure out interesting parameters of run (eg, total file length)
         sql=['SELECT * FROM sCellFile',...
              ' WHERE cellid="',cellid,'"',...
              ' AND runclassid=',num2str(batchdata.runclassid),...
              ' AND stimspeedid=',num2str(batchdata.stimspeedid),...
              ' AND stimfmtcode=',num2str(stimfmtcode),...
              ' AND respfmtcode=',num2str(respfmtcode),...
              ' ORDER BY resplen'];
         cellfiledata=mysql(sql);
         bigdata(expidx,domidx,resdata(runidx).masterid).totlen=...
             sum(cat(1,cellfiledata.resplen));
         bigdata(expidx,domidx,resdata(runidx).masterid).spikes=...
             sum(cat(1,cellfiledata.spikes));
         bigdata(expidx,domidx,resdata(runidx).masterid).meanrate=...
             bigdata(expidx,domidx,resdata(runidx).masterid).spikes./...
             bigdata(expidx,domidx,resdata(runidx).masterid).totlen ...
             .* 1000/14;
         %fprintf('%s (%d/%d): rate=%.1f Hz\n',...
         %        cellid,expidx,domidx,...
         %        bigdata(expidx,domidx,resdata(runidx).masterid).meanrate);
      end
   end
end

% don't need this any more?
clear preddata

fprintf('bigdata matrix loaded expcount=%d domcount=%d\n',expcount,domcount);

goodcells=find(squeeze(sum(sum(cellmaster,1),2)));
cellcount=length(goodcells);
cellmaster=cellmaster(:,:,goodcells);
bigdata=bigdata(:,:,goodcells);

ii=min(find(cellmaster(1,1,:)));
nlcount=9;
nluse=1:nlcount;

% predxc: cell X expclass X in-model X cnfclass X out-model
predxc=ones(cellcount,expcount,domcount,expcount,nlcount).*nan;
predp=ones(cellcount,expcount,domcount,expcount,nlcount);
meanrate=zeros(cellcount,expcount).*nan;
totlen=ones(cellcount,expcount).*nan;
spikes=ones(cellcount,expcount).*nan;
celllist=cell(cellcount,1);
for expidx=1:expcount,
   for domidx=1:domcount,
      tgoodidx=find(cellmaster(expidx,domidx,:));
      tgoodidx=tgoodidx(:)';
      for cellidx=tgoodidx,
         cnfcount=size(bigdata(expidx,domidx,cellidx).predxc,1);
         celllist{cellidx}=bigdata(expidx,domidx,cellidx).cellid;
         predxc(cellidx,expidx,domidx,1:cnfcount,:)=...
             reshape(bigdata(expidx,domidx,cellidx).predxc(:,1,1,nluse),...
                     [1 1 1 cnfcount nlcount]);
         if ~isempty(bigdata(expidx,domidx,cellidx).predp),
            predp(cellidx,expidx,domidx,1:cnfcount,:)=...
                reshape(bigdata(expidx,domidx,cellidx).predp(:,1,1,nluse),...
                        [1 1 1 cnfcount nlcount]);
         end
         if domidx==1,
            totlen(cellidx,expidx)=bigdata(expidx,domidx,cellidx).totlen;
            spikes(cellidx,expidx)=bigdata(expidx,domidx,cellidx).spikes;
            meanrate(cellidx,expidx)=bigdata(expidx,domidx,cellidx).meanrate;
         end
      end
   end
end
%longrevidx=find(totlen(:,1,1)>1000);

MINLEN=500;
PMIN=0.01;

expstr={'pix','fft','pow','lin2'};
instr={'pix','pow','psf'};
cnfstr=expstr;
outstr={'lin','rec','full','fitrec1','fitrec2','fitrec3',...
        'fulrec1','fulrec2','fulrec3'};

figure(1);
clf
exprange1=[1 2; 2 3];
exprange2=[2 3; 4 4];
cnfrange=[1 1; 1 1];
nlrange=[5 5; 5 5];

rowcount=size(exprange1,1);
colcount=size(exprange1,2);
for ii=1:length(exprange1(:)),
   subplot(rowcount,colcount,ii);
   n1=[expstr{exprange1(ii)}];
   n2=[expstr{exprange2(ii)}];
   set1=predxc(:,1,exprange1(ii),cnfrange(ii),nlrange(ii));
   set2=predxc(:,1,exprange2(ii),cnfrange(ii),nlrange(ii));
   p1=predp(:,1,exprange1(ii),cnfrange(ii),nlrange(ii));
   p2=predp(:,1,exprange2(ii),cnfrange(ii),nlrange(ii));
   goodidx=find(~isnan(set1) & ~isnan(set2));
   sigidx=2-(p1(goodidx)<PMIN | p2(goodidx)<PMIN);
   
   plotcomp(set1(goodidx),set2(goodidx),n1,n2,[-0.4 1.0 -0.4 1.0],sigidx);
end

set(gcf,'PaperPosition',[0.25 0.25 8 10.5],'PaperOrientation','portrait');
drawnow


