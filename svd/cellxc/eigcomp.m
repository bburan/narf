% eigcomp.m:  compute eigenvectors per fixation of review response data
%
% created SVD 10/08/03 -- ripped off of xcfittune
%
% function [eigmatrix,eigval]=eigcomp(cellid)
%
function [eigmatrix,eigval]=eigcomp(cellid)

dbopen;

% figure out what files to use. batch 24 is review
batch=24;
[cellfiledata,times,batchdata]=cellfiletimes(cellid,batch);

checkpath='/auto/k2/share/data/movfixdat/';
eigpath='/auto/k2/share/data/eigdata/';

movThresh=0.75;

STARTLAG=1;
STOPLAG=20;
opts.quiet=1;
opts.cattrials=1;
opts.usenan=1;

bresp=[];

for respfileidx=1:length(cellfiledata),
   
   % find fixations
   [FirstFr,LastFr]=movgetcheck(checkpath,cellid,...
                                cellfiledata(respfileidx).stimfile,movThresh);
   
   % load respfile
   respfile=[cellfiledata(respfileidx).path,...
             cellfiledata(respfileidx).respfile];
   raw_resp=respload(respfile,'',1,1,1);
   
   %if FirstFr(1)<3,
   %   FirstFr=FirstFr(2:end);
   %   LastFr=LastFr(2:end);
   %end
   %[rpsth,rraw,okframeidx]=respgetepochs(raw_resp,FirstFr-2,LastFr,1,opts);
   
   [rpsth,rraw,okframeidx]=respgetepochs(raw_resp,FirstFr,LastFr,1,opts);
   
   if STOPLAG > size(rpsth,1),
      extrapad=STOPLAG-size(rpsth,1);
   else
      extrapad=0;
   end
   
   rpsth(find(isnan(rpsth)))=-1;
   rpsth=padpsth(rpsth,extrapad);
   rpsth(find(rpsth==-1))=nan;
   
   rpsth=rpsth(STARTLAG:STOPLAG,:);
   
   bresp=[bresp rpsth];
   
end

nanidx=sum(isnan(bresp),1)==0;
bresp=bresp(:,nanidx);

figure(1);
[eigmatrix,eigval]=pca(bresp,4);
eigval=diag(eigval);

for ii=1:size(eigmatrix,1),
   if sum(eigmatrix(:,ii))<0,
      eigmatrix(:,ii)=-eigmatrix(:,ii);
   end
end


outfile=sprintf('%s%s.eig.mat',eigpath,cellid);
save(outfile,'eigmatrix','eigval');



   
   
