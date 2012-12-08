% function r=tnlpart2(cellid,batch)
%
% load results from resfile and do separable time fit
%
% r=0 if no entries found in db, =1 otherwise
%
function r=tnlpart2(runidx,batch)

if ~exist('batch','var'),
   batch=80;
   fprintf('assuming batch=80\n');
end

dbopen;
if isnumeric(runidx),
   sql=['SELECT * from sRunData WHERE id=',num2str(runidx)];
   rundata=mysql(sql);
   cellid=rundata.cellid;
else
   cellid=runidx;
   goodbatch=zeros(1,length(batch));
   sql=['SELECT * from sRunData WHERE cellid="',cellid,'"',...
        ' AND batch=',num2str(batch)];
   rundata=mysql(sql);
end

if length(rundata)==0,
   disp('no entries found in db!');
   if nargout>0,
      r=0;
   end
   return
end

outfile=[rundata(1).respath,rundata(1).resfile,'.gz'];
fprintf('outfile: %s\n',outfile);

if strcmp(outfile(end-2:end),'.gz'),
   zload(outfile);
else
   load(outfile);
end



% rest parms to match full time RC
params.maxlag=[-6 13];
params.resploadcmd='respload';
params.resploadparms={'',1,1,1};
params.stimloadcmd='loadimfile';
params.stimloadparms={0,0,16,0,1};

[cellfiledata,times,tparms]=cellfiletimes(cellid,26);
params.times=times;
starttimes=times(1).start;
stoptimes=times(1).stop;
fitfile=times(2).fileidx;
fitstartframe=times(2).start;
fitstopframe=times(2).stop;
predfile=times(3).fileidx;
predstartframe=times(3).start;
predstopframe=times(3).stop;

attcount=size(times(1).start,2);  % ignore attcount>1 for the time being
attidx=1;

xcloadfiles;

% pull out spatial kernel...
h=strf(1).h(:,3);
mS=strf(1).mS;
savestrf=strf;

hsep=[h h.*(h>0) h.*(h<0)];


% project stim onto kernel, positve only and negative only

movlen=size(stim,1);
stimbak=stim;
sproj=(stim-repmat(mS',[movlen 1]))*hsep;

% do rc on just time proj
firstseg=1;
stim=sproj(:,1);
spacecount=1;
params.sfscount=1;
xccore;


% do rc on pos/neg sep projections
firstseg=1;
stim=sproj(:,2:3);
spacecount=2;
params.sfscount=1;
params.maxlag=[-6 13];
xccore;

for bootidx=1:params.resampcount,
   for sfsidx=1:size(H,3),
      H(:,:,sfsidx,1,bootidx)=...
          sSA2(:,:,1,bootidx)^-1 * H(:,:,sfsidx,1,bootidx);
   end
end

clear strf
params.nloutparm=1;  % only do linear fit here
xcfit2;

plot(strf.h');

keyboard

