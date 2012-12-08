
clear xcparms

cellid='r0301';
speed=60;
runclassid=2;
fmt='pixel';

bootcount=10;

cellfiledata=dbgetscellfile('cellid',cellid,'speed',speed,...
                            'runclassid',runclassid,'fmt',fmt);
sql=['SELECT * FROM gCellMaster WHERE id=',...
     num2str(cellfiledata(1).masterid)];
celldata=mysql(sql);

xcparms.stimfiles={};
xcparms.respfiles={};
for ii=1:length(cellfiledata),
   xcparms.stimfiles{ii}=[cellfiledata(ii).stimpath,...
                    cellfiledata(ii).stimfile];
   xcparms.respfiles{ii}=[cellfiledata(ii).path,...
                    cellfiledata(ii).respfile];
end

[ll,iconside]=imfileinfo(xcparms.stimfiles{1});
scalepix=round(iconside(1)./celldata.rfsize*16);

xcparms.stimloadcmd='loadimfile';
xcparms.stimloadparms={scalepix,33,16};
if 1,
   xcparms.stimfiltercmd='movpower';
   xcparms.stimfilterparms={0,0,0,1,0};
   xcparms.kernfmt='pfft';
else
   xcparms.stimfiltercmd='';
   xcparms.kernfmt='space';
end
xcparms.outfile='/auto/k5/david/tmp/cellxc2out.mat';
xcparms.sfsstep=6;
xcparms.sfscount=30;
xcparms.sffiltsigma=5;
xcparms.smoothtime=2;
xcparms.docellfit2=0;

% figure out times
tparms=xcparms;
tparms.predfrac=0.0;
tparms.fitfrac=0.0;
times=xcfilefracs(tparms);

filelens=times(1).stop-times(1).start+1;
filelens(times(3).fileidx)=filelens(times(3).fileidx)+...
    times(3).stop-times(3).start+1;
totlen=sum(filelens);

bootstep=totlen/bootcount;
fstep=cumsum(filelens);
for bootidx=1:bootcount,
   bootparms(bootidx).times=times;
   valstart=
end


keyboard



xcparms.predfrac=0.1;
xcparms.fitfrac=0.1;
cellxcnodb(xcparms);


