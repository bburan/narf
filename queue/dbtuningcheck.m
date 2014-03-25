function dbtuningcheck(cellid);
    
FORCE_TORC=1;

%baphy_set_path;
dbopen;

tt=dbReadTuning(cellid);
if isfield(tt,'strfsparse'),
    disp([cellid ' existing tuning:']);
    tt
else
    %disp('no tuning in db');
end

topo_data=struct();
cc=1;
topo_data(cc).cellid=cellid;

sql=['SELECT sCellFile.*,gCellMaster.depth,'...
     ' gPenetration.wellposition, wellimfile FROM sCellFile',...
     ' INNER JOIN gCellMaster ON gCellMaster.id=sCellFile.masterid',...
     ' INNER JOIN gPenetration ON gPenetration.id=gCellMaster.penid',...
     ' WHERE sCellFile.cellid="',cellid,'"',...
     ' AND runclassid in (1)',...
     ' ORDER BY runclassid LIMIT 1'];
cellfiledata=mysql(sql);

sql=['SELECT sCellFile.* FROM sCellFile',...
     ' WHERE sCellFile.cellid="',cellid,'"',...
     ' AND runclassid in (2)',...
     ' ORDER BY rawid DESC LIMIT 1'];
ftcfiledata=mysql(sql);
     
if ~isempty(ftcfiledata) && ~(FORCE_TORC || ~isempty(cellfiledata)),
    
    topo_data(cc).ftc_parmfile=[ftcfiledata.stimpath ftcfiledata.stimfile];
    topo_data(cc).ftc_parmfile=nslmakelocal(topo_data(cc).ftc_parmfile);
    topo_data(cc).ftc_spikefile=[ftcfiledata.path ftcfiledata.respfile];
    topo_data(cc).ftc_spikefile=nslmakelocal(topo_data(cc).ftc_spikefile);
    topo_data(cc).channum=ftcfiledata.channum;
    topo_data(cc).unit=ftcfiledata.unit;
    topo_data(cc).area=ftcfiledata.area;
    [topo_data(cc).ftc_bf,topo_data(cc).ftc_lat]=...
        ftc_tuning(topo_data(cc).ftc_spikefile,topo_data(cc).channum,...
                   topo_data(cc).unit);
    fprintf('ftc: bf: %d lat: %d\n',...
            topo_data(cc).ftc_bf,topo_data(cc).ftc_lat);
    
    % save cell information to gSingleCell
    if topo_data(cc).ftc_lat>0,
        thisbf=topo_data(cc).ftc_bf;
        thislat=topo_data(cc).ftc_lat;
        thissrcid=ftcfiledata.rawid;
    else
        thisbf=0;
        thislat=0;
        thissrcid=0;
    end
    
    tt=struct();
    tt.bf=thisbf;  % ftc (or tor if no ftc)
    tt.lat=thislat;  % ftc (or tor if no ftc)
    tt.torlat=0;
    tt.torbf=0;
    tt.bw=0;
    tt.dur=0;
    tt.rawid=thissrcid;
    tt.snr=0;
    tt.linpred=0;
    tt.strfsparse=0;
    
    disp([cellid 'new data:'])
    tt
    
    sql=sprintf(['UPDATE gSingleCell set bf=%d,bw=%.2f,latency=%d,',...
                 'duration=%d,area="%s" WHERE cellid="%s"'],...
                thisbf,0,tt.lat,0,...
                topo_data(cc).area,topo_data(cc).cellid);
    %mysql(sql);
    %dbWriteTuning(topo_data(cc).cellid,tt);
    %topo_data(cc).bad=0;
    
elseif ~isempty(cellfiledata),
    
    topo_data(cc).tor_parmfile=[cellfiledata.stimpath cellfiledata.stimfile];
    topo_data(cc).tor_parmfile=nslmakelocal(topo_data(cc).tor_parmfile);
    topo_data(cc).tor_spikefile=[cellfiledata.path cellfiledata.respfile];
    topo_data(cc).tor_spikefile=nslmakelocal(topo_data(cc).tor_spikefile);
    topo_data(cc).channum=cellfiledata.channum;
    topo_data(cc).unit=cellfiledata.unit;
    topo_data(cc).area=cellfiledata.area;
    
    % figure out which well image used for this pen
    topo_data(cc).wellimfile=cellfiledata.wellimfile;
    
    % record x,y coordinates for this site
    pos=strsep(char(cellfiledata.wellposition),'+');
    if length(pos)>cellfiledata.channum,
        xx=strsep(pos{cellfiledata.channum},',',0);
        topo_data(cc).x=xx{1};
        topo_data(cc).y=xx{2};
    else
        topo_data(cc).x=0;
        topo_data(cc).y=0;
    end
    
    [topo_data(cc).bf,topo_data(cc).bw,topo_data(cc).lat, ...
     topo_data(cc).offlat,topo_data(cc).snr,topo_data(cc).linpred,...
     topo_data(cc).strf]=...
        tor_tuning(topo_data(cc).tor_parmfile,topo_data(cc).tor_spikefile,...
                   topo_data(cc).channum,topo_data(cc).unit);
    
    fprintf('tor: bf: %d bw: %.2f lat: %d dur: %d\n',...
            topo_data(cc).bf,topo_data(cc).bw,topo_data(cc).lat, ...
            topo_data(cc).offlat-topo_data(cc).lat);
    
    if ~isempty(ftcfiledata),
        topo_data(cc).ftc_parmfile=[ftcfiledata.stimpath ftcfiledata.stimfile];
        topo_data(cc).ftc_parmfile=nslmakelocal(topo_data(cc).ftc_parmfile);
        topo_data(cc).ftc_spikefile=[ftcfiledata.path ftcfiledata.respfile];
        topo_data(cc).ftc_spikefile=nslmakelocal(topo_data(cc).ftc_spikefile);
        [topo_data(cc).ftc_bf,topo_data(cc).ftc_lat]=...
            ftc_tuning(topo_data(cc).ftc_spikefile,topo_data(cc).channum,...
                       topo_data(cc).unit);
        fprintf('ftc: bf: %d lat: %d\n',...
                topo_data(cc).ftc_bf,topo_data(cc).ftc_lat);
    else
        topo_data(cc).ftc_bf=0;
        topo_data(cc).ftc_lat=0;
        disp('no ftc data');
    end
    
    strf=topo_data(cc).strf;
    
    strfsparse=max(abs(strf(:)))./std(strf(:));
    fprintf('strf: snr: %.2f linsnr: %.2f sparse: %.2f\n',...
            topo_data(cc).snr,topo_data(cc).linpred,strfsparse);
    % save cell information to gSingleCell
    if topo_data(cc).ftc_bf>0,
        thisbf=topo_data(cc).ftc_bf;
        thislat=topo_data(cc).ftc_lat;
        thissrcid=ftcfiledata.rawid;
    else
        thisbf=topo_data(cc).bf;
        thislat=topo_data(cc).lat;
        thissrcid=cellfiledata.rawid;
    end
    
    tt=struct();
    tt.bf=thisbf;  % ftc (or tor if no ftc)
    tt.lat=thislat;  % ftc (or tor if no ftc)
    tt.torlat=topo_data(cc).lat;
    tt.torbf=topo_data(cc).bf;
    tt.bw=topo_data(cc).bw;
    tt.dur=topo_data(cc).offlat-topo_data(cc).lat;
    tt.rawid=thissrcid;
    tt.snr=topo_data(cc).snr;
    tt.linpred=topo_data(cc).linpred;
    tt.strfsparse=strfsparse;
    
    disp([cellid 'new data:'])
    tt
    
    
    sql=sprintf(['UPDATE gSingleCell set bf=%d,bw=%.2f,latency=%d,',...
                 'duration=%d,area="%s" WHERE cellid="%s"'],...
                thisbf,topo_data(cc).bw,thislat,...
                topo_data(cc).offlat-topo_data(cc).lat,...
                topo_data(cc).area,topo_data(cc).cellid);
    %mysql(sql);
    dbWriteTuning(topo_data(cc).cellid,tt,1);
    %topo_data(cc).bad=0;
         
else
    
    % disp('no data');
end
