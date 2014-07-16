function res=bnb_tuning(cellid,runclass)
    
    dbopen;
    
    if ~exist('runclass','var'),
        runclass='BNB';
    end
    
    res=[];
    cfd=dbgetscellfile('cellid',cellid,'runclass',runclass);
    usesorted=1;
    if isempty(cfd),
        sql=['SELECT gDataRaw.*,gSingleCell.channum',...
             ' FROM gDataRaw INNER JOIN gSingleCell',...
             ' ON gDataRaw.masterid=gSingleCell.masterid',...
             ' WHERE runclass="',runclass,'"',...
             ' AND gSingleCell.cellid="',cellid,'"',...
             ' AND not(gDataRaw.bad)'];
        rawdata=mysql(sql);
        if ~isempty(rawdata),
            cfd=struct;
            cfd.stimpath=rawdata(1).resppath;
            cfd.stimfile=rawdata(1).parmfile;
            cfd.channum=rawdata(1).channum;
            cfd.unit=1;
        end
        usesorted=0;
    end
    
    if isempty(cfd),
        return;
    end
    
    parmfile=[cfd(1).stimpath cfd(1).stimfile];
    
    fprintf('bnb_tuning: cell %s file %s\n',cellid,basename(parmfile));
    LoadMFile(parmfile);
    
    options.PreStimSilence=exptparams.TrialObject.ReferenceHandle.PreStimSilence;
    options.PostStimSilence=exptparams.TrialObject.ReferenceHandle.PostStimSilence;
    options.rasterfs=100;
    options.channel=cfd(1).channum;
    options.unit=cfd(1).unit;
    options.usesorted=usesorted;
    options.datause='Ref Only';
    [r,tags]=raster_load(parmfile,options.channel,options.unit,options);
    
    fset=[];
    for ii=1:length(tags),
        tt=strsep(tags{ii},',',1);
        fset=cat(1,fset,str2num(tt{2}));
    end
    [unique_freq,~,mapidx]=unique(fset);
    freqcount=length(unique_freq);
    
    TrialCount=sum(~isnan(r(1,:)));
    MaxFreq=round(min(TrialCount./4,30));
    
    if freqcount>MaxFreq,
        newidx=round(linspace(1,freqcount+1,MaxFreq+1));
        newfreq=zeros(MaxFreq,1);
        for jj=1:MaxFreq,
            newfreq(jj)=round(2.^mean(log2(unique_freq(newidx(jj):(newidx(jj+1)-1)))));
            ff=find(ismember(fset,unique_freq(newidx(jj):(newidx(jj+1)-1))));
            fset(ff)=newfreq(jj);
        end
        [unique_freq,~,mapidx]=unique(fset);
        freqcount=length(unique_freq);
    end
    
    b=zeros(freqcount,1);
    p=zeros(freqcount,1);
    pe=zeros(freqcount,1);
    bbase=1:round(options.PreStimSilence.*options.rasterfs+1);
    bb=round(options.PreStimSilence.*options.rasterfs+1):...
       round(size(r,1)-options.PostStimSilence.*options.rasterfs);
    baseline=nanmean(nanmean(r(1:bb(1),:)));
    
    for ii=1:freqcount,
        tr=r(bbase,:,find(mapidx==ii));
        b(ii)=nanmean(tr(:));
        tr=r(bb,:,find(mapidx==ii));
        p(ii)=nanmean(tr(:))-b(ii)+baseline;
        pe(ii)=nanstd(tr(:))./sqrt(sum(~isnan(tr(:))));
    end
    
    ps=gsmooth(p,freqcount./40);
    
    mm=min(find(ps==max(ps)));
    if p(mm)-2.*pe(mm)>baseline,
        bf=unique_freq(mm);
        
        % resample to find lo and hi bands
        xx=1:0.5:freqcount;
        prs=interp1(1:freqcount,ps,xx);
        pers=interp1(1:freqcount,pe,xx);
        freqrs=round(2.^interp1(1:freqcount,log2(unique_freq),xx));
        mm2=min(find(prs==max(prs)));
        
        lidx=mm2;
        while lidx>1 && prs(lidx-1)-2.*pers(lidx)>baseline &&...
                prs(lidx-1)>=prs(mm2)./2,
            lidx=lidx-1;
        end
        hidx=mm2;
        while hidx<length(prs) && prs(hidx+1)-2.*pers(hidx+1)>baseline &&...
                prs(hidx+1)>=prs(mm2)./2
            hidx=hidx+1;
        end
        blo=freqrs(lidx);
        bhi=freqrs(hidx);
        
    else
        bf=0;
        blo=0;
        bhi=0;
    end
    
    if strcmpi(runclass,'BNB'),
        res.bnbbf=bf;
        res.bnblof=blo;
        res.bnbhif=bhi;
    elseif strcmpi(runclass,'FTC'),
        res.bf=bf;
        res.lof=blo;
        res.hif=bhi;
    end
    
    figure(1);
    clf
    subplot(2,1,1);
    raster_plot(parmfile,r,tags,gca,options);
    
    subplot(2,1,2);
    plot([0 freqcount+1],[0 0]+baseline,'k--');
    hold on
    errorbar(p,pe);
    if bf>0,
        plot(mm,p(mm),'ro');
    end
    hold off
    
    stimulusvalues=unique_freq;
    stimuluscount=length(stimulusvalues);
    set(gca,'XLim',[0 stimuluscount+1]);  % make x axis look nice
    if max(stimulusvalues)>4000,
        stimulusKHz=round(stimulusvalues./100)./10;
    elseif max(stimulusvalues)>100
        stimulusKHz=round(stimulusvalues./10)./100;
    else
        stimulusKHz=stimulusvalues;
    end
    stimlabelidx=1:2:stimuluscount;
    set(gca,'XTick',stimlabelidx,'XTickLabel',stimulusKHz(stimlabelidx));
    
    %axis([0 freqcount+1 aa(3:4)]);
    if bf>0,
        ht=title(sprintf('%s Rep%d BF=%.0f (%d-%d)',...
                         cellid,size(r,2),bf,blo,bhi));
    else
        ht=title(sprintf('%s Rep%d BF=nan',...
                         cellid,size(r,2)));
    end
    set(ht,'Interpreter','none');
    set(gcf,'Name',sprintf('%s(%d)',basename(parmfile),size(r,2)));
    
    drawnow

