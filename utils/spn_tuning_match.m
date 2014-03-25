function res=spn_tuning_match(cellid,batchid,verbose);
    
    if ~exist('verbose','var'),
        verbose=1;
    end
    
    dbopen;
    
    cfd=dbbatchcells(batchid,cellid);
    rawid=cfd(1).rawid;
    baphyparms=dbReadData(rawid);
    
    tt=dbReadTuning(cellid);
    if ~isfield(tt,'bnbbf'),
        tbnb=bnb_tuning(cellid);
        if isempty(tbnb),
            tbnb.bnbbf=[];
        end
        dbWriteTuning(cellid,tbnb,1);
        tt=dbReadTuning(cellid);
    end
    if ~isfield(tt,'torbf'),
        dbtuningcheck(cellid);
        tt=dbReadTuning(cellid);
    end
    if ~isfield(tt,'lof'),
        tftc=bnb_tuning(cellid,'FTC');
        if isempty(tftc),
            tftc.lof=[];
            tftc.hif=[];
        end
        dbWriteTuning(cellid,tftc,1);
        tt=dbReadTuning(cellid);
    end
        
    blo=[];
    if ~isempty(tt.lof) && tt.lof>0,
        if verbose,
            fprintf('FTC:  BF=%.0f  (%.0f-%.0f)\n',...
                    tt.bf,tt.lof,tt.hif);
        end
        
        blo=tt.lof;
        bhi=tt.hif;
    end
    if isfield(tt,'torbf') && ~isempty(tt.torbf) && tt.torbf>0,
        if tt.snr>0.12,
            % require min SNR to avoid spurious tuning
            blo=round(2.^(log2(tt.torbf)-tt.bw./2));
            bhi=round(2.^(log2(tt.torbf)+tt.bw./2));
            if verbose,
                fprintf('TOR:  BF=%.0f  (%.0f-%.0f)\n',...
                        tt.torbf,blo,bhi);
            end
        end
    end
    if ~isempty(tt.bnbbf) && tt.bnbbf>0,
        if verbose,
            fprintf('BNB:  BF=%.0f  (%.0f-%.0f)\n',...
                    tt.bnbbf,tt.bnblof,tt.bnbhif);
        end
        blo=tt.bnblof;
        bhi=tt.bnbhif;
    end
    
    if isempty(blo),
        disp('cannot calculate overlap');
        res=-1;
        return
    end
    
    lblo=log2(blo);
    lbhi=log2(bhi);
    bandcount=length(baphyparms.Ref_LowFreq);
    
    for bandidx=1:bandcount,
    
        stimlo=baphyparms.Ref_LowFreq(bandidx);
        stimhi=baphyparms.Ref_HighFreq(bandidx);
        lslo=log2(stimlo);
        lshi=log2(stimhi);
        if verbose,
            fprintf('SPN range band %d: %d-%d\n',bandidx,stimlo,stimhi);
        end
        
        % fraction overlap==portion of lof-hif encompassed by band
        
        if lshi<lblo || lslo>lbhi,
            res(bandidx)=0;
        elseif lblo>=lslo && lbhi<=lshi,
            res(bandidx)=1;
        elseif lblo<lslo && lbhi<=lshi,
            res(bandidx)=1-(lslo-lblo)./(lbhi-lblo);
        elseif lblo>=lslo && lbhi>lshi,
            res(bandidx)=1-(lbhi-lshi)./(lbhi-lblo);
        else
            res(bandidx)=1-(lshi-lslo)./(lbhi-lblo);
        end
    end
    
    
