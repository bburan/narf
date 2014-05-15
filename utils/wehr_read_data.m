function [r,s,onsettimes,secPerSegment,esequence]=wehr_read_data(f,vdim,SR);

global NARF_DEBUG
envpath='/auto/data/daq/wehr/soundfiles/sourcefiles/';


fname=basename(f);

if vdim==0,
    vstring='i-clamp';
elseif vdim==1,
    vstring='v-clamp_E';
elseif vdim==2,
    vstring='v-clamp_I';
end

if strcmp(computer,'PCWIN'),
    ppdir=tempdir;
elseif exist('/auto/data/tmp/tstim/','dir')
    ppdir=['/auto/data/tmp/tstim/'];
else
    ppdir=['/tmp/'];
end
cachefilename=sprintf('%sload_stim_resps_wehr_%s_%s_fs%d.mat',...
                      ppdir,fname,vstring,SR);
if exist(cachefilename,'file') && (isempty(NARF_DEBUG) || ~NARF_DEBUG),
    fprintf('loading cached file %s\n',cachefilename);
    load(cachefilename);
else
    load(f);
    
    pp=fileparts(f);
    envfile=strrep(out.epochfilenames{1},'\',filesep);
    envfile=strrep(basename(envfile),'sourcefile','envfile');
    e=load([envpath envfile]);
    
    SRint=e.EnvSamplingRate;
    secPerSegment=size(e.EnvSet,1)./length(e.sequence)./SRint;
    
    s0=zeros(SRint,size(e.EnvSet,2));
    
    % stim/resp data will be loaded into s/r matrices
    % SR Hz sampling rate for both
    s=[];r=[];onsettimes=[];
    setcount=size(out.M1,1);
    for ii=1:setcount,
        
        % average response across repetitions
        if size(out.M1,4)>1,
            % v-clamp data
            % dim 2 index specifies which v-vclamp mode
            %ri=squeeze(mean(out.M1(ii,vdim,1:out.nreps(ii,vdim),:),3));
            ri=squeeze(out.mM1(ii,vdim,:));
            if vdim==1,
                % flip sign for E current so that up="more"
                ri=-ri;
            end
            % load low-res stimulus just to make sure it's aligned with
            % the envelope
            dsi=squeeze(mean(out.M1stim(ii,1,1:out.nreps(ii),:),3));
            sentencesperoffset=size(out.sequences,4);
        else
            % i-clamp data
            ri=squeeze(mean(out.M1(ii,:,:),2));
            %ri=squeeze(out.mM1(ii,:))';
            dsi=squeeze(mean(out.M1stim(ii,1:out.nreps(ii),:),2));
            sentencesperoffset=size(out.sequences,3);
        end
        ri=resample(ri,SRint,out.samprate);
        dsi=resample(dsi,SRint,out.samprate);
        
        % number of SRint hz samples per segment (segment==sentence?)
        offset=(secPerSegment*SRint)*sentencesperoffset;
        %ii
        %keyboard
        if ii<size(out.M1,1) || offset.*ii<size(e.EnvSet,1),
            si=[s0;e.EnvSet(offset.*(ii-1)+(1:offset),:);s0];
        else
            si=[s0;e.EnvSet((offset*(ii-1)):end,:);s0];
        end
        
        if SR<SRint,
            ri=resample(ri,SR,SRint);
            dsi=resample(dsi,SR,SRint);
            si0=si;
            si=[];
            for jj=1:size(si0,2);
                tsi=resample(si0(:,jj),SR,SRint);
                si=[si tsi];
            end
        end
        
        % trim possible artifact, last 0.75 sec of r
        ri=ri(1:(end-round(0.75*SR)),:);
        
        % make sure stim and resp are the same length
        if length(ri)>length(si),
            ri=ri(1:length(si),:);
            dsi=dsi(1:length(si),:);
        else
            si=si(1:length(ri),:);
            dsi=dsi(1:length(ri),:);
        end
        
        %remove possible artifacts from 0.1 sec onset
        ri=ri(round(0.1*SR)+1:end);
        si=si(round(0.1*SR)+1:end,:);
        dsi=dsi(round(0.1*SR)+1:end);
        
        if 1,
            % remove trend by subtracting min-filtered version
            % of signal
            N=round(SR./2);
            sy=gsmooth(ri(:),N./2);
            y=zeros(size(sy(:)));
            for ii=1:length(ri(:));
                xx=max(1,ii-N);
                yy=min(length(sy(:)),ii+N);
                y(ii)=min(sy(xx:yy));
            end
            y=gsmooth(y,N./2);
            
            if NARF_DEBUG,
                if isempty(NARF_DEBUG_FIGURE),
                    NARF_DEBUG_FIGURE=figure;
                end
                sfigure(NARF_DEBUG_FIGURE);
                clf
                subplot(2,1,1);
                plot([ri(:) sy ri(:)-y y]);
                axis tight
                subplot(2,1,2);
                ty=ri(:)-y;
                ts=nanmean(si,2);
                plot([ts./max(ts) ty./max(ty)]);
                axis tight
                pause(0.1);
                %keyboard
            end
            ri=ri(:)-y;
            ri=sqrt(abs(ri)).*sign(ri);
            ri=ri;
        elseif 0,
            % remove linear trend
            ri=detrend(ri);
        else
            % remove trend by substracting 5th order polynomial fit
            order = 5;
            %p = polyfit((1:numel(ri))', ri(:), order);
            p = polyfit((1:numel(ri))', gsmooth(ri(:),50), order);
            %plot([ri(:)  gsmooth(ri(:),100) polyval(p, (1:numel(ri))')]);
            
            if NARF_DEBUG,
                if isempty(NARF_DEBUG_FIGURE),
                    NARF_DEBUG_FIGURE=figure;
                end
                sfigure(NARF_DEBUG_FIGURE);
                clf
                plot([ri(:) polyval(p, (1:numel(ri))') ...
                      ri(:) - polyval(p, (1:numel(ri))')]);
                keyboard;
            end
            
            ri = ri(:) - polyval(p, (1:numel(ri))');
        end
        
        onsettimesi=(0:floor(length(si)./(secPerSegment*SR)-1))*...
            (secPerSegment.*SR)+1+round(length(s0)./SRint*SR)-10;%-f_idx;
        onsettimes=[onsettimes onsettimesi+length(r(:))];
        
        if length(si)<size(s,1),
            sdiff=size(s,1)-length(si);
            si=cat(1,si,ones(sdiff,size(si,2)).*nan);
            ri=cat(1,ri,ones(sdiff,size(ri,2)).*nan);
        end
        si=permute(si,[1 3 2]);
        
        r=[r ri];
        s=[s si];
    end
    fprintf('Saving cache file %s\n',cachefilename);
    esequence=e.sequence;
    save(cachefilename,'r','s','onsettimes','secPerSegment',...
         'esequence');
end


