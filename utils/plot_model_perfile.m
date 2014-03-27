function plotpath = plot_model_perfile(stack,xxx,meta,sel_results)
    
    global META XXX STACK NARF_SAVED_IMAGES_PATH;
    
    SPECIAL_FORMAT=1;
    
    if length(sel_results)>1,
        disp('can only plot one model. truncating selection');
        sel_results=sel_results(1);
    end
    
    midx=1;
    load_model(char(sel_results(midx).modelpath));
    calc_xxx(1);
    
    maxplots=0;
    jpgnames={};
    ap{midx} = [];
    modelcount=1;
    for ii = 1:length(STACK)
        m = STACK{ii}{1};
        if isfield(m, 'auto_plot') && ~isempty(m.auto_plot),
            ap{midx}(end+1) = ii;
        end
        modelcount=max(modelcount,length(STACK{ii}));
    end
    maxplots = max(maxplots,length(ap{midx})*2);
    
    
    fh=figure;
    %p=get(gcf,'Position');
    %set(gcf,'Position',[p(1)-p(3) p(2) p(3).*3 p(4)]);
    
    savefilecodes=XXX{1}.filecodes;
    ufilecodes=unique(XXX{1}.filecodes);
    
    % figure out max channel of strf
    [mods, mod_idxs] = find_modules(STACK, 'fir_filter');
    fir_sum=[];
    for ii=1:length(mods{1}),
        fir_sum=cat(2,fir_sum,mods{1}{ii}.coefs);
    end
    mm=std(fir_sum,0,2);
    maxchan=min(find(mm==max(mm)));
    firnchan=length(mm);
    depinchan=1;depoutchan=1;
    if maxchan==1
        depinchan=1, depoutchan=1
    elseif maxchan==2 && firnchan==2,
        depinchan=2, depoutchan=2
    elseif maxchan==2 && firnchan>2,
        depinchan=1, depoutchan=2
    elseif maxchan==3 && firnchan==4,
        depinchan=2, depoutchan=3,
    elseif maxchan==3 && firnchan==6,
        depinchan=1, depoutchan=1,
    end
    savetraces=cell(length(ufilecodes),1);
    
    nplots=length(ap{midx});
    clf
    for ii=1:nplots,
        idx = ap{midx}(ii);
        
        for mm = 1:modelcount,
            ts=STACK;
            for tt=1:length(ts),
                if length(ts{tt})>1,
                    ts{tt}=ts{tt}(mm);
                end
            end
            m = ts{idx}{1};
            
            ufileidx=max(find(strcmp(XXX{1}.filecodes,ufilecodes{mm})));
            fit_file=XXX{1}.training_set{ufileidx};
            test_file=XXX{1}.test_set{ufileidx};
            XXX;
            XXX{1}.filecodes=ufilecodes(mm);
            
            sfigure(fh);
            subplot(maxplots,modelcount,((ii-1).*2.*modelcount+mm));
            
            %fns = fieldnames(XXX{idx+1}.dat);
            sel.stimfile = test_file;
            sel.chan_idx = 1;
            sel.stim_idx = 1;
            
            % Call the auto-plot function
            try,
                m.auto_plot(sel, ts(1:idx), XXX(1:idx+1));
            catch
                disp('auto plot error');
            end
            XXX{1}.filecodes=savefilecodes;
            
            if ii==1 && mm==2,
                title(basename(sel_results.modelpath),'Interpreter','none');
            end
            
            %keyboard
            stiminchan=1;
            stimoutchan=1;
            if strcmpi(m.name,'depression_filter_bank'),
                stiminchan=depinchan;
                stinoutchan=depoutchan;
            elseif strcmpi(m.name,'fir_filter'),
                stiminchan=maxchan;
            end
            
            subplot(maxplots,modelcount,(((ii-1).*2+1).*modelcount+mm));
            stimin=XXX{idx}.dat.(test_file).stim(:,1,stiminchan);
            stimout=XXX{idx+1}.dat.(test_file).stim(:,1,stimoutchan);
            respout=XXX{idx+1}.dat.(test_file).respavg(:,1);
            
            fs=STACK{1}{1}.raw_resp_fs;
            
            plen=min(round(STACK{1}{1}.raw_resp_fs.*2.0),size(stimin,1));
            if SPECIAL_FORMAT,
                stimin=nanmean(stimin(1:plen,:),2);
                stimout=nanmean(stimout(1:plen,:),2);
                respout=respout(1:plen);
            end
            
            if strcmpi(m.name,'depression_filter_bank'),
                savetraces{mm}=[stimin stimout];
            elseif strcmpi(m.name,'fir_filter'),
                stimout(stimout<0)=0;
                stimin=stimin./max(stimin(:)).*max(stimout(:));
                savetraces{mm}=[savetraces{mm} stimout respout];
            end
            
            tt=(1:size(stimin,1))./fs;
            plot(tt,nanmean(stimout,2),'r','LineWidth',1);
            hold on
            plot(tt,nanmean(stimin,2),'k','LineWidth',1);
            axis tight;
            hold off
            if (midx==1),
                ylabel('resp');
            end
        end
    end
    
    for mm=1:modelcount,
        subplot(maxplots,modelcount,(maxplots-1).*modelcount+mm);
        dd=savetraces{mm};
        ddall=mean(cat(3,savetraces{:}),3);
        
        dd(:,4)=gsmooth(dd(:,4),1);
        dd(:,1)=dd(:,1)./max(ddall(:,1)).*max(ddall(:,3));
        dd(:,2)=dd(:,2)./max(ddall(:,2)).*max(ddall(:,3));
        plot(dd);
    end
    fullpage landscape
    
    return
    
% $$$         % plot pred/resp PSTHs
% $$$         subplot(maxplots,modelcount,((maxplots-2).*modelcount+midx));
% $$$         fs=STACK{1}{1}.raw_resp_fs;
% $$$         % standard baphy data, plot val
% $$$         ts=XXX{end}.test_set{1};
% $$$         dat=XXX{end}.dat.(ts);
% $$$         rr=dat.respavg(:,1);
% $$$         pp=dat.stim(:,1);
% $$$             
% $$$         plen=min(round(STACK{1}{1}.raw_resp_fs.*3.5),...
% $$$                  size(rr,1));
% $$$         if SPECIAL_FORMAT,
% $$$             rr=rr(1:plen,:);
% $$$             pp=pp(1:plen,:);
% $$$         end
% $$$             
% $$$         tt=(1:size(rr,1))./fs;
% $$$         plot(tt,nanmean(rr,2),'r','LineWidth',2);
% $$$         hold on
% $$$         plot(tt,nanmean(pp,2),'k','LineWidth',2);
% $$$         axis tight;
% $$$         hold off
% $$$         if (midx==1),
% $$$             ylabel('resp');
% $$$         end
% $$$         
% $$$         % plot stim in bottom row
% $$$         subplot(maxplots,modelcount,((maxplots-1).*modelcount+midx));
% $$$         ts=XXX{2}.test_set{1};
% $$$         dat=XXX{2}.dat.(ts);
% $$$         sstim=squeeze(dat.stim(:,1,:));
% $$$         sstim=sstim(1:plen,:);
% $$$         plot(tt,sstim,'LineWidth',2);
% $$$         if (midx==1),
% $$$             xlabel('time (sec)');
% $$$             ylabel('stim');
% $$$         end
% $$$     end
% $$$     
    % TEXT AT TOP
    if ~isfield(META, 'batch')
        META.batch = 0;
    end
    
    % Print the text at the top
    sfigure(fh);
    %axtt = axes('Parent', fig , 'Units','pixels', ...
    %            'Position', [lb nplots*ph+bb w-lb-rb th-bb]);
    subplot(maxplots,1,1);
    set(gca,'xtick',[]); set(gca,'xticklabel',[]);
    set(gca,'ytick',[]); set(gca,'yticklabel',[]);
    ts=sprintf('Batch:     %d\nCellid:    %s\nModel:     %s\nTrain Set: %s\nTest Set:  %s\nTrain r:   %.3f\nTest r:    %.3f {\\pm} %.3f', ...
               META.batch, strrep(XXX{1}.cellid,'_','\_'),...
               strrep(META.modelname,'_','\_'), ...
               strrep([XXX{1}.training_set{:}],'_','\_'),...
               strrep([XXX{1}.test_set{:}],'_','\_'), ...
               XXX{end}.score_train_corr, ...
               XXX{end}.score_test_corr, ...
               XXX{end}.score_test_floorcorr);
    ax_text  = text('FontSize', 6,'Position', [0.0, 0.45], ...
                    'String', ts);
    axis off
    
     
    h=get(gcf,'Children');
    hh=get(h,'Children');
    hh=cat(1,hh{:});
    hhh=get(hh,'Children');
    hhh=cat(1,hhh{:});
    hh=cat(1,h,hh,hhh);
    
    for hhi=1:length(hh),
        if strcmpi(get(hh(hhi),'Type'),'text'),
            set(hh(hhi),'FontSize',6);
        end
        if strcmpi(get(hh(hhi),'Type'),'axes'),
            %set(hh(hhi),'XTickLabelMode','manual');
            set(get(hh(hhi),'XLabel'),'FontSize',6);
            %set(hh(hhi),'YTickLabelMode','manual');
            set(get(hh(hhi),'YLabel'),'FontSize',6);
            set(get(hh(hhi),'ZLabel'),'FontSize',6);
            set(get(hh(hhi),'Title'),'FontSize',6);
        end
   end
   fullpage landscape
   
   if SPECIAL_FORMAT
       disp('to save eps:');
       fprintf('print -depsc -f%d %s.modelcomp.eps\n',fh,XXX{1}.cellid);
       
       if ~isempty(jpgnames),
           fprintf('%s\n',jpgnames{:});
       end
   else
       disp('to save pdf:');
       fprintf('print -dpdf -f%d %s.modelcomp.pdf\n',fh,XXX{1}.cellid);
   end