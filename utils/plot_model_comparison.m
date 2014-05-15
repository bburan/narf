function plotpath = plot_model_comparison(stack,xxx,meta,sel_results)
    
    global META XXX STACK NARF_SAVED_IMAGES_PATH;
    
    SPECIAL_FORMAT=1;
    PLOT_STIM=0;
    
    if length(sel_results)>5,
        disp('max comparison is 5 models, truncating list');
        sel_results=sel_results(1:5);
    end
    
    modelcount=length(sel_results);
    
    stacks=cell(modelcount,1);
    xxxs=cell(modelcount,1);
    metas=cell(modelcount,1);
    maxplots=0;
    jpgnames={};
    for midx=1:modelcount,
        load_model(char(sel_results(midx).modelpath));
        stacks{midx}=STACK;
        xxxs{midx}=XXX;
        metas{midx}=META;

        % Scan through the STACK looking for things to .auto_plot
        ap{midx} = [];
        for ii = 1:length(STACK)
            m = STACK{ii}{1};
            if isfield(m, 'auto_plot') && ~isempty(m.auto_plot) && ...
                    ~strcmpi(m.name,'correlation') && ...
                    ~strcmpi(m.name,'mean_squared_error'),
                ap{midx}(end+1) = ii;
            end
        end
    end
    %length(ap{midx})
    maxplots = length(ap{midx})+2+PLOT_STIM;

    
    fh=figure;
    p=get(gcf,'Position');
    set(gcf,'Position',[p(1)-p(3) p(2) p(3).*3 p(4)]);
    
    for midx=1:modelcount,
        STACK=stacks{midx};
        XXX=xxxs{midx};
        META=metas{midx};
        calc_xxx(1);
        nplots=length(ap{midx});
        
        % Call the auto-plot functions
        for ii = 1:nplots
            idx = ap{midx}(ii);
            m = STACK{idx}{1};
            
            plotfn = m.auto_plot;
            
            sfigure(fh);
            subplot(maxplots,modelcount,(ii.*modelcount+midx));
            %ax = axes('Parent', fig, 'Units', 'pixels', ...
            %          'Position', [lb (nplots-ii)*ph+bb w-lb-rb ph-bb]);
            
            fns = fieldnames(XXX{idx+1}.dat);
            sel.stimfile = fns{1};
            sel.chan_idx = 1;
            sel.stim_idx = 1;
            plotfn(sel, STACK(1:idx), XXX(1:idx+1));
            
            if 0 && SPECIAL_FORMAT && strcmp(m.name,'fir_filter'),
                fhs=figure;
                mmm=max(abs(m.coefs(:)));
                imagesc(m.coefs,[-mmm mmm]);
                axis off
                axis image xy
                jpgnames{length(jpgnames)+1}=...
                    sprintf('fig_print(%d,''%s.strfcomp.%s.png'',[],1,10)',...
                            fhs,XXX{1}.cellid,META.modelname);
            end
        end
        
        % TEXT AT TOP
        if ~isfield(META, 'batch')
            META.batch = 0;
        end
        
        % Print the text at the top
        sfigure(fh);
        %axtt = axes('Parent', fig , 'Units','pixels', ...
        %            'Position', [lb nplots*ph+bb w-lb-rb th-bb]);
        subplot(maxplots,modelcount,midx);
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
        subplot(maxplots,modelcount,((maxplots-1-PLOT_STIM).*modelcount+midx));
        fs=STACK{1}{1}.raw_resp_fs;
        if strcmpi(STACK{1}{1}.name,'load_stim_resps_wehr'),
            for setidx=1:length(XXX{end}.test_set),
                ts=XXX{end}.test_set{setidx};
                dat=XXX{end}.dat.(ts);
                ss=find(diff(isnan(dat.respavg(:)))==-1);
                dd=find(diff(isnan(dat.respavg(:)))==1);
                maxlen=max(dd-ss);
                pp=zeros(maxlen,length(dd));
                rr=zeros(maxlen,length(dd));
                for repidx=1:length(dd),
                    pp(:,repidx)=dat.stim((-maxlen+1:0)+dd(repidx));
                    rr(:,repidx)=dat.respavg((-maxlen+1:0)+dd(repidx));
                end
                if setidx==1,
                    os=nanmax(rr(:))/10;
                end
                rr=rr+os;
                pp=pp+os;
                if setidx==2,
                    rr=-rr;
                    pp=-pp;
                end
                tt=(1:size(rr,1))./fs;
                if setidx==1,
                    plot(tt,nanmean(rr,2),'r');
                else
                    plot(tt,nanmean(rr,2),'b');
                end
                hold on
                plot(tt,nanmean(pp,2),'k');
            end
            axis tight;
            hold off
            if (midx==1),
                xlabel('time (sec)');
                ylabel('<--I E-->');
            end
        else
            % standard baphy data, plot val
            ts=XXX{end}.test_set{1};
            dat=XXX{end}.dat.(ts);
            rr=dat.respavg(:,1);
            pp=dat.stim(:,1);
            
            plen=min(round(STACK{1}{1}.raw_resp_fs.*3.5),...
                     size(rr,1));
            if SPECIAL_FORMAT,
                rr=rr(1:plen,:);
                pp=pp(1:plen,:);
            end
            
            tt=(1:size(rr,1))./fs;
            plot(tt,nanmean(rr,2),'r','LineWidth',2);
            hold on
            plot(tt,nanmean(pp,2),'k','LineWidth',2);
            axis tight;
            hold off
            if (midx==1),
                ylabel('resp');
            end
            
            % plot stim in bottom row
            if PLOT_STIM,
                subplot(maxplots,modelcount,((maxplots-1).*modelcount+midx));
                ts=XXX{2}.test_set{1};
                dat=XXX{2}.dat.(ts);
                sstim=squeeze(dat.stim(:,1,:));
                sstim=sstim(1:plen,:);
                plot(tt,sstim,'LineWidth',2);
                if (midx==1),
                    xlabel('time (sec)');
                    ylabel('stim');
                end
            end

        end
    end
    
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