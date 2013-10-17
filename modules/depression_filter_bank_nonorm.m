function m = depression_filter_bank_nonorm(args)
% apply bank of synaptic depression filters to each input channel

    global NARF_DEBUG
    
% Module fields that must ALWAYS be defined
m = [];
m.mdl = @depression_filter_bank_nonorm;
m.name = 'depression_filter_bank_nonorm';
m.fn = @do_depression_filter;
m.pretty_name = 'Depression Filter Bank No Norm';
m.editable_fields = {'num_channels', 'strength', 'tau', 'gain',...
                    'per_channel', 'input', 'time', 'output' };
m.isready_pred = @isready_always;

% Module fields that are specific to THIS MODULE
m.num_channels = 2;     % is this a useful parameter?
m.per_channel = 0;
m.strength = [0 0.2];   % as fraction of stimulus max magnitude
m.tau = [0 100];        % in ms
m.gain = [1 1];        % in ms
m.strength_norm = 100;  % should be approximately stimmax
m.tau_norm = 1000;      % 1000 means ms, 10 means 10ths of sec, 1 means sec
m.raw_stim_freq = 100000;
m.input = 'stim';
m.time = 'stim_time';
m.output = 'stim';

% Optional fields
m.is_splittable = true;
m.auto_plot = @do_depression_cartoon;
m.auto_init = @auto_init_depression_filter_bank;

m.plot_gui_create_fn = @create_chan_selector_gui;
m.plot_fns = {};
m.plot_fns{1}.fn = @do_plot_all_default_outputs;
m.plot_fns{1}.pretty_name = 'Filtered Stimuli (All)';
m.plot_fns{2}.fn = @do_plot_single_default_output;
m.plot_fns{2}.pretty_name = 'Filtered Stimuli (Single)';
m.plot_fns{3}.fn = @do_plot_channels_as_heatmap;
m.plot_fns{3}.pretty_name = 'Filtered Stimuli (Heatmap)';
m.plot_fns{4}.fn = @do_depression_cartoon;
m.plot_fns{4}.pretty_name = 'Depression cartoon';

% Overwrite the default module fields with arguments 
if nargin > 0
    m = merge_structs(m, args);
end

% ------------------------------------------------------------------------
% INSTANCE METHODS

function mdl = auto_init_depression_filter_bank(stack, xxx)
    % NOTE: Unlike most plot functions, auto_init functions get a 
    % STACK and XXX which do not yet have this module or module's data
    % added to them. 
        disp('depression_filter_bank init');
        mdl=m;
        
        sf=fieldnames(xxx{end}.dat);
        sf=sf{1};
        [T, S, N] = size(xxx{end}.dat.(sf).(mdl.input));
        if mdl.per_channel,
            mdl.num_channels=length(mdl.tau)./N;
        else
            mdl.num_channels=length(mdl.tau);
        end
    end

% Finally, define the 'methods' of this module, as if it were a class
function x = do_depression_filter(mdl, x, stack, xxx)
    
    if ~isfield(mdl,'tau_norm'),
        mdl.tau_norm=1000;
    end
    if ~isfield(mdl,'strength_norm'),
        mdl.strength_norm=100;
    end
    
    rawfs=stack{1}{1}.raw_stim_fs;
    
    % bookkeeping for per-channel depression
    sf=fieldnames(x.dat);
    sf=sf{1};
    [T, S, N] = size(x.dat.(sf).(mdl.input));
    if isfield(mdl,'per_channel') && mdl.per_channel,
        num_channels=length(mdl.tau)./N;
    else
        num_channels=length(mdl.tau);
    end
    gainvector=abs(repmat(mdl.gain(1:num_channels),[N 1]));
    gainvector(gainvector<1)=exp(gainvector(gainvector<1)-1);
    gainvector(gainvector>1)=sqrt(gainvector(gainvector>1));
    gainvector=gainvector(:);
    
    % Now actually compute depression
    for sf = fieldnames(x.dat)', sf = sf{1};
        [T, S, N] = size(x.dat.(sf).(mdl.input));
        if isfield(mdl,'per_channel') && mdl.per_channel && num_channels<1,
            error('channel count does not work with per_channel==1');
        end
        ret = zeros(T, S, N*num_channels);
        for s = 1:S
            stim_in=squeeze(x.dat.(sf).(mdl.input)(:,s,:))';
            % require positive threshold
            stim_in=stim_in.*(stim_in>0);
            
            if isfield(mdl,'per_channel') && mdl.per_channel
                depresp=[];
                for jj=1:N,
                    tresp=depression_bank(...
                        stim_in(jj,:), mdl.strength(jj)./mdl.strength_norm,...
                        mdl.tau(jj) .* rawfs / mdl.tau_norm, 1);
                    depresp=cat(1,depresp,tresp);
                end
            else
               depresp=depression_bank(...
                   stim_in, mdl.strength ./ mdl.strength_norm,...
                   mdl.tau .* rawfs ./ mdl.tau_norm, 1);
            end
            if any(gainvector~=1),
                depresp=depresp .* ...
                    repmat(gainvector,[1 size(depresp,2) size(depresp,3)]);
            end
            depresp=permute(depresp',[1 3 2]);
            ret(:,s,:) = depresp;
        end
        
        x.dat.(sf).(mdl.output) = ret; 
    end
end

% % Plot the filter responses
function do_depression_cartoon(sel, stack, xxx)
    
    mdl = stack{end}{1};
    mdls = stack{end};
    fs=stack{1}{1}.raw_stim_fs;
    tau=[];
    strength=[];
    
    x = struct();
    x.dat.demo.stim=[zeros(fs./2,1);ones(fs,1)./2;zeros(fs,1);
                     ones(round(fs./10),1);zeros(round(fs./10),1);
                     ones(round(fs./10),1);zeros(7.*fs./10,1)];
    if mdl.per_channel,
        x.dat.demo.stim=repmat(x.dat.demo.stim,[1 1 length(mdl.tau)]);
    end
    x.dat.demo.(mdl.input)=x.dat.demo.stim;
    
    data=[];
    for jj=1:length(mdls),
        xfiltered=do_depression_filter(mdls{jj}, x, stack, xxx);
        xfiltered.dat.demo.stim=xfiltered.dat.demo.(mdl.input);
        
        data=cat(1,data,...
                 [x.dat.demo.stim(:,1)+1 squeeze(xfiltered.dat.demo.stim)]);
        data((end-9):end,:,:)=nan;
    end
    timeaxis=(1:size(data,1))'./fs;
    plot(timeaxis,data);
    
    axis([timeaxis([1 end])' -0.1 2.1]);
    for jj=1:length(mdls),
        text(1.5+(jj-1)*3.5,2,'stim','VerticalAlign','top');
        legstr={};
        for ii=1:length(mdls{jj}.tau);
            legstr{ii}=sprintf('%d:u=%.1f,tau=%.3f',...
                               ii,mdls{jj}.strength(ii),...
                               mdls{jj}.tau(ii)./mdls{jj}.tau_norm);
        end
        text(1.5+(jj-1)*3.5,1,legstr,'VerticalAlign','top');
    end
    
end
 
% function do_plot_filtered_stim(stack, xxx)
%     mdl = stack{end};
%     x = xxx{end};
% 
%     % Find the GUI controls
%     [sf, stim_idx, ~] = get_baphy_plot_controls(stack);
%     chan_idx = get(load_mod.plot_gui.selected_stim_chan_popup, 'Value');
%     dat = x.dat.(sf);
%     stepsize=load_mod.stimulus_channel_count;
%     plot(dat.(mdl.time), ...
%          squeeze(dat.(mdl.output)(:,stim_idx,chan_idx:stepsize:end)));
%     
%     axis tight;
% end
% 
% % Plot heatmap of the filter responses
% function do_plot_filtered_stim_heatmap(stack, xxx)
%     mdl = stack{end};
%     x = xxx{end};
%     
%    % Find the GUI controls   
%     [sf, stim_idx, ~] = get_baphy_plot_controls(stack);
%     dat = x.dat.(sf);
%     chancount=size(dat.(mdl.output),4);
%     imagesc(dat.(mdl.time),1:chancount, ...
%          squeeze(dat.(mdl.output)(:,stim_idx,:))');
%     
%     axis tight;
%     axis xy
% end

end

