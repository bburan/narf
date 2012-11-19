function params = gammatone_filterbank(args)

% Default Parameters
params = [];
params.bank_min_freq = 100;
params.bank_max_freq = 30000;
params.num_gammatone_filters = 32;
params.t_win = 0.025;              % Gamma filter window size [seconds]
params.t_hop = 1/resp_fs;          % Time between gamma window centers [seconds]
params.stim_sample_freq = 100000;
params.align_phase = false;        % Align the gamma filters to be in phase?
params.pretty_name = 'Gammatone Filter Bank';
params.editable_fields = {'bank_min_freq', 'bank_max_freq', ...
                          'num_gammatone_filters', 'align_phase', 'sampfs'};
           
% Overwrite the defaults with the arguments iff they are in editable_fields
if nargin == 1
    fns = fieldnames(args);
    for idx = 1:length(fns);
        if any(ismember(params.editable_fields, fns{idx}))
            params.(fns{idx}) = args.(fns{idx});
        end
    end
end

% Not needed anymore
% mktimeaxis = @(x, fs) ( [0:(1/fs):(length(x) - 1)/fs]);
% t_resp = mktimeaxis(resp,resp_fs);     % Response time axis

% -----------------------------------------------------------------------
% Values derived mathematically from the above. 


% ------------------------------------------------------------------------
% INSTANCE METHODS FOLLOW

function filtered_x = do_gammatone_filter(parm, x)
    
    [gamma_resp, gamma_envs, gamma_frqs] = gammatonebank(x, ...
                      parm.bank_min_freq, parm.bank_max_freq, ...
                      parm.num_gammatone_filters, parm.stim_sample_freq, ...
                      parm.alignphase);                  
    
                  
                  
    N_win = round(t_win * stim_fs);        % Number of samples in each window
    N_hop = round(t_hop * stim_fs);        % Number of samples to hop center by
    N_cols = 1 + floor((length(stim) - N_win) / N_hop); % Total number of windows

    % Rectify and smooth the gamma response for visualization
    gamma_sqr = [gamma_resp.^2];        
    gamma_pow = zeros(N_gfs, N_cols);
    for i = 1:N_cols
        gamma_pow(:,i) = 20*log10(sqrt(mean(gamma_sqr(:,(i-1)*N_hop + [1:N_win]),2)));
        %gamma_pow(:,i) = sqrt(mean(gamma_sqr(:,(i-1)*N_hop + [1:N_win]),2));
    end

    %              
    %filtered_x=[];
    %for i = 1:length(parm.low_freqs)
%        tmp = filter(parm.coefs{i}{1}, parm.coefs{i}{2}, x,[],2);
%        filtered_x = cat(3, filtered_x, tmp); 
%    end
end
 
function do_plot_gammatone_filter_as_colormap(parm)    
    subplot('Position', [0.05 0.6 0.9 0.18]);
    imagesc(gamma_pow); axis xy; 
    % caxis(clipOutliers(gamma_pow,1));
    % colorbar;
    setAxisLabelCallback(gca, @(x)(1000*t_hop*x), 'X');
    setAxisLabelCallback(gca, @(y)round(gamma_frqs(y)), 'Y');
    % set(gca,'XLim',[0 length(gamma_pow)]);
    xlabel('Time [ms]');
    ylabel('Frequency [Hz]');
    if parm = align_phase
        h=title('Gammatone Filter Response (Phase locked filters)');
    else
        h=title('Gammatone Filter Response');
    end
    
end

% ------------------------------------------------------------------------
% Put the instance methods in the params struct
params.freq_resp_plot_fn = @do_plot_elliptic_bandpass_filter_bank_frq_resp;
params.preproc_fn = @do_elliptic_filter;

end
