function plot_examples_for_paper(batch, cellids, modelnames)
global STACK META XXX;

modelnames = {'fb18ch100_lognn_fir15_siglog100_fit05h_fit05c'...
    'fb18ch100_lognn_wc01_fir15_siglog100_fit05h_fit05c', ...
    'fb18ch100_lognn_wcg01_fir15_siglog100_fit05h_fit05c', ...
    'fb18ch100_lognn_wcg01_ap3z1_dexp_fit09c', ...
    'fb18ch100_lognn_wc02_fir15_siglog100_fit05h_fit05c', ...
    'fb18ch100_lognn_wcg02_fir15_siglog100_fit05h_fit05c', ...
    'fb18ch100_lognn_wcg02_ap3z1_dexp_fit09c', ...
    'fb18ch100_lognn_wc03_fir15_siglog100_fit05h_fit05c', ...
    'fb18ch100_lognn_wcg03_fir15_siglog100_fit05h_fit05c', ...
    'fb18ch100_lognn_wcg03_ap3z1_dexp_fit09c', ...
    'fb18ch100_lognn_wc04_fir15_siglog100_fit05h_fit05c', ...
    'fb18ch100_lognn_wcg04_fir15_siglog100_fit05h_fit05c', ...
    'fb18ch100_lognn_wcg04_ap3z1_dexp_fit09c'};

n_tall = 5;
n_wide = 3;
subplot = @(m,n,p) subtightplot (m, n, p, [0.01 0.05], [0.1 0.01], [0.1 0.01]);

for ii = 1:length(cellids)
    
    % Model path
    mp = ['/auto/data/code/saved_models/' num2str(batch) '/' ...
        cellids{ii} '/' num2str(batch) '_' cellids{ii} '_'] ;
  
    figure('Name', cellids{ii}, ...
           'NumberTitle', 'off', 'Position', [20 50 650 1000]);       
   
	% Plot TORC-estimated STRF
    subplot(n_tall,n_wide,3);    
    [cfd, ~, ~] = dbgetscellfile('cellid', cellids{ii});
    for jj = 1:length(cfd);
        % TODO: Replace magic number  1 with better description of TOR files
        if (cfd(jj).runclassid == 1)
            strf_offline2([cfd(jj).stimpath cfd(jj).stimfile], ...
                [cfd(jj).path cfd(jj).respfile], ...
                cfd(jj).channum, cfd(jj).unit);
        end
    end
    xlabel(''); set(gca,'XtickLabel',[]);
    
    % Plot full-channel STRF
    subplot(n_tall,n_wide,1);
    load_model([mp modelnames{1} '.mat']);
    calc_xxx(1);
    [mods, mods_idx] = find_modules(STACK, 'fir_filter');
    m = mods{end};
    idx = mods_idx{end};
    m = m{1};
    m.auto_plot({}, STACK(1:idx), XXX(1:idx));
    textLoc(sprintf('r_test=%0.4f', XXX{end}.score_test_corr),'East', 'Interpreter', 'none' );
    ylabel(''); xlabel(''); set(gca,'XtickLabel',[],'YtickLabel',[]);
        
    for jj = 2:length(modelnames)
        load_model([mp modelnames{jj} '.mat']);
        calc_xxx(1);
        
        [m, idx] = find_modules(STACK, 'pole_zeros', 1);
        
        if isempty(m)
            [mods, mods_idx] = find_modules(STACK, 'fir_filter');
            if isempty(mods)
                continue;
            end
            m = mods{end};
            idx = mods_idx{end};
        end
        
        m = m{1};
        subplot(n_tall, n_wide, jj+2)
        m.auto_plot({}, STACK(1:idx), XXX(1:idx));
        textLoc(sprintf('r_test=%0.4f', XXX{end}.score_test_corr),'East', 'Interpreter', 'none');
        ylabel('');
        xlabel('');
        set(gca,'XtickLabel',[],'YtickLabel',[]);
    end
end

end