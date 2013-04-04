function batch_241_table_maker()

batch = 241;
cellids = {'por018b-c1', ...
          'por018c-c1', ...
          'por019b-a1', ...
          'por022a-a1', ...
          'por022a-c1', ...
          'por023b-b1', ...
          'por027a-a1', ...
          'por027b-b2', ...
          'por027b-c1'};
    
% Plot a heat map of stuff
function fh = phm(M, xl, yl, the_title) 
    fh = figure('Position', [0 0 2000 1200]); clf;
    heatmap(M, xl, yl, '%2.0f', 'TickAngle', 90,...
                'ShowAllTicks', true, 'TickFontSize', 12);
    set(gca,'Position',[.05 .3 .9 .65])
    title(the_title);
end
      
m_all_cellarray = {};
    
for ii = 1:length(cellids),
    cellid = cellids{ii};
    
    m_all = db_get_models(batch, cellid, {'npnlx'});
        
    % Append to our global cell array
    for jj = 1:length(m_all)
        m_all_cellarray{end+1} = m_all(jj);
        m_all_cellarray{end}.modelname = char(m_all_cellarray{end}.modelname);
        m_all_cellarray{end}.modelpath = char(m_all_cellarray{end}.modelpath);
        m_all_cellarray{end}.figurefile = char(m_all_cellarray{end}.figurefile);
    end
      
%    M = [];
%     
%     % Group models by second tokens
%     for jj = 1:length(m_all)
%         m = m_all(jj);
%         toks = tokenize_string(char(m.modelname));
%         tok = toks{2}{1};
%         if isfield(M, tok)
%             M.(tok) = cat(1, M.(tok), m.r_test);
%         else
%             M.(tok) = [m.r_test];
%         end
%     end
%           
%     % Print a table entry showing the average of each grouping
%     fns = sort(fieldnames(M));
%     fprintf('%s\n', cellid);
%     for jj = 1:length(fns)
%         f = fns{jj};
%         fprintf('\t%s: %.2f<%.2f<%.2f\n', f, min(M.(f)), nanmean(M.(f)), max(M.(f)));
%         %[Y, X] = hist(M.(f));
%     end
    
    
end  

[M, xl, yl] = select_summaries(m_all_cellarray, ...
    @(c) 100*getfield(c, 'r_test'), ...
    @(cc) rotate_tokens(cc{1}.modelname, 2));
fh = phm(M, xl, yl, 'Test Set R, Sorted by 2nd token');
                                

end
