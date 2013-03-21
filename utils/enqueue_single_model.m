%function enqueue_single_model(modulekeys, batch, cellid, training_set, ...
%                              test_set, force_rerun[=0])
%
% if force_rerun -- add jobs to queue even if results already exist
% in NarfResults
%
function enqueue_single_model(modulekeys, batch, cellid, training_set, ...
                              test_set, filecodes, force_rerun)
    
    if ~exist('force_rerun','var'),
        force_rerun=0;
    end
    
    s = ['fit_single_model(' repl_write(modulekeys) ...
         ', ' num2str(batch) ...     
         ', ' repl_write(cellid) ...
         ', ' repl_write(training_set) ...
         ', ' repl_write(test_set) ... 
         ', ' repl_write(filecodes) ... 
         ')'];
    
    tmp = cellfun(@(n) sprintf('%s_', n), modulekeys, 'UniformOutput', false);
    modelname = strcat(tmp{:});
    modelname = modelname(1:end-1); % Remove trailing underscore
            
    note=sprintf('%s/%d/%s',cellid,batch,modelname);
    
    sql=['SELECT * FROM NarfResults WHERE cellid="',cellid,'"',...
         ' AND batch=',num2str(batch),...
         ' AND modelname="',modelname,'"'];
    rdata=mysql(sql);
    
    if ~isempty(rdata) && ~force_rerun,
        fprintf('NarfResults for %s/%d/%s already exist, skipping\n',...
                cellid,batch,modelname);
        return
    end
    
    sql=['SELECT * FROM tQueue WHERE parmstring="',s,'"'];
    qdata=mysql(sql);
    
    if ~isempty(qdata) && qdata(1).complete<=0
        fprintf('Incomplete queue entry for %s/%d/%s already exists, skipping.\n',...
                cellid,batch,modelname);
        return
    elseif ~isempty(qdata) && qdata(1).complete==2,
        fprintf('Dead queue entry for %s/%d/%s already exists, skipping.\n',...
                cellid,batch,modelname);
        return
    elseif ~isempty(qdata) && qdata(1).complete==1,
        fprintf('Resetting existing queue entry for %s/%d/%s\n',...
                cellid,batch,modelname);
        dbsetqueue(qdata(1).id,0,0);
    else
        fprintf('Adding new queue entry for %s/%d/%s\n',...
                cellid,batch,modelname);
        %fprintf('%s', s);
        dbaddqueuemaster(s,note);
    end
end