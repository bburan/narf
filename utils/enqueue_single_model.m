function enqueue_single_model(batch, cellid, modulekeys, ...
                              training_set, test_set, filecodes, force_rerun)
% function enqueue_single_model(batch, cellid, modulekeys, ...
%                             training_set, test_set, filecodes, force_rerun)
%
% Enqueues a single model onto the distributed job queueing system so that
% a model can be processed on a bunch of machines.
%
% ARGUMENTS:
%    batch  
%    cellid
%    modulekeys     {'env100', 'log2', 'etc'}
%    training_set
%    test_set
%    force_rerun    When true, jobs are added to queue even if a matching
%                   model is found to exist in NarfResults mysql table.
%
% RETURNS: nothing
%    

if ~exist('force_rerun','var'),
    force_rerun=0;
end

s = ['fit_single_model(' ...
    '' num2str(batch) ...
    ', ' write_readably(cellid) ...
    ', ' write_readably(modulekeys) ...
    ', ' write_readably(training_set) ...
    ', ' write_readably(test_set) ...
    ', ' write_readably(filecodes) ...
    ')'];

% Build the modelname
tmp = cellfun(@(n) sprintf('%s_', n), modulekeys, 'UniformOutput', false);
modelname = strcat(tmp{:});
modelname = modelname(1:end-1); % Remove trailing underscore

note = sprintf('%s/%d/%s',cellid,batch,modelname);

sql = ['SELECT * FROM NarfResults WHERE cellid="',cellid,'"',...
    ' AND batch=',num2str(batch),...
    ' AND modelname="',modelname,'"'];
rdata = mysql(sql);

if ~isempty(rdata) && ~force_rerun,
    fprintf('NarfResults for %s already exist, skipping\n', note);
    return
end

sql=['SELECT * FROM tQueue WHERE note="', note, '"'];
qdata=mysql(sql);

if ~isempty(qdata) && qdata(1).complete<=0
    fprintf('Incomplete entry for %s already exists, skipping.\n',note);
    if ~strcmp(note,qdata(1).note),
        disp('Fixing truncated tQueue.note');
        sql=['UPDATE tQueue set note="',note,'"',...
             ' WHERE id=',num2str(qdata(1).id)];
        mysql(sql);
    end
    
elseif ~isempty(qdata) && qdata(1).complete==2,
    fprintf('Dead queue entry for %s already exists. Resetting.\n', note);
    sql=['UPDATE tQueue SET complete=0, progress=0',...
         ' WHERE id=',num2str(qdata(1).id)];
    [res,r]=mysql(sql);
    
elseif ~isempty(qdata) && qdata(1).complete==1,
    fprintf('Resetting existing queue entry for %s\n', note);
    sql=['UPDATE tQueue SET complete=0, progress=0',...
         ' WHERE id=',num2str(qdata(1).id)];
    [res,r]=mysql(sql);
    
else
    fprintf('Adding new queue entry for %s\n', note);
    dbaddqueuemaster(s,note);
end
