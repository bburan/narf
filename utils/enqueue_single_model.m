function enqueue_single_model(batch, modulekeys, fitterkeys, cellid,...
                              training_set, test_set, force_rerun)
% function enqueue_single_model(batch, modulekeys, fitterkeys, cellid,...
%                             training_set, test_set, force_rerun)
%
% Enqueues a single model onto the distributed job queueing system so that
% a model can be processed on a bunch of machines.
%
% ARGUMENTS:
%    batch  
%    modulekeys     {'env100', 'log2', 'etc'}
%    fitterkeys     {'prior34', 'lsq'}
%    cellid
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

s = ['fit_single_model(' repl_write(modulekeys) ...
    ', ' repl_write(modulekeys) ...
    ', ' num2str(batch) ...
    ', ' repl_write(cellid) ...
    ', ' repl_write(training_set) ...
    ', ' repl_write(test_set) ...
    ')'];

% Build the modelname
tmp = cellfun(@(n) sprintf('%s_', n), modulekeys, 'UniformOutput', false);
modelname = strcat(tmp{:});
modelname = modelname(1:end-1); % Remove trailing underscore

% Build the fittername
tmp = cellfun(@(n) sprintf('%s_', n), fitterkeys, 'UniformOutput', false);
fittername = strcat(tmp{:});
fittername = fittername(1:end-1);

note = sprintf('%s/%d/%s/%s',cellid,batch,modelname, fittername);

sql = ['SELECT * FROM NarfResults WHERE cellid="',cellid,'"',...
    ' AND batch=',num2str(batch),...
    ' AND modelname="',modelname,'"', ...
    ' AND fittername="',fittername,'"'];
rdata = mysql(sql);

if ~isempty(rdata) && ~force_rerun,
    fprintf('NarfResults for %s already exist, skipping\n', note);
    return
end

sql=['SELECT * FROM tQueue WHERE parmstring="', s, '"'];
qdata=mysql(sql);

if ~isempty(qdata) && qdata(1).complete<=0
        fprintf('Incomplete queue entry for %s already exists, skipping.\n',...
                note);
        if ~strcmp(note,qdata(1).note),
            disp('fixing truncated tQueue.note');
            sql=['UPDATE tQueue set note="',note,'"',...
                 ' WHERE id=',num2str(qdata(1).id)];
            mysql(sql);
        end
    return
elseif ~isempty(qdata) && qdata(1).complete==2,
    fprintf('Dead queue entry for %s already exists, skipping.\n', note);
    return
elseif ~isempty(qdata) && qdata(1).complete==1,
    fprintf('Resetting existing queue entry for %s\n', note);
    dbsetqueue(qdata(1).id,0,0);
else
    fprintf('Adding new queue entry for %s\n', note);
    dbaddqueuemaster(s,note);
end
