function ret = request_celldb_batch(batch)
% Exploiting all the code written by Stephen, find the cellids, training,
% and test sets that correspond with a particular batch number. Returns a
% cell array of structs with three fields per struct:
%      .cellid
%      .training_set
%      .test_set

ret = {};

if ~exist('batch'),
   disp('syntax error: request_celldb_batch(batch) parameters required');
   return
end

dbopen;
sql = ['SELECT * from sRunData WHERE batch=', num2str(batch)];
rundata = mysql(sql);

if length(rundata) == 0,
   fprintf('batch number %d not found!\n', batch);
   return
end

models = {'linear_fit_spn_depression'};

for nn = 1:length(rundata)
    cellid = rundata(nn).cellid;
    
    fprintf('CELLID=%s BATCH=%d\n', cellid, batch);

    % figure out what files to use for what stage of the analysis
    [cellfiledata, times, params] = cellfiletimes(cellid, rundata(nn).batch);
    
    train_set={};
    test_set={};
    
    for ii=1:length(times(1).fileidx),
        if times(1).stop(ii)>times(1).start(ii),
            train_set{end+1}=basename(params.stimfiles{times(1).fileidx(ii)});
        end
    end
    
    for ii=1:length(times(3).fileidx),
        if times(3).stop(ii)>times(3).start(ii),
            test_set{end+1}=basename(params.stimfiles{times(3).fileidx(ii)});
        end
    end
    
    ret{nn} = [];
    ret{nn}.cellid = cellid;
    ret{nn}.training_set = train_set;
    ret{nn}.test_set = test_set;
    
end

