% function celldb_wrapper(cellid,batch)
%
% wrapper for narf code that pulls train/test files from celldb
% based on batch data.
%
%
% created SVD 1/11/13 -- ripped off of cellxcqueue
%
function celldb_wrapper(cellid,batch)

close all
drawnow

disp('celldb_wrapper.m:');
baphy_set_path;
narf_set_path;
global NARF_PATH STACK XXX mdls;
mdls = scan_directory_for_modules();

if ~exist('batch','var'),
   disp('syntax error: spn_wrapper(cellid,batch) parameters required');
   return
end

dbopen;
sql=['SELECT * from sRunData WHERE cellid="',cellid,'"',...
     ' AND batch=',num2str(batch)];
rundata=mysql(sql);

if length(rundata)==0,
   fprintf('(cellid,batch)=(%s,%d) not found!\n',cellid,batch);
   return
end

fprintf('CELLID=%s BATCH=%d\n',cellid,batch);

% figure out what files to use for what stage of the analysis
[cellfiledata,times,params]=cellfiletimes(cellid,rundata.batch);

models = {'linear_fit_spn_depression'};

close all

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

global NARF_PATH STACK XXX;

savepath = [NARF_PATH filesep 'saved_models' filesep cellid];
if ~exist(savepath,'dir'),
    mkdir(fileparts(savepath),cellid);
end

% Build the model and train it 
linear_fit_spn_depression(cellid, train_set);

open_narf_gui;
hnarf=gcf;
drawnow;

STACK{4}.fit_fields={'coefs','baseline'};
recalc_xxx(3);
fit_with_lsqcurvefit;

close(hnarf);
open_narf_gui;
hnarf=gcf;
drawnow;

linpredmin=min(min(XXX{5}.dat.(train_set{1}).stim));
linpredmax=max(max(XXX{5}.dat.(train_set{1}).stim));

STACK{5} = mdls.nonparm_nonlinearity;
%phi3=max(max(XXX{5}.dat.(train_set{1}).respavg))./6;
%phi1=linpredmin+(linpredmax-linpredmin)./10;
%phi2=0.001;
%STACK{5} = mdls.nonlinearity.mdl(struct('phi', [phi1 phi2 phi3 0], ...
%                                        'nlfn', @sigmoidal));
%phi2=linpredmin+(linpredmax-linpredmin)./10;
%phi1=0.01;
%STACK{5} = mdls.nonlinearity.mdl(struct('phi', [phi1 phi2], ...
%                                        'nlfn', @exponential));
STACK{6} = mdls.correlation;

%STACK{4}.fit_fields = {};
%STACK{5}.fit_fields = {'phi'};
%fit_with_lsqcurvefit();
recalc_xxx(1); 

% Add the test set AFTER the training has occured, so that it is
% not needlessly computed by default by the system. Recompute so
% that the training and test scores are updated.
XXX{1}.test_set = test_set;
recalc_xxx(1); 

% Save the file name 
filename = sprintf('%s/%s__%f__%f__bat%d.mat', ...
                   savepath, cellid, ...
                   XXX{end}.score_train_corr, ...
                   XXX{end}.score_test_corr, ...
                   batch);

save_model_stack(filename, STACK, XXX);


% Plot the results
close(hnarf);
open_narf_gui;




