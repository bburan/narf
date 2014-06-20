
function m = mutual_information(args)
% A function to estimate the mutual information between two signals.
%
% Returns a function module 'm' which implements the MODULE interface.
% See docs/modules.org for more information.

% Module fields that must ALWAYS be defined
m = [];
m.mdl = @mutual_information;
m.name = 'mutual_information';
m.fn = @do_mutual_information;
m.pretty_name = 'Mutual information';
m.editable_fields = {'input1', 'input2', 'nbins', 'train_score', 'test_score'};
m.isready_pred = @isready_always;

% Module fields that are specific to THIS MODULE
m.input1 = 'stim';
m.input2 = 'respavg';
m.nbins  = 100; 
m.train_score = 'score_train_mi';
m.test_score = 'score_test_mi';
m.output = 'score_mi';  

% Overwrite the default module fields with arguments 
if nargin > 0
    m = merge_structs(m, args);
end

% Optimize this module for tree traversal  
m.required = {m.input1, m.input2};   % Signal dependencies
m.modifies = {m.output, m.train_score, m.test_score};  % These signals are modified
                      
% Optional fields
m.plot_fns = {};
m.plot_fns{1}.fn = @do_plot_mutual_information_inputs;
m.plot_fns{1}.pretty_name = 'Mutual Information';

function x = do_mutual_information(mdl, x)            

    % Compute the training set mutual_information, ignoring nans
    p = flatten_field(x.dat, x.training_set, mdl.input1);
    q = flatten_field(x.dat, x.training_set, mdl.input2);     
    vals = excise([p q]);
    [mi_train, n2train] = mi(vals(:,1), vals(:,2), mdl.nbins);
    x.(mdl.train_score) = mi_train; 
    x.trainvals = vals;
    x.trainn2 = n2train;
    
    % Compute the test set mutual_information, ignoring nans
    ptest = flatten_field(x.dat, x.test_set, mdl.input1);
    qtest = flatten_field(x.dat, x.test_set, mdl.input2); 
    vals = excise([ptest qtest]);
    [mi_test, n2test] = mi(vals(:,1), vals(:,2), mdl.nbins);
    x.(mdl.test_score) = mi_test;
    x.testvals = vals;
    x.testn2 = n2test;
end

function do_plot_mutual_information_inputs(sel, stack, xxx)
    [mdls, xins, xouts] = calc_paramsets(stack, xxx(1:end)); 
    xout = xxx{end};
    mdl = mdls{1};    
    
    imagesc(xout.testn2);
    textLoc(sprintf(' Est MI: %f\n Val MI : %f', ...
            xout.(mdls{1}.train_score), xout.(mdl.test_score)), 'NorthWest');

    do_xlabel('Bin');
    do_ylabel('Bin');
    
end


end
