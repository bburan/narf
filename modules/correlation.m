function m = correlation(args)
% A function to linearly correlate two signals. 
% 
% The signals are assumed to be two dimensional matrices:
%  Dimension 1 is time
%  Dimension 2 is the stimuli #
%
% Returns a function module 'm' which implements the MODULE interface.
% See docs/modules.org for more information.

% Module fields that must ALWAYS be defined
m = [];
m.mdl = @correlation;
m.name = 'correlation';
m.fn = @do_correlation;
m.pretty_name = 'Correlation';
m.editable_fields = {'input1', 'input2', 'time', 'train_score', 'test_score'};
m.isready_pred = @isready_always;

% Module fields that are specific to THIS MODULE
m.input1 = 'stim';
m.input2 = 'respavg';
m.input_resp = 'resp';
m.time   = 'stim_time';
m.train_score = 'score_train_corr';
m.test_score = 'score_test_corr';
m.test_r_ceiling = 'score_test_ceilingcorr';
m.test_r_floor = 'score_test_floorcorr';
m.output = 'corr_score';  

% Overwrite the default module fields with arguments 
if nargin > 0
    m = merge_structs(m, args);
end

% Optional fields
m.plot_fns = {};
m.plot_fns{1}.fn = @do_plot_correlation_inputs;
m.plot_fns{1}.pretty_name = 'Correlation';

function x = do_correlation(mdl, x, stack, xxx)            
    % Compute the training set correlation, ignoring nans
        
    p = flatten_field(x.dat, x.training_set, mdl.input1);
    q = flatten_field(x.dat, x.training_set, mdl.input2); 
    
    R = corrcoef(excise([p q]));
    if isnan(R)
        x.(mdl.train_score) = NaN; 
        x.(mdl.output) = -2;
    else
        x.(mdl.train_score) = R(2,1);
        x.(mdl.output) = 1 - R(2,1);
    end
    
    % Compute the test set correlation, ignoring nans
    ptest = flatten_field(x.dat, x.test_set, mdl.input1);
    qtest = flatten_field(x.dat, x.test_set, mdl.input2); 
    R = corrcoef(excise([ptest qtest]));
    if isnan(R)
        x.(mdl.test_score) = NaN;
    else
        x.(mdl.test_score) = R(2,1);
    end
            
    if ~isfield(mdl, 'test_r_ceiling')
        mdl.test_r_ceiling = 'score_test_ceilingcorr';
    end
    if ~isfield(mdl, 'test_r_floor')
        mdl.test_r_floor = 'score_test_floorcorr';
    end    
    
    % TODO: 
    %  Come back and make this more general so that we compute respavg from
    %  resp and this module goes back to having just two inputs. (a vector
    %  and a scalar)
    
    % Extract resps as a matrix, not a vector    
    %if isfield(x.dat., 'resp')
        
        respmat = [];
        for ii = 1:length(x.test_set),
            [A, B, C] = size(x.dat.(x.test_set{ii}).resp);
            M = reshape(x.dat.(x.test_set{ii}).resp, A*B, C);
            if ~isempty(respmat) && C > size(respmat, 2)
                respmat = cat(1, [respmat nan(size(respmat,1), M-size(respmat,2))] , M);
            elseif C < size(respmat, 2)
                respmat = cat(1, respmat, [M nan(size(M,1), size(respmat,2) - M)]);
            else
                respmat = cat(1, respmat, M);
            end
        end
        
        % Remove rows with all NaNs
        idxs = find(~isnan(qtest) & ~isnan(ptest));                
        tmp = [ptest qtest respmat];
        tmp = tmp(idxs, :);
        tmp = excise(tmp')';
        
        pe = tmp(:,1);
        qe = tmp(:,2);
        re = tmp(:,3:end);
        
        x.(mdl.test_r_ceiling) = rceiling(pe, re);
        [~, x.(mdl.test_r_floor)] = rfloor(pe, qe);
        
%     else
%         tmp = excise([ptest qtest]);
%         pe = tmp(:,1);
%         qe = tmp(:,2);        
%         
%         x.(mdl.test_r_ceiling) = NaN;
%         [~, x.(mdl.test_r_floor)] = rfloor(pe, qe);
%     end
%     
end

function do_plot_correlation_inputs(sel, stack, xxx)
    [mdls, xins, xouts] = calc_paramsets(stack, xxx(1:end-1)); 
    xout = xxx{end};
    mdl = mdls{1};
    
    do_plot_scatter(sel, {xxx}, mdls{1}.input1, mdls{1}.input2, 500);
    
    if ~isfield(mdl, 'test_r_ceiling')
        mdl.test_r_ceiling = 'score_test_ceilingcorr';
    end
    if ~isfield(mdl, 'test_r_floor')
        mdl.test_r_floor = 'score_test_floorcorr';
    end    
    textLoc(sprintf(' Train r: %f\n Test r : %f ({\\pm}%f)\n Ceiling r: %f', ...
                    xout.(mdl.train_score), xout.(mdl.test_score),...
                    xout.(mdl.test_r_floor), xout.(mdl.test_r_ceiling)),...
            'NorthWest');
    
    % If 
    
    
    do_xlabel('Prediction [Hz]');
    do_ylabel('Actual Response [Hz]');
    
end


end