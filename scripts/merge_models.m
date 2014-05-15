function fitset=merge_models(batch, cellids, modelnames)

global XXX STACK META


% this version makes the usual scaling of weights so the sum is equal to 1
    function x = format_coef(x)
        x = abs(x);
        x = x/sum(x);
    end

% this version constrains the coef to be higher than (roughly) 1/2n
%     function x = format_coef(x)
%         x = abs(x);
%         x = x/sum(x);
%         x = x + 1/length(x);
%         x = x/sum(x);
%     end

    function score = getfit(x, exp_train, mdl_train)
        predic = zeros(size(exp_train));
        x = format_coef(x);
        for m=1:size(mdl_train,2),
            predic = predic + mdl_train(:,m)*x(m);
        end
        score = nanmean((predic - exp_train).^2);
        score = score / (nanvar(exp_train)+(score==0));
    end


if nargout==0,
    [FileName,PathName,~] = ...
        uiputfile(['model_merge_', num2str(batch), '.csv'],'Save to...');
    if ~FileName,
        disp('Canceled save.');
        return;
    else
        fprintf('Saving csv to %s%s\n',PathName,FileName);
        fid = fopen([PathName FileName],'w');
        fprintf(fid,'CELL_ID MODEL TRAINING_SCORE TESTING_SCORE TRAINING_CORR TESTING_CORR COEF\n');
    end
end


stack_extractor = @(x) x;
meta_extractor = @(x) x;
[stacks, metas, x0s,preds] = load_model_batch(batch, cellids, modelnames, ...
    stack_extractor, meta_extractor);
cellcount=length(cellids);
modelcount=length(modelnames);


optimization_set = 'training_set';
testing_set = 'test_set';
% optimization_set = 'test_set';
% testing_set = 'training_set';


for cc=1:cellcount,
    fprintf('********\nDoing cell %d\n********\n',cc);
    
    % STEP 1:
    % we extract the original recording along with the model predictions
    train_nmse = zeros(modelcount,1);
    train_cor = zeros(modelcount,1);
    for mm=1:modelcount,
        fprintf('********\nDoing model %d\n********\n',mm);
        
        XXX={x0s{mm,cc}};
        STACK=stacks{mm,cc};
        META=metas{mm,cc};
        calc_xxx(1);
        xxx = XXX{end};
        if mm==1,
            % we get the original stimulus only in the query for the first model
            exp_train = flatten_field(xxx.dat, xxx.(optimization_set), 'respavg');
            mdl_train = zeros(length(flatten_field(xxx.dat, xxx.(optimization_set), 'stim')), modelcount);
     %%%%%       merge_train = zeros(length(flatten_field(xxx.dat, xxx.(optimization_set), 'stim')), 1);
        end
        
        if (size(mdl_train(:,mm)) == size(flatten_field(xxx.dat, xxx.(optimization_set), 'stim'))),
            mdl_train(:,mm) = flatten_field(xxx.dat, xxx.(optimization_set), 'stim');
            train_score = nanmean((mdl_train(:,mm) - exp_train).^2);
            train_nmse(mm) = train_score / (nanvar(exp_train)+(train_score==0));
            R = corrcoef(excise([mdl_train(:,mm) exp_train]));
            train_cor(mm) = R(2,1);
        else
            fprintf('\n\n#######################\nBUG in model %s !!!\n#######################\n\n', META.modelname);
        end
    end
    
    
    % STEP 2
    % we optimize the weights so to maximize the merged nmse in training
    test_nmse = zeros(modelcount,1);
    test_cor = zeros(modelcount,1);
    
    % the following gets a good guess for the start and optimize from there
%     starts = eye(modelcount);
%     best_start = format_coef(repmat(1/modelcount,modelcount,1));
%     best_perf = getfit(best_start, exp_train, mdl_train);
%     for start=1:modelcount,
%         s = format_coef(starts(start,:))
%         perf = getfit(s, exp_train, mdl_train);
%         if (perf < best_perf)
%             best_start = s;
%             best_perf = perf;
%         end
%     end
%     coefs = fminsearch(@(x) getfit(x, exp_train, mdl_train), best_start);
%     coefs = format_coef(coefs);
    
    % the following optimizes from different good starting guesses
    % and keeps the overall best solution
    starts = eye(modelcount);
    best_coefs = fminsearch(@(x) getfit(x, exp_train, mdl_train), repmat(1/modelcount,modelcount,1));
    best_coefs = format_coef(best_coefs);
    best_perf = getfit(best_coefs, exp_train, mdl_train);
    for start=1:modelcount,
        coefs = fminsearch(@(x) getfit(x, exp_train, mdl_train), starts(start,:));
        coefs = format_coef(coefs);
        perf = getfit(coefs, exp_train, mdl_train);
        if (perf < best_perf)
            best_coefs = coefs;
            best_perf = perf;
        end
    end
    coefs = best_coefs;
    
    
        % the following gets a good guess for the start and optimize from there
%     coefs = fminsearch(@(x) getfit(x, exp_train, mdl_train), repmat(1/modelcount,modelcount,1));
%     coefs = format_coef(coefs);

    for mm=1:modelcount,
        XXX={x0s{mm,cc}};
        STACK=stacks{mm,cc};
        META=metas{mm,cc};
        calc_xxx(1);
        xxx = XXX{end};
        if mm==1
            % we get the original stimulus only in the query for the first model
            exp_test = flatten_field(xxx.dat, xxx.(testing_set), 'respavg');
            mdl_test = zeros(length(flatten_field(xxx.dat, xxx.(testing_set), 'stim')), modelcount);
            merge_test = zeros(length(flatten_field(xxx.dat, xxx.(testing_set), 'stim')), 1);
            merge_train = zeros(length(flatten_field(xxx.dat, xxx.(optimization_set), 'stim')), 1);
        end
        mdl_test(:,mm) = flatten_field(xxx.dat, xxx.(testing_set), 'stim');
        test_score = nanmean((mdl_test(:,mm) - exp_test).^2);
        test_nmse(mm) = test_score / (nanvar(exp_test)+(test_score==0));
        R = corrcoef(excise([mdl_test(:,mm) exp_test]));
        test_cor(mm) = R(2,1);
        
        merge_test = merge_test + coefs(mm) * mdl_test(:,mm);
        merge_train = merge_train + coefs(mm) * mdl_train(:,mm);
    end
    
    train_score = nanmean((merge_train - exp_train).^2);
    train_nmse2 = train_score / (nanvar(exp_train)+(train_score==0));
    R = corrcoef(excise([merge_train exp_train]));
    train_cor2 = R(2,1);
    
    test_score = nanmean((merge_test - exp_test).^2);
    test_nmse2 = test_score / (nanvar(exp_test)+(test_score==0));
    R = corrcoef(excise([merge_test exp_test]));
    test_cor2 = R(2,1);
    
    for mm = 1:modelcount
        fprintf(fid,'%s %s %f %f %f %f %f\n', cellids{cc} , modelnames{mm}, train_nmse(mm), test_nmse(mm), train_cor(mm), test_cor(mm), coefs(mm));
    end
    fprintf(fid,'%s merged_model %f %f %f %f 0\n', cellids{cc}, train_nmse2, test_nmse2, train_cor2, test_cor2);
    
end

fclose(fid);


end

