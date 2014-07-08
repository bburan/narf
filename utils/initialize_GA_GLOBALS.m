function initialize_GA_GLOBALS(phi_init)

global GA_MaskFittable GA_LowerBound GA_UpperBound GA_Bounded STACK;

% This helper function initializes all the global variables containing the
% constraints.
% !! caveat: it was written before the switch to flatXXX and flatSTACK.
% Everything should work fine as long as there is no split module.


% for all modules, we check the size of fit_constraints.lower and
% fit_constraints.upper and extend it to match fit_fields if
% relevant
for ii = 1:length(STACK)
    if isfield(STACK{ii}{1}, 'fit_fields') && ...
            isfield(STACK{ii}{1}, 'fit_constraints')
        for i = 1:length(STACK{ii}{1}.fit_fields)
            field = cell2mat(STACK{ii}{1}.fit_fields(i));
            [idx, con] = get_constraint_index(STACK{ii}{1}.fit_constraints, field);
            if idx
                for bound = {'lower', 'upper'}
                    bound = cell2mat(bound);
                    if isfield(con, bound) && ~isequal(size(STACK{ii}{1}.(field)), size(con.(bound)))
                        if isequal(size(con.(bound)), [1 1])
                            STACK{ii}{1}.fit_constraints{idx}.(bound) = ...
                                con.(bound)*ones(size(STACK{ii}{1}.(field), 1), size(STACK{ii}{1}.(field), 2));
                        else
                            error('Wrong constraint size')
                        end
                    end
                end
            end
        end
    end
end

% initialization of the MaskFittable and VarName
mask_rows = 0;
fittable_nb = 0;
for ii = 1:length(STACK)
    if isfield(STACK{ii}{1}, 'fit_fields')
        mask_rows = mask_rows + 1;
        for jj = 1:numel(STACK{ii}{1}.fit_fields)
            fittable_nb = fittable_nb + numel(STACK{ii}{1}.fit_fields(jj));
        end
    end
end

GA_MaskFittable = zeros(mask_rows, fittable_nb);
GA_LowerBound = -Inf(1, fittable_nb);
GA_UpperBound = Inf(1, fittable_nb);
GA_Bounded = false(1, fittable_nb);
mask_rowi = 1;
variable_index = 1;
for ii = 1:length(STACK)
    if isfield(STACK{ii}{1}, 'fit_fields')
        l = 0;
        for i = 1:length(STACK{ii}{1}.fit_fields)
            variable = cell2mat(STACK{ii}{1}.fit_fields(i));
            n_params_in_field = size(STACK{ii}{1}.(variable),1) * ...
                size(STACK{ii}{1}.(variable),2);
            if isfield(STACK{ii}{1}, 'fit_constraints')
                for j = 1:length(STACK{ii}{1}.fit_constraints)
                    if strcmp(STACK{ii}{1}.fit_constraints{j}.var, variable)
                        k_stack = 1;
                        for k = (variable_index+l):(variable_index+l+n_params_in_field-1)
                            if isfield(STACK{ii}{1}.fit_constraints{j}, 'lower')
                                GA_LowerBound(k) = STACK{ii}{1}.fit_constraints{j}.lower(k_stack);
                            end
                            if isfield(STACK{ii}{1}.fit_constraints{j}, 'upper')
                                GA_UpperBound(k) = STACK{ii}{1}.fit_constraints{j}.upper(k_stack);
                                if (GA_LowerBound(k) ~= -Inf)
                                    GA_Bounded(k) = true;
                                end
                            end
                            k_stack = k_stack + 1;
                        end
                    end
                end
            end
            l = l + n_params_in_field;
        end
        GA_MaskFittable(mask_rowi,variable_index:(variable_index+l-1)) = ii;
        variable_index = variable_index+l;
        mask_rowi = mask_rowi + 1;
    end
end

if sum(~GA_Bounded)
    fprintf('Warning: one or more module did not set lower+upper bounds.\n Partially unbounded optimization may be buggy but will proceed.\n');
end

end


