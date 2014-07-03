function switch_fittable(choice)
global STACK;

% This is a helper function designed to be called with 'all', 'default' or
% any string. This function will add/remove variables in fit_fields
% depending on what the 'fitter' field of the 'fit_constraints'.


if strcmp(choice, 'all')
    % we allow the optimization of all fields
    for ii = 1:length(STACK)
        if isfield(STACK{ii}{1}, 'fit_fields') && ...
                isfield(STACK{ii}{1}, 'fit_constraints')
            for i = 1:length(STACK{ii}{1}.fit_constraints)
                field = STACK{ii}{1}.fit_constraints{i}.var;
                % turn on any variable (if required)
                matching_field = strcmp(field, STACK{ii}{1}.fit_fields);
                if sum(matching_field) == 0
                    % no match => turn on
                    STACK{ii}{1}.fit_fields = {STACK{ii}{1}.fit_fields{:} field};
                end
            end
        end
    end
elseif strcmp(choice, 'default')
    % we allow only the optimization of fields without any
    % alternative fitter
    for ii = 1:length(STACK)
        if isfield(STACK{ii}{1}, 'fit_fields') && ...
                isfield(STACK{ii}{1}, 'fit_constraints')
            for i = 1:length(STACK{ii}{1}.fit_constraints)
                field = STACK{ii}{1}.fit_constraints{i}.var;
                if isfield(STACK{ii}{1}.fit_constraints{i}, 'fitter')
                    % there is a fitter indicated => turn off
                    % only if matching variable
                    matching_field = strcmp(field, STACK{ii}{1}.fit_fields);
                    % turn the variable off if necessary
                    if sum(matching_field) > 0
                        % match => turn off
                        STACK{ii}{1}.fit_fields = {STACK{ii}{1}.fit_fields{~matching_field}};
                    end
                else
                    % there is no fitter indicated => turn on
                    % any variable
                    matching_field = strcmp(field, STACK{ii}{1}.fit_fields);
                    if sum(matching_field) == 0
                        % no match => turn on
                        STACK{ii}{1}.fit_fields = {STACK{ii}{1}.fit_fields{:} field};
                    end
                end
            end
        end
    end
else
    % we allow only the optimization of fields with the specified
    % alternative fitter
    for ii = 1:length(STACK)
        if isfield(STACK{ii}{1}, 'fit_fields') && ...
                isfield(STACK{ii}{1}, 'fit_constraints')
            for i = 1:length(STACK{ii}{1}.fit_constraints)
                field = STACK{ii}{1}.fit_constraints{i}.var;
                if isfield(STACK{ii}{1}.fit_constraints{i}, 'fitter')
                    % there is a fitter indicated => turn on
                    % only if matching variable
                    matching_field = strcmp(field, STACK{ii}{1}.fit_fields);
                    if strcmp(choice, STACK{ii}{1}.fit_constraints{i}.fitter)
                        % turn the variable on if necessary
                        if sum(matching_field) == 0
                            % no match => turn on
                            STACK{ii}{1}.fit_fields = {STACK{ii}{1}.fit_fields{:} field};
                        end
                    else
                        % turn the variable off
                        if sum(matching_field) > 0
                            STACK{ii}{1}.fit_fields = {STACK{ii}{1}.fit_fields{~matching_field}};
                        end
                    end
                else
                    % there is no fitter indicated => turn off
                    % any variable
                    matching_field = strcmp(field, STACK{ii}{1}.fit_fields);
                    if sum(matching_field) > 0
                        % match & on => turn off
                        STACK{ii}{1}.fit_fields = {STACK{ii}{1}.fit_fields{~matching_field}};
                    end
                end
            end
        end
    end
end
end

