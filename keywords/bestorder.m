function bestorder()
% JL - June 2015

% This script detects if there is a Lindeberg kernel in the STACK - if so,
% it tries fitting with all orders in ([0,1,2] x [0,1,2]), to eventually
% find out which order describes best the neuron

%

global STACK

nmse();

    function [term_cond, term_score,  n_iters, options] = quick_search()
        [term_cond, term_score,  n_iters, options] = ...
            fit_boo('StopAtStepNumber', 1);
    end

term_score = Inf;
% best_order = [-Inf -Inf];

initial_stack = STACK;

order_x=0;
order_t=0;


% search: try all different orders
for order_x=0:2
%     for order_t=0:2
        % Store the previous state of the STACK and revert to default
        cached_stack = STACK;
        STACK = initial_stack;
        
        % Randomize variables in the fit fields
        for kk = 1:length(STACK)
            if ~isfield(cached_stack{kk}{1}, 'fit_fields')
                continue;
            else
                for ll = 1:length(STACK{kk}{:})
                    % different random initialization depending on the filter
                    if (strcmp(STACK{kk}{ll}.name, 'lindeberg_filter'))
                        % LINDEBERG FILTER
                        mdl = STACK{kk}{ll};
                        params = zeros(7,1);
                        params(1) = mdl.num_dims.*rand(1); %xshift
                        params(2) = (mdl.num_coefs/4).*rand(1); %tshift
                        params(3) = (0.1 + (mdl.num_dims/2-0.1).*rand(1)); %s
                        params(4) = (0.5 + (mdl.num_coefs/3-0.5).*rand(1)); % tau
                        params(5) = 2.*(rand(1)-0.5); % v in [-1;1]
                        params(6) = 10.*(rand(1)-0.5); % norm_factor in [-5;5]
                        STACK{kk}{ll}.('lincoefs') = params;
                        STACK{kk}{ll}.('baseline') = 0;
                        STACK{kk}{ll}.('order_x') = order_x;
                        STACK{kk}{ll}.('order_t') = order_t;
                    elseif (strcmp(STACK{kk}{ll}.name, 'lindeberg_spectral'))
                        % LINDEBERG SPECTRAL-ONLY
                        STACK{kk}{ll}.order = order_x;
                        % default initialization:
                        STACK{kk}{ll}.bf  = 1;
                        STACK{kk}{ll}.s   = 1;
                        STACK{kk}{ll}.add_factor = 0;
                        STACK{kk}{ll}.norm_factor = 1;
                        % random initialization:
%                         STACK{kk}{ll}.bf = STACK{kk}{ll}.num_channels.*rand(1); %xshift
%                         STACK{kk}{ll}.s = (0.1 + (STACK{kk}{ll}.num_channels/2-0.1).*rand(1)); %s
%                         STACK{kk}{ll}.norm_factor = 10.*(rand(1)-0.5); % norm_factor in [-5;5]
%                         STACK{kk}{ll}.add_factor = 0; % add_factor is set to 0
                    end
                end
            end
        end
        
        prefitrnd(1);
        fit05c();
        [~,new_term_score,~,~] = quick_search();
        
        
        fprintf('======================================================================\n');
        fprintf('======================================================================\n');
        fprintf('=================       Finished search at order %d, %d      ==========\n', order_x, order_t);
        fprintf('=================       Score = %f          ==========\n', new_term_score);

        
        if term_score < new_term_score,
            % no improvement => revert
            STACK = cached_stack;
            fprintf('=================                                           ==========\n');
            fprintf('=================                                           ==========\n');
            fprintf('=================        Reverting to old order          ==========\n');
            fprintf('======================================================================\n');
            fprintf('======================================================================\n');
            
        else
            % improvement => keep new (do nothing but update the best score)
            term_score = new_term_score;
            fprintf('=================                                           ==========\n');
            fprintf('=================                                           ==========\n');
            fprintf('=================           Keeping new order            ==========\n');
            fprintf('======================================================================\n');
            fprintf('======================================================================\n');
            
        end
%     end
end

end
