function [LinearConstr,Iterate,nineqcstr,neqcstr,ncstr] = ga_constrValidate(nonlcon,Iterate,Aineq,bineq,Aeq,beq,lb,ub,type,options)
%constrValidate validate nonlinear constraint function and create LinearConstr
%   structure for linear constraints
%
%   Private to GA

%   Copyright 2005-2009 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2009/10/10 20:09:02 $

nineqcstr = size(Aineq,1);
neqcstr   = size(Aeq,1);
ncstr     = nineqcstr + neqcstr;
% Validate nonlinear constraint function
if ~isempty(nonlcon)
    % Evaluate nonlinear constraints for the first time
    try
        [cineq,ceq] = nonlcon(Iterate.x');
        Iterate.cineq = zeros(numel(cineq),1);
        Iterate.ceq = zeros(numel(ceq),1);
        Iterate.cineq(:) = cineq;
        Iterate.ceq(:) = ceq;
    catch userFcn_ME
        gads_ME = MException('globaloptim:constrvalidate:confunCheck', ...
            'Failure in initial user-supplied nonlinear constraint function evaluation.');
            userFcn_ME = addCause(userFcn_ME,gads_ME);
            rethrow(userFcn_ME)
    end
    c = [Iterate.cineq;Iterate.ceq];
    if ~all(isreal(c) & isfinite(c))
        error('globaloptim:constrvalidate:confunNotReal','Constraint function must return real value.');
    end
end

if ~strcmpi(type,'unconstrained')
    % Check linear constraint satisfaction at stopping?
    linconCheck = false;
    LinearConstr = struct('Aineq',Aineq,'bineq',bineq,'Aeq',Aeq,'beq',beq,'lb',lb,'ub',ub, ...
    'feasibleX',Iterate.x,'type',type);
    % Mutation function for constrained GA
    mutationFcn = func2str(options.MutationFcn);
    if ~strcmpi(mutationFcn,'mutationadaptfeasible')
        linconCheck = true;
        % Warn for using 'mutationuniform' or 'mutationgaussian' mutation
        % functions
        if any(strcmpi(mutationFcn,{'mutationuniform','mutationgaussian'}))
           warning('globaloptim:constrvalidate:unconstrainedMutationFcn','You are using ''%s'' mutation function for constrained minimization.\nSolution may be infeasible; use ''@mutationadaptfeasible'' function for constrained minimization.',mutationFcn)
        end
    end
    LinearConstr.linconCheck = linconCheck;
else
    LinearConstr = struct('Aineq',Aineq,'bineq',bineq,'Aeq',Aeq,'beq',beq,'lb',lb,'ub',ub, ...
    'feasibleX',Iterate.x,'type',type,'linconCheck',false);
end
