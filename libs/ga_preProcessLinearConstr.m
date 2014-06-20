function [XOUT,Aineq,Bineq,Aeq,Beq,LB,UB,msg,exitflag] = ...
    preProcessLinearConstr(XOUT,Aineq,Bineq,Aeq,Beq,LB,UB,numberOfVariables,type,verbosity)
%PREPROCESSLINEARCONSTR validates dimension of constraint matrices, removes
%   redundancy in them, and finds initial feasible point.
%
%   private to PATTERNSEARCH and GA.

%   Copyright 2003-2009 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2009/10/10 20:09:42 $

% Initialize
msg = '';
exitflag = 1;

% If problem type is unconstrained return here
if strcmpi(type,'unconstrained')
    return;
end

% If called from PFMINBND, do this bound check only
if strcmpi(type,'boundconstraints')
    A = eye(numberOfVariables);
    % Check the box constraints (bounds) first.
    lbound   = A*XOUT-LB >= 0;
    ubound   = A*XOUT-UB <= 0;
    feasible = all(lbound) && all(ubound);
    if ~feasible
        XOUT(~lbound) = LB(~lbound);
        XOUT(~ubound) = UB(~ubound);
    end
    return;
end

% ------------Only linear constrained problems reach here--------------
% Initialize
tol = sqrt(eps);
Xin = XOUT;

% We allow row or column vectors for Beq and Bineq
Bineq = Bineq(:);
Beq   = Beq(:);
% Set the constraints up: defaults and check size
[nineqcstr,n] = size(Aineq);
[neqcstr,m]   = size(Aeq);
if ~isempty(Aineq)
    if ~isequal(length(Bineq),nineqcstr)
        error('globaloptim:preprocesslinearconstr:inconsistentAineqAndBineq','The number of rows in A must be the same as the length of b.')
    elseif ~isequal(numberOfVariables,n)
        error('globaloptim:preprocesslinearconstr:inconsistentAineqAndX0','The number of columns in A must be the same as the length of X0.')
    end
elseif ~isempty(Bineq)
    error('globaloptim:preprocesslinearconstr:emptyAineqNotBineq','The constraint matrices A and b are not consistent.')
end
if ~isempty(Aeq)
    if ~isequal(length(Beq),neqcstr)
        error('globaloptim:preprocesslinearconstr:inconsistentAeqAndBeq','The number of rows in Aeq must be the same as the length of beq.')
    elseif ~isequal(numberOfVariables,m)
        error('globaloptim:preprocesslinearconstr:inconsistentAeqAndX0','The number of columns in Aeq must be the same as the length of X0.')
    end
elseif ~isempty(Beq)
    error('globaloptim:preprocesslinearconstr:emptyAeqNotBeq','The constraint matrices Aeq and beq are not consistent.')
end

% Remove dependent constraint, if any and find a basic solution
[XOUT,Aineq,Bineq,Aeq,Beq,msg,how,exitflag]= eqnsolv(XOUT,Aineq,Bineq,Aeq,Beq,LB,UB,verbosity);

% Is initial point feasible?
if strcmp(how,'infeasible')
    % Equalities are inconsistent, so return original X
    XOUT = Xin;
    return
end

%------------- Find a feasible point-------------
Abox = eye(numberOfVariables);
% Check the bound constraints
lbound   = Abox*XOUT-LB>=0;
ubound   = Abox*XOUT-UB<=0;
feasible = all(lbound) && all(ubound);
if ~feasible
    XOUT(~lbound) = LB(~lbound);
    XOUT(~ubound) = UB(~ubound);
end

% Now add the linear constraints too and check it
feasible = isTrialFeasible(XOUT,Aineq,Bineq,Aeq,Beq,[],[],tol);

% Not a feasible point? find an initial feasible point using LP (Two
% approaches)
if ~feasible
    % Find a feasible initial point using linprog active-set algorithm
    [XOUT,unused,success] = linprog([],Aineq,Bineq,Aeq,Beq,LB,UB,XOUT, ...
        optimset('Simplex','off','LargeScale','off','Display','off'));
    feasible = isTrialFeasible(XOUT,Aineq,Bineq,Aeq,Beq,LB,UB,tol);
    if success <=0 || ~feasible
        % Add a slack variable and find a feasible initial point
        [XOUT,success] = initialfeasible(XOUT,numberOfVariables,Aineq,Bineq,Aeq,Beq,LB,UB);
        if success > 0
            feasible = isTrialFeasible(XOUT,Aineq,Bineq,Aeq,Beq,LB,UB,tol);
        end
    end
    
    % Add additional linprog algorithm
    if success <= 0 || ~feasible
        % Find a feasible initial point using lipsol algorithm
        [XOUT,unused,success] = linprog([],Aineq,Bineq,Aeq,Beq,LB,UB,XOUT, ...
            optimset('Display','off'));
        if success > 0
            feasible = isTrialFeasible(XOUT,Aineq,Bineq,Aeq,Beq,LB,UB,tol);
        end
    end
    if success <= 0 || ~feasible
        % Find a feasible initial point using Simplex algorithm
        [XOUT,unused,success] = linprog([],Aineq,Bineq,Aeq,Beq,LB,UB,XOUT, ...
            optimset('Simplex','on','LargeScale','off','Display','off'));
        if success > 0        
            feasible = isTrialFeasible(XOUT,Aineq,Bineq,Aeq,Beq,LB,UB,tol);
        end
    end
    
    % Quit now
    if success <= 0 || ~feasible
        XOUT = Xin;
        msg = sprintf('Could not find a feasible initial point.');
        exitflag = -2;
    end
end

