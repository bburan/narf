% Perform the nonlinear fitting routine, starting at start_depth and
% recalculating until the end of the XXX datastructure. 
% 
function prediction = my_mse(phi, start_depth)
    global XXX;
    unpack_fittables(phi);
    %fprintf('my_mse() is recalcing XXX{%d}...\nPhi = ', start_depth);
    %disp(phi');
    fprintf('.');
    recalc_xxx(start_depth);
    prediction = XXX{end}.dat.(XXX{1}.training_set{1}).stim(:);
end