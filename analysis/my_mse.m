% Perform the nonlinear fitting
function prediction = my_mse(phi, stack_depth),
    global XXX;
    unpack_fittables(phi);
    fprintf('my_mse() is recalcing XXX{%d}...\nPhi = ', stack_depth);
    disp(phi');
    recalc_xxx(stack_depth);
    prediction = XXX{stack_depth+1}.dat.por024a19_p_SPN.stim(:);
end