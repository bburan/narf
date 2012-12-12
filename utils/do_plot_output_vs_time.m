function do_plot_output_vs_time(stack, xxx)
    mdl = stack{end};
    x = xxx{end};

    do_plot_time_series(stack, xxx, mdl.input_time, mdl.output);
    
end