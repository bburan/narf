function custom_demo(stack,xxx,meta);
    global XXX STACK META
    
    recalc_xxx(1);
    
    infostr=sprintf('Cellid: %s\nModelname: %s\nr_fit: %.3f\nr_test: %.3f',...
                    XXX{1}.cellid,META.modelname,XXX{end-1}.score_train_corr,...
                    XXX{end-1}.score_test_corr);
    msgbox(infostr,XXX{1}.cellid);
    