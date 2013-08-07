function narf_strf_plot(STACK,XXX,META,sel_results);
    
    firmod=find_modules(STACK,'fir_filter');
    
    coefs=firmod{1}{1}.coefs;
    rasterfs=STACK{1}{1}.raw_resp_fs;
    numchannels=size(coefs,1);
    ff=round(2.^linspace(log2(125),log2(16000),numchannels));
    
    fh=figure;
    plotastrf(coefs,4,ff,rasterfs);
    
    title(sprintf('%s - %.3f',XXX{1}.cellid,sel_results.r_test));
    
    outpath=['/auto/data/code/saved_images/narf_strf_plot/' ...
             num2str(sel_results.batch) '/'];
    drawnow;
    print(['-f',num2str(fh)],'-depsc',[outpath sel_results.cellid '.strf.eps']);
