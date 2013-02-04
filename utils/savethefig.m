function fh = savethefig(fh, filename, subdir)
    global NARF_SAVED_ANALYSIS_PATH;
    if ~exist([NARF_SAVED_ANALYSIS_PATH filesep subdir]')
        mkdir([NARF_SAVED_ANALYSIS_PATH filesep subdir]);
    end
    pngfile = [NARF_SAVED_ANALYSIS_PATH filesep subdir filesep filename];
    set(gcf,'PaperPositionMode','auto');
    set(gcf,'InvertHardcopy','off');
    %set(0,'defaultAxesFontName', 'Arial');
    %set(0,'defaultTextFontName', 'Arial');
    print(fh, pngfile, '-dpng');
    close(fh);
    fh = nan;
end