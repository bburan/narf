function fh = savethefig(fh, filename)
    global NARF_SAVED_ANALYSIS_PATH;
    pngfile = [NARF_SAVED_ANALYSIS_PATH filesep filename];
    set(gcf,'PaperPositionMode','auto');
    set(gcf,'InvertHardcopy','off');
    %set(0,'defaultAxesFontName', 'Arial');
    %set(0,'defaultTextFontName', 'Arial');
    print(fh, pngfile, '-dpng');
    close(fh);
    fh = nan;
end