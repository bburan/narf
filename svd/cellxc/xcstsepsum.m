
figure

%% FIRST - regular vs space-time separable within class

subplot(3,3,1);
[xc,celllist,p]=batchcomp([23;86],[2;4],[1;1],0);
title(sprintf('review full v stsep w/in p<%.2f',p(end)));

subplot(3,3,2);
[xc,celllist,p]=batchcomp([29;87],[2;1],[2;4],0);
title(sprintf('gratrev full v stsep w/in p<%.2f',p(end)));

%% SECOND - regular vs rev-time pred rev

subplot(3,3,3);
[xc,celllist,p]=batchcomp([29;86],[2;2],[1;1],0);
title(sprintf('gratrev sfft full vs revtime p<%.2f',p(end)));

%% THIRD - rev-time predicting review cross-class comp

subplot(3,3,4);
[xc,celllist,p]=batchcomp([86;86],[2;4],[1;1],0);
title(sprintf('grspace revtime v review stsep p<%.2f',p(end)));

subplot(3,3,5);
[xc,celllist,p]=batchcomp([83;83],[2;1],[1;1],0);
title(sprintf('grpos revtime v revpos stsep p<%.2f',p(end)));

subplot(3,3,6);
[xc,celllist,p]=batchcomp([86;83],[4;1],[1;1],0);
title(sprintf('rev stsep v revpos p<%.2f',p(end)));

subplot(3,3,7);
[xc,celllist,p]=batchcomp([86;85],[4;1],[1;1],0);
title(sprintf('rev stsep v revneg p<%.2f',p(end)));

subplot(3,3,8);
[xc,celllist,p]=batchcomp([87;84],[1;1],[4;2],0);
title(sprintf('gr-pos v gr-st w/in p<%.2f',p(end)));

subplot(3,3,9);
[xc,celllist,p]=batchcomp([86;83],[2;2],[1;1],0);
title(sprintf('grspace v gr-pos revtime p<%.2f',p(end)));

%subplot(3,3,9);
%[xc,celllist,p]=batchcomp([86;85],[2;2],[1;1],0);
%title(sprintf('grspace v gr-neg revtime p<%.2f',p(end)));




set(gcf,'PaperPosition',[0.25 0.25 8 10.5],'PaperOrientation','portrait');

return


figure

%% FIRST - regular vs space-time separable within class

subplot(3,3,1);
[xc,celllist,p]=batchcomp([23;86],[2;1],[1;4],0);
title(sprintf('review full v stsep w/in p<%.2f',p(end)));

subplot(3,3,2);
[xc,celllist,p]=batchcomp([32;88],[2;1],[3;4],0);
title(sprintf('natrev full v stsep w/in p<%.2f',p(end)));

subplot(3,3,3);
[xc,celllist,p]=batchcomp([29;87],[2;1],[2;4],0);
title(sprintf('gratrev full v stsep w/in p<%.2f',p(end)));

%% SECOND - regular vs space-time separable predicitng review

subplot(3,3,4);
[xc,celllist,p]=batchcomp([23;86],[2;1],[1;4],0);
title(sprintf('review sfft full vs revtime p<%.2f',p(end)));

subplot(3,3,5);
[xc,celllist,p]=batchcomp([32;86],[2;1],[1;3],0);
title(sprintf('natrev sfft full vs revtime p<%.2f',p(end)));

subplot(3,3,6);
[xc,celllist,p]=batchcomp([29;86],[2;1],[1;2],0);
title(sprintf('gratrev sfft full vs revtime p<%.2f',p(end)));

%% THIRD - rev-time predicting review cross-class comp

subplot(3,3,7);
[xc,celllist,p]=batchcomp([86;86],[1;1],[2;4],0);
title(sprintf('gratrev v review sfft stsep p<%.2f',p(end)));

subplot(3,3,8);
[xc,celllist,p]=batchcomp([86;86],[1;1],[3;4],0);
title(sprintf('natrev v review sfft stsep p<%.2f',p(end)));

subplot(3,3,9);
[xc,celllist,p]=batchcomp([86;86],[1;1],[2;3],0);
title(sprintf('gratrev v natrev sfft stsep p<%.2f',p(end)));


set(gcf,'PaperPosition',[0.25 0.25 8 10.5],'PaperOrientation','portrait');
