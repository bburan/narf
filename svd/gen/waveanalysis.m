function waveanalysis(storo)
%************** routine to look at waveform data in isolation: October 20, 2006
%format: waveanalysis(storo)
%  storo.waveforms = 179x32 matrix of 179 recorded neuron waveforms, each
%                    waveform is 32 samples representing 800us (25us per
%                    sample)
%  storo.class = classifications of units from the experiment
%                0 => non-visually responsive unit, 99 of them in total
%                1 => unusual waveform shape, 6 of them in total
%                2 => low rate, unable to classify ISI, 16 units
%                3 => narrow waveform, inhibitory unit?, 16 in total
%                4 => broad waveform, non-bursty firing, 27 in total
%                5 => broad waveform, bursty firing, 15 in total
% what this routine should do:
%   make a plot first that is similar to what John showed at SFN,
%    the plot should include 1) excluded waveform shapes, 2) the
%    waveforms of all included units, 3) this histogram of durations
%    for only the subset of visually responsive neurons where ISI was
%    classified (the histogram of waveforms John showed in his talk)
%  then make the same plots including the entire population

%************** spline interpolate your data **************************
N = size(storo.waveforms,1);
dTT = 25;
dT = 1.25;  % discretize to 5 micro-secs
for ii = 1:N
   yy = spline(1:32,storo.waveforms(ii,:),[1:(dT/dTT):32]);
   owaveforms(ii,:) = yy;
end
storo.waveforms = owaveforms;
T = size(storo.waveforms,2);  % time samples of waveforms
%**********************************************************************

%************ normalize all waveform amplitudes to size of 1 ******
for i = 1:N
    wave = storo.waveforms(i,:);
    maxo = max(wave);
    mino = min(wave);
    wave = wave / (maxo-mino);
    storo.waveforms(i,:) = wave;
end
%*************** compute durations of waveforms, peak to valley ***
durations = [];
for i = 1:N
  wave = storo.waveforms(i,:);
  ypeak = find( wave == max(wave) );
  ymin = find( wave == min(wave) );
  wdur = (ypeak(1)-ymin(1))*dT;  % time in u-secs from valley to peak
  durations = [durations wdur];
end
storo.durations = durations;
%******* compute which waveforms should be excluded *************
allwaves = storo.waveforms;
excluded = [];
%*********** criteria to identify some non-standard waveform shape
tooearly = 1; %floor( 200 / dT );  % do not consider peaks before 200ms                 
for i = 1:N
    wave = storo.waveforms(i,:);
    wave = storo.waveforms(i,:);
    y = find( wave == min(wave(tooearly:T)) );
    tmin = y(1)*dT;
    y = find( wave == max(wave(tooearly:T)) );
    tmax = y(1)*dT;
    if ((tmin<200) | (tmin>400) | (tmax<300))
         excluded = [excluded i];
    end
end
storo.excluded = excluded;

figure(1);
if (0)  % change this to zero if you want to skip to full analysis
  %********** plot all the excluded waveforms **********
  subplot(3,2,1);
  X = (1:T) .* dT;  %time scale in u-secs
  for i = 1:size(excluded,2)
    plot(X,storo.waveforms(excluded(i),:),'k-'); hold on;
  end
  xlabel('Time (u-sec)');
  ylabel('Normalized Amplitude');
  title('Excluded Waveforms');
  %********** plot all visually responsive, ISI categorized, included waveforms ******
  subplot(3,2,3);
  X = (1:T) .* dT;  %time scale in u-secs
  for i = 1:N
    if (storo.class(i) >= 3)
        plot(X,storo.waveforms(i,:),'k-'); hold on;
    end
  end
  xlabel('Time (u-sec)');
  ylabel('Normalized Amplitude');
  title('Included Waveforms');
  
  %********** plot a histogram of all visually responsive, and
  %********** ISI categorized, included waveforms ******
  subplot(3,2,5);
  inunits{1} = find( (storo.class == 3) );
  inunits{2} = find( (storo.class == 4) );
  inunits{3} = find( (storo.class == 5) );
  %*************** compute histograms on these classes
  vx = 75:dT:525;

  duros = [storo.durations(inunits{1}) storo.durations(inunits{2}) storo.durations(inunits{3})];
  [dip,p,xlow,xup] = HartigansDipSignifTest(duros,10000);
  disp(sprintf('d:%f p:%f ',dip,p));
  input('check');

  ahist = hist(storo.durations(inunits{1}),vx);
  bhist = hist(storo.durations(inunits{2}),vx);
  chist = hist(storo.durations(inunits{3}),vx);
  mato = [ahist; bhist ; chist]';
  colormap([[1,0,0];[0,0,1];[0.6,0.6,0]]);
  bar(vx,mato,1.0,'stacked'); shading flat; hold on
  %***********************************************
  V = axis;
  axis([80 510 0 V(4)]);
  xlabel('Wave Duration');
  ylabel('# of Cells');
  title(sprintf('Dip: %6.2f p:%6.4f',dip,p));
  %************** how many in each group ************
  H = text(120,0.9*V(4),sprintf('N=%d',size(inunits{1},1)));
  set(H,'Color',[1,0,0]);
  set(H,'Fontsize',8);
  H = text(120,0.8*V(4),sprintf('N=%d',size(inunits{2},1)));
  set(H,'Color',[0,0,1]); 
  set(H,'Fontsize',8);
  H = text(120,0.7*V(4),sprintf('N=%d',size(inunits{3},1)));
  set(H,'Color',[0.6,0.6,0]);
  set(H,'Fontsize',8);
end

if (1)
  %********** plot all the excluded waveforms **********
  subplot(3,2,2);
  X = (1:T) .* dT;  %time scale in u-secs
  for i = 1:size(excluded,2)
    plot(X,storo.waveforms(excluded(i),:),'k-'); hold on;
  end
  xlabel('Time (u-sec)');
  ylabel('Normalized Amplitude');
  title('Excluded Waveforms');
  %********** plot all visually responsive, ISI categorized, included waveforms ******
  subplot(3,2,4);
  allincluded = [];
  X = (1:T) .* dT;  %time scale in u-secs
  for i = 1:N
    if (ismember(i,excluded) == 0)
        plot(X,storo.waveforms(i,:),'k-'); hold on;
        allincluded = [allincluded i];
    end
  end
  xlabel('Time (u-sec)');
  ylabel('Normalized Amplitude');
  title('Included Waveforms');
  
  %********** plot a histogram of all visually responsive, and
  %********** ISI categorized, included waveforms ******
  subplot(3,2,6);
  %***********************************
  duros = [storo.durations(allincluded)];
  [dip,p,xlow,xup] = HartigansDipSignifTest(duros,10000);
  disp(sprintf('Hartigans Test dip-value:%f p:%f',dip,p));
  input('Press return to continue');
  %*************** compute histograms on these classes
  vx = 75:dT:525;
  ahist = hist(storo.durations(allincluded),vx);
  mato = [ahist]';
  bar(vx,mato,1.0,'stacked'); shading flat; hold on;
  %***********************************************
  V = axis;
  axis([80 510 0 V(4)]);
  xlabel('Wave Duration');
  ylabel('# of Cells');
  title(sprintf('All Units, Dip: %6.3f p:%6.4f',dip,p));
   %************** how many in each group ************
  H = text(120,0.8*V(4),sprintf('N=%d',size(allincluded,2)));
  set(H,'Color',[0,0,0]);
  set(H,'Fontsize',8);
end

return;