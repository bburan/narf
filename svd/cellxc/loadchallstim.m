function mov = loadchallstim(sImFile,startframe,endframe);

load(sImFile,'stim');

if ~exist('stim','var'),
   load(sImFile,'vstim');
   if exist('vstim','var'),
      stim=vstim;
   else
      disp('loadchallstim.m: no stim or vstim???')
      keyboard
   end
end

if exist('vstim','var'),
   clear vstim
end

if exist('endframe','var') & endframe>0,
   if endframe>size(stim,3),
      endframe=size(stim,3);
   end
   mov=stim(:,:,startframe:endframe);
else
   mov=stim;
end
