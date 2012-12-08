%epochdemo.m
%
%Created: SVD 9/18/01 (from wavtest.m)
%
%
%movgetepochs.m
%readdfile.m
%d2psth.m
%respgetepochs.m

close all

%stimfile='/david/data/R150B/test.review.size.imsm';
%respfile='/auto/k3/david/data/R150B/R150B.review.size.clown.1.d';
%stimfile='/auto/k3/david/data/R170A/test.review.size.clown.imsm';
%respfiles={'/auto/k3/david/data/R170A/R170A.sizereview.clown.1.d',...
%           '/auto/k3/david/data/R170A/R170A.sizereview.clown.2.d'};
stimfile='/auto/k3/david/data/R110A/test.review.canoe.imsm';
respfiles={'/auto/k3/david/data/R110A/R110A.review.canoe.1.d',...
           '/auto/k3/david/data/R110A/R110A.review.canoe.2.d'};
%stimfile='/auto/k3/david/data/R169B/test.review.size.canoe.imsm';
%respfiles={'/auto/k3/david/data/R169B/R169B.sizereview.canoe.1.d',...
%           '/auto/k3/david/data/R169B/R169B.sizereview.canoe.2.d'};

binsperframe=1;   % 14 means 1 ms resolution, 1 means 14 ms bins
svsizemult=1;     % take trials matching these sv parameters
svcolor=0;
svstimclass=0;
compact=0;        % want all the separate trials

% determine times of start and stop in each "fixation" of a review movie
[firstframes,lastframes,mov]=movgetepochs(stimfile);
patches=mov(:,:,firstframes);

% load responses
for ii=1:length(respfiles),
   if strcmp(respfiles{ii}(length(respfiles{ii})-1:length(respfiles{ii})),'.d'),
      r=readdfile(respfiles{ii});
      tr=d2psth(r,binsperframe,svsizemult,svcolor,svstimclass,compact);
   else
      tr=respload(respfiles{ii});
   end
   if ii==1,
      raw_resp=tr;
   else
      raw_resp=[raw_resp tr(:,2:size(tr,2))];
   end
end

% extract corresponding responses during each fixation
[rpsth,rraw]=respgetepochs(raw_resp,firstframes,lastframes,binsperframe);


