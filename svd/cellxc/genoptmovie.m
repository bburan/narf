% function [indexpath,indexfile]=genoptmovie(kern,movformat,framecount);
%
% created SVD 6/5/02
%
function [indexpath,indexfile]=genoptmovie(kern,movformat,framecount);

rcsetstrings;

VER='synth0.1';
STIMPATH='/auto/pcp1/exp/movies/optimal/';
%STIMPATH='/auto/h1/exp/stimarchive/REESE14/optimal/';
%STIMPATH='/auto/lsd1/david/stimarchive/REESE13/optimal/';
MEANFRAME=8;  % approx 8 Hz
cellid=movformat.cellid;

if ~exist('framecount'),
   framecount=60;
else
   framecount=framecount/2; % hacky adjust framecount to match ben's
end
fprintf('Generating %d frames (%d opt, %d random)\n',...
	framecount*2,framecount,framecount);

% generate an indexpath that doesn't overwrite and existing one
newdir=[cellid,'-',runclassstr{movformat.cellfiledata(1).runclassid+1},...
        '-',movformat.cellfiledata(1).stimfilefmt,'-synth-1/'];
idxcount=1;
while exist([STIMPATH,newdir],'dir')
   idxcount=idxcount+1;
   newdir=[cellid,'-',runclassstr{movformat.cellfiledata(1).runclassid+1},...
           '-',movformat.cellfiledata(1).stimfilefmt,'-synth-',...
           num2str(idxcount),'/'];
end
mkdir(STIMPATH,newdir);

indexpath=[STIMPATH,newdir];
indexfile='indexopt';

[mov,transmov]=kern2opt(kern,movformat,framecount);

movshuff=zeros(size(mov));
transmovshuff=zeros(size(transmov));
for ii=1:framecount,
   kernshuff=shuffle(kern);
   [movshuff(:,:,ii),transmovshuff(:,ii)]=kern2opt(kernshuff,movformat,1);
end

framedur=floor(rand(framecount*2,1)*3)+MEANFRAME-1;

fid=fopen([indexpath,indexfile],'w');
fprintf(fid,'# stim %s\n# pix %d\n# framecount %d\n',...
        VER,movformat.stimwindowsize,...
        sum(framedur));

for ii=1:framecount,
   fn=sprintf('image_%.4d.pgm',ii*2-1);
   pgmWrite(mov(:,:,ii),[indexpath,fn]);
   fprintf(fid,'%s %d\n',fn,framedur(ii*2-1));

   fn=sprintf('image_%.4d.pgm',ii*2);
   pgmWrite(movshuff(:,:,ii),[indexpath,fn]);
   fprintf(fid,'%s %d\n',fn,framedur(ii*2));
end

% save transformed stim to imsm file
imfile=[newdir(1:(end-1)),'.imsm'];
spacecount=size(transmov,1);
transmovout=reshape([transmov; transmovshuff],spacecount,framecount*2);
writeimfile(transmovout,[indexpath,imfile],4);

fprintf('indexpath: %s\nindexfile: %s\n',indexpath,indexfile);





