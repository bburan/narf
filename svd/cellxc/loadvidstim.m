% function mov = loaddmsstim(sImFile,startframe,stopframe,scale_to_pix, 
%                           mask_to_pix,crop_to_pix,smallfmt,dosmooth)
%
% only load individual flashed stimulus frames
% hacked out of loadimfile 8/7/03
%
function mov = loadvidstim(sImFile,startframe,stopframe,varargin)

respfile=[sImFile(1:end-4),'resp.mat'];
[r,frameidx]=resploadvidpfth(respfile);

if exist('stopframe','var') & stopframe>0,
   frameidx=frameidx(1:stopframe);
end
if exist('startframe','var') & startframe>0,
   frameidx=frameidx(startframe:end);
end

mov=loadimframes(sImFile,frameidx,varargin{:});
