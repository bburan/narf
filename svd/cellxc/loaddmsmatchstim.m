% function mov = loaddmsmatchstim(sImFile,startframe,stopframe,scale_to_pix, 
%                                 mask_to_pix,crop_to_pix,smallfmt,dosmooth)
%
% only load individual flashed stimulus frames
% hacked out of loadimfile 8/7/03
%
function mov = loaddmsmatchstim(sImFile,startframe,stopframe,varargin)

respfile=[sImFile(1:end-4),'resp.mat'];
[r,picid,frameidx]=resploaddmsmatch(respfile);

if exist('startframe','var') & startframe>0,
   frameidx=frameidx(find(frameidx>=startframe));
end
if exist('stopframe','var') & stopframe>0,
   frameidx=frameidx(find(frameidx<=stopframe));
end

mov=loadimframes(sImFile,frameidx,varargin{:});

