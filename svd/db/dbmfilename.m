% function [mfilename,evpfilename,rawid]=dbmfilename(globalparams,runclass)
%
% inputs:
% globalparams.siteid is string name of site.
% globalparams.Mode is 'Behaviour', 'Physiology', or 'Behaviour & Physiology'
% runclass is abbreviated name of stim/task, eg FTC, TOR, DMS, etc.
%
% outputs:
% appropriate mfilename and evpfilename for next expt at this site
%
% created SVD 2005-11-21
%
function [mfilename,evpfilename,rawid]=dbmfilename(globalparams,runclass)

dbopen;

[rawdata,site,celldata,rcunit,rcchan]=dbgetsite(globalparams.SiteID);

if length(site)==0,
    globalparams.penname=globalparams.SiteID(1:end-1);
    sql=sprintf('SELECT * FROM gPenetration WHERE penname="%s"',...
                globalparams.penname);
    pdata=mysql(sql);
    
    if length(pdata)==0,
        %yn=questdlg(sprintf('Penetration %s doesn''t exist in cellDB. Create it?',globalparams.penname), ...
        %            'cellDB','Yes','No','Yes');
        yn='Yes';
        if isempty(yn) | strcmp(yn,'Yes'),
            globalparams.penid=dbcreatepen(globalparams);
        else
            globalparams.penid=0;
        end
    else
       globalparams.penid=pdata(1).id;
    end
    
    %yn=questdlg(sprintf('Site %s doesn''t exist in cellDB. Create it?',...
    %                    globalparams.SiteID), ...
    %            'cellDB','Yes','No','Yes');
    yn='Yes';
    if isempty(yn) | strcmp(yn,'Yes'),
       globalparams.masterid=dbcreatesite(globalparams);
    end
else
   globalparams.penid=site(1).penid;
   globalparams.masterid=site(1).id;
end

% guess passive for now
if strcmp(globalparams.Physiology,'No'),
    babbr='t';
elseif strcmp(globalparams.Physiology,'Yes -- Passive'),
    babbr='p';
elseif strcmp(globalparams.Physiology,'Yes -- Behavior'),
    babbr='a';
end

rawcount=length(rawdata);

if babbr=='t',
    rawdata=mysql(['SELECT gDataRaw.* FROM gDataRaw,gCellMaster',...
        ' where penid=',num2str(globalparams.penid),...
        ' AND gDataRaw.masterid=gCellMaster.id']);
    rawcount=length(rawdata);
    mfilename=sprintf('%s%s%s%s_%s_%s_%d.m',globalparams.outpath,...
                      globalparams.Ferret,filesep,...
                      globalparams.Ferret,datestr(now,'yyyy_mm_dd'), ...
                      runclass,rawcount+1);
    evpfilename=sprintf('%s%s%s%s_%s_%s_%d.evp',globalparams.outpath,...
                        globalparams.Ferret,filesep,...
                        globalparams.Ferret,datestr(now,'yyyy_mm_dd'), ...
                        runclass,rawcount+1);
else
    mfilename=sprintf('%s%s%s%s%02d_%s_%s.m',globalparams.outpath,...
                      globalparams.Ferret,filesep,...
                      globalparams.SiteID,rawcount+1,babbr,runclass);
    evpfilename=sprintf('%s%s%s%s%02d_%s_%s.evp',globalparams.outpath,...
                        globalparams.Ferret,filesep,...
                        globalparams.SiteID,rawcount+1,babbr,runclass);
end

if 1,
    prompt={'Save to:                                                         _'};
    name='m-file name';
    numlines=1;
    defaultanswer={mfilename};
    answer=inputdlg(prompt,name,numlines,defaultanswer);
    if ~isempty(answer),
        mfilename=answer{1};
        rawid=dbcreateraw(globalparams,runclass,mfilename,evpfilename);
    else
        mfilename='';
        rawid=-1;
    end
else
    yn=questdlg(sprintf('File %s doesn''t exist in cellDB. Create it?',...
        mfilename), ...
        'cellDB','Yes','No','Cancel','Yes');
    if isempty(yn) | strcmp(yn,'Yes'),
        rawid=dbcreateraw(globalparams,runclass,mfilename,evpfilename);
    elseif strcmp(yn,'Cancel'),
        rawid=-1;
    else
        rawid=0;
    end
end
