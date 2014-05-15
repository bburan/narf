function [celllist,files]=wehr_db()

envpath='/auto/data/daq/wehr/soundfiles/sourcefiles/';
basepath='/auto/users/svd/data/wehr/SpNoise_Data/';
basepathcc='/auto/data/daq/wehr/data/Current Clamp Files (SPNoise)/';
basepathvc='/auto/data/daq/wehr/data/Voltage Clamp Files (SPNoise)/';
basepathcc2='/auto/data/daq/wehr/data/NewData_041414/SPNoise CurrentClamp/';
basepathvc2='/auto/data/daq/wehr/data/NewData_041414/SPNoise VoltageClamp/';
basepathvc3='/auto/data/daq/wehr/data/NewData_041414/SPNoise VoltageClamp/ExcCurrentsOnly/';

files=struct();
files.out121012_002_002.file=[basepathcc 'out121012-002-002.mat'];
files.out121012_002_002.respfmt=0;
files.out121412_007_002.file=[basepathcc 'out121412-007-002.mat'];
files.out121412_007_002.respfmt=0;
files.out121812_002_003.file=[basepathcc 'out121812-002-003.mat'];
files.out121812_002_003.respfmt=0;
files.out121812_002_005.file=[basepathcc 'out121812-002-005.mat'];
files.out121812_002_005.respfmt=0;
files.out121812_004_003.file=[basepathcc 'out121812-004-003.mat'];
files.out121812_004_003.respfmt=0;
files.out121812_004_005.file=[basepathcc 'out121812-004-005.mat'];
files.out121812_004_005.respfmt=0;
files.out122012_003_001.file=[basepathcc 'out122012-003-001.mat'];
files.out122012_003_001.respfmt=0;
files.out012213_007_002.file=[basepathcc 'out012213-007-002.mat'];
files.out012213_007_002.respfmt=0;
files.out031113_002_002.file=[basepathcc 'out031113-002-002.mat'];
files.out031113_002_002.respfmt=0;
files.out040213_003_002.file=[basepathcc 'out040213-003-002.mat'];
files.out040213_003_002.respfmt=0;
files.out040213_004_001.file=[basepathcc 'out040213-004-001.mat'];
files.out040213_004_001.respfmt=0;
files.out051113_003_004.file=[basepathcc 'out051113-003-004.mat'];
files.out051113_003_004.respfmt=0;
files.out051113_003_009.file=[basepathcc 'out051113-003-009.mat'];
files.out051113_003_009.respfmt=0;
files.out051113_004_003.file=[basepathcc 'out051113-004-003.mat'];
files.out051113_004_003.respfmt=0;
files.out061913_002_002.file=[basepathcc 'out061913-002-002.mat'];
files.out061913_002_002.respfmt=0;
files.out062113_005_001.file=[basepathcc 'out062113-005-001.mat'];
files.out062113_005_001.respfmt=0;

files.out121112_004_003_E.file=[basepathvc 'out121112-004-003.mat'];
files.out121112_004_003_E.respfmt=1;
files.out121112_004_003_I.file=[basepathvc 'out121112-004-003.mat'];
files.out121112_004_003_I.respfmt=2;
files.out121112_004_005_E.file=[basepathvc 'out121112-004-005.mat'];
files.out121112_004_005_E.respfmt=1;
files.out121112_004_005_I.file=[basepathvc 'out121112-004-005.mat'];
files.out121112_004_005_I.respfmt=2;
files.out122012_002_001_E.file=[basepathvc 'out122012-002-001.mat'];
files.out122012_002_001_E.respfmt=1;
files.out122012_002_001_I.file=[basepathvc 'out122012-002-001.mat'];
files.out122012_002_001_I.respfmt=2;
files.out122012_007_001_E.file=[basepathvc 'out122012-007-001.mat'];
files.out122012_007_001_E.respfmt=1;
files.out122012_007_001_I.file=[basepathvc 'out122012-007-001.mat'];
files.out122012_007_001_I.respfmt=2;

files.out022313_005_002_E.file=[basepathvc 'out022313-005-002.mat'];
files.out022313_005_002_E.respfmt=1;
files.out022313_005_002_I.file=[basepathvc 'out022313-005-002.mat'];
files.out022313_005_002_I.respfmt=2;

files.out022313_008_001_E.file=[basepathvc 'out022313-008-001.mat'];
files.out022313_008_001_E.respfmt=1;
files.out022313_008_001_I.file=[basepathvc 'out022313-008-001.mat'];
files.out022313_008_001_I.respfmt=2;

files.out052213_003_001_E.file=[basepathvc 'out052213-003-001.mat'];
files.out052213_003_001_E.respfmt=1;
files.out052213_003_001_I.file=[basepathvc 'out052213-003-001.mat'];
files.out052213_003_001_I.respfmt=2;

files.out052213_004_001_E.file=[basepathvc 'out052213-004-001.mat'];
files.out052213_004_001_E.respfmt=1;
files.out052213_004_001_I.file=[basepathvc 'out052213-004-001.mat'];
files.out052213_004_001_I.respfmt=2;

files.out071113_003_002_E.file=[basepathvc 'out071113-003-002.mat'];
files.out071113_003_002_E.respfmt=1;
files.out071113_003_002_I.file=[basepathvc 'out071113-003-002.mat'];
files.out071113_003_002_I.respfmt=2;

files.out071113_004_001_E.file=[basepathvc 'out071113-004-001.mat'];
files.out071113_004_001_E.respfmt=1;
files.out071113_004_001_I.file=[basepathvc 'out071113-004-001.mat'];
files.out071113_004_001_I.respfmt=2;

files.out082613_002_002_E.file=[basepathvc 'out082613-002-002.mat'];
files.out082613_002_002_E.respfmt=1;
files.out082613_002_002_I.file=[basepathvc 'out082613-002-002.mat'];
files.out082613_002_002_I.respfmt=2;

files.out082613_003_001_E.file=[basepathvc 'out082613-003-001.mat'];
files.out082613_003_001_E.respfmt=1;
files.out082613_003_001_I.file=[basepathvc 'out082613-003-001.mat'];
files.out082613_003_001_I.respfmt=2;

files.out082713_005_002_E.file=[basepathvc 'out082713-005-002.mat'];
files.out082713_005_002_E.respfmt=1;
files.out082713_005_002_I.file=[basepathvc 'out082713-005-002.mat'];
files.out082713_005_002_I.respfmt=2;

files.out082813_002_003_E.file=[basepathvc 'out082813-002-003.mat'];
files.out082813_002_003_E.respfmt=1;
files.out082813_002_003_I.file=[basepathvc 'out082813-002-003.mat'];
files.out082813_002_003_I.respfmt=2;

files.out121812_003_002_E.file=[basepathvc 'out121812-003-002.mat'];
files.out121812_003_002_E.respfmt=1;
files.out121812_003_002_I.file=[basepathvc 'out121812-003-002.mat'];
files.out121812_003_002_I.respfmt=2;

files.out021714_002_001.file=[basepathcc2 'out021714-002-001.mat'];
files.out021714_002_001.respfmt=0;
files.out021714_004_002.file=[basepathcc2 'out021714-004-002.mat'];
files.out021714_004_002.respfmt=0;
files.out021714_005_001.file=[basepathcc2 'out021714-005-001.mat'];
files.out021714_005_001.respfmt=0;
files.out021714_005_003.file=[basepathcc2 'out021714-005-003.mat'];
files.out021714_005_003.respfmt=0;
files.out021914_003_004.file=[basepathcc2 'out021914-003-004.mat'];
files.out021914_003_004.respfmt=0;
files.out031914_002_002.file=[basepathcc2 'out031914-002-002.mat'];
files.out031914_002_002.respfmt=0;
files.out031914_002_005.file=[basepathcc2 'out031914-002-005.mat'];
files.out031914_002_005.respfmt=0;
files.out032014_002_002.file=[basepathcc2 'out032014-002-002.mat'];
files.out032014_002_002.respfmt=0;
files.out032014_004_002.file=[basepathcc2 'out032014-004-002.mat'];
files.out032014_004_002.respfmt=0;
files.out032014_004_006.file=[basepathcc2 'out032014-004-006.mat'];
files.out032014_004_006.respfmt=0;

files.out011714_002_005_E.file=[basepathvc2 'out011714-002-005.mat'];
files.out011714_002_005_E.respfmt=1;
files.out011714_002_005_I.file=[basepathvc2 'out011714-002-005.mat'];
files.out011714_002_005_I.respfmt=2;
files.out012814_002_003_E.file=[basepathvc2 'out012814-002-003.mat'];
files.out012814_002_003_E.respfmt=1;
files.out012814_002_003_I.file=[basepathvc2 'out012814-002-003.mat'];
files.out012814_002_003_I.respfmt=2;
files.out012914_002_004_E.file=[basepathvc2 'out012914-002-004.mat'];
files.out012914_002_004_E.respfmt=1;
files.out012914_002_004_I.file=[basepathvc2 'out012914-002-004.mat'];
files.out012914_002_004_I.respfmt=2;
files.out021714_004_003_E.file=[basepathvc2 'out021714-004-003.mat'];
files.out021714_004_003_E.respfmt=1;
files.out021714_004_003_I.file=[basepathvc2 'out021714-004-003.mat'];
files.out021714_004_003_I.respfmt=2;
files.out032014_004_003_E.file=[basepathvc2 'out032014-004-003.mat'];
files.out032014_004_003_E.respfmt=1;
files.out032014_004_003_I.file=[basepathvc2 'out032014-004-003.mat'];
files.out032014_004_003_I.respfmt=2;
%files._E.file=[basepathvc2 ''];
%files._E.respfmt=1;
%files._I.file=[basepathvc2 ''];
%files._I.respfmt=2;

% E only cells
files.out011714_002_002_E.file=[basepathvc3 'out011714-002-002.mat'];
files.out011714_002_002_E.respfmt=1;
files.out012014_004_003_E.file=[basepathvc3 'out012014-004-003.mat'];
files.out012014_004_003_E.respfmt=1;
files.out021714_002_003_E.file=[basepathvc3 'out021714-002-003.mat'];
files.out021714_002_003_E.respfmt=1;
files.out031814_002_004_E.file=[basepathvc3 'out031814-002-004.mat'];
files.out031814_002_004_E.respfmt=1;
files.out031914_002_004_E.file=[basepathvc3 'out031914-002-004.mat'];
files.out031914_002_004_E.respfmt=1;

celllist=fieldnames(files);

