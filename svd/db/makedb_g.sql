# is this a comment line?

USE test;

DROP TABLE IF EXISTS gAnimal;
DROP TABLE IF EXISTS gPenetration;
DROP TABLE IF EXISTS gCellMaster;
DROP TABLE IF EXISTS gSingleCell;
DROP TABLE IF EXISTS gSingleRaw;
DROP TABLE IF EXISTS gDataRaw;
DROP TABLE IF EXISTS gRunClass;
DROP TABLE IF EXISTS gUserPrefs;
DROP TABLE IF EXISTS sCellFile;

CREATE TABLE gAnimal (
    id 	                integer         unsigned auto_increment not null, 

    animal              varchar(50),
    sex                 varchar(10)     default 'f',
    eartag              varchar(50),
    cellprefix          varchar(50),
    birthday            date,
    caretaker           varchar(50),
    arrivalweight       integer,
    pullweight          integer,
    poleweight          integer,
    onschedule          integer         default 0,
    retired             integer         default 0,
    medical             text,
    notes               text,
    
    addedby             varchar(50),
    info                varchar(255),
    lastmod             timestamp,

    primary key         (id)
);


CREATE TABLE gPenetration (
    id 	                integer         unsigned auto_increment not null, 
    penname             varchar(50)     not null,

    animal              varchar(15),
    well                integer         default 0,
    
    pendate             varchar(50),
    who                 varchar(50),
    fixtime             varchar(50),
    water               double(6,6),
    weight              double(6,6),
    ear                 varchar(10),
    numchans            integer,
    racknotes           text,
    speakernotes        text,
    probenotes          text,
    electrodenotes      text,
    crap                integer,
    training            integer         default 0,
    
    impedance           integer,
    impedancenotes      text,
    stability           integer,
    stabilitynotes      text,
    eye                 varchar(10),
    mondist             double(6,6),
    etudeg              double(6,6),

    descentnotes        text,
    
    addedby             varchar(50),
    info                varchar(255),
    lastmod             timestamp,

    primary key         (id)
);

CREATE TABLE gCellMaster (
    id 	                integer         unsigned auto_increment not null, 
    siteid              varchar(15)     not null,
    cellid              varchar(15)     not null,
    penid               integer         not null, 
    penname             varchar(50),
    
    animal              varchar(15),
    well                integer         default 0,
    training            integer         default 0,
    
    depth               varchar(255),
    umperdepth          double          default 0.5,
    findtime            varchar(50),
    polarity            varchar(10),
    comments            text,
    crap                integer         default 0,
    
    descentnotes        text,
    
    handplot            text,
    area                varchar(255),
    rfppd               double          default 0,
    rfsource            varchar(255),
    rfsize              integer         default 0,
    xoffset             integer         default 0,
    yoffset             integer         default 0,
    eyecal              varchar(255),
    quality             integer         default 0,
    latency             integer         default 0,
    
    addedby             varchar(50),
    info                varchar(255),
    lastmod             timestamp,

    primary key         (id),
    index               cellididx (cellid)
);

CREATE TABLE gSingleCell (
    id 	                integer         unsigned auto_increment not null, 
    siteid              varchar(15)     not null,
    cellid              varchar(15)     not null,
    masterid            integer         not null,
    penid               integer         not null,
    rawid               integer         not null,
    
    channel             varchar(2),
    channum             integer         default 1,
    unit                integer         default 1,

    area                varchar(50),
    handplot            text,
    quality             integer         default 0,
    crap                integer         default 0,
    latency             integer,
    cf                  integer,
    bw                  integer,

    rfsource            varchar(255),
    rfsize              integer         default 0,
    xoffset             integer         default 0,
    yoffset             integer         default 0,
    
    addedby             varchar(50),
    info                varchar(255),
    lastmod             timestamp,
    
    primary key         (id),
    index               cellididx (cellid),
    index               masteridx (masterid)
);

CREATE TABLE gSingleRaw (
    id 	                integer         unsigned auto_increment not null, 
    cellid              varchar(15)     not null,
    masterid            integer         not null,
    singleid            integer         not null,
    penid               integer         not null,
    rawid               integer         not null,
    
    channel             varchar(2),
    unit                integer,
    
    crap                integer         default 0,
    isolation           double,
    
    addedby             varchar(50),
    info                varchar(255),
    lastmod             timestamp,
    
    primary key         (id),
    index               singleidx (singleid),
    index               cellididx (cellid),
    index               masteridx (masterid)
);


CREATE TABLE gDataRaw (
    id                  integer         unsigned auto_increment not null, 
    cellid              varchar(15)     not null,
    masterid            integer         not null,
    
    runclassid          integer         not null,
    runclass            varchar(255),    
    stimspeedid         double,
    task                varchar(255),
    behavior            varchar(255),
    stimclass           varchar(255),
    training            integer         default 0,
    bad                 integer         default 0,
    
    parmfile            varchar(255),
    respfileevp         varchar(255),
    respfileraw         varchar(255),
    respfile            varchar(255),
    reps                integer,

    stimpath            varchar(255),
    stimfile            varchar(255),
   
    timejuice           double,
    parameters          text,
    comments            text,
    corrtrials          integer,
    trials              integer,

    resppath            varchar(255),
    matlabfile          varchar(255),
    eyecalfile          varchar(255),
    plexonfile          varchar(255),
    maxrate             integer,
    fixtime             integer,
    seclength           integer,
    time                varchar(50),
    isolation           integer,
    snr                 double,
    syncpulse           integer,
    monitorfreq         double,
    stimconf            integer,
    healthy             integer,
    eyewin              double,
    
    addedby             varchar(50),
    info                varchar(255),
    lastmod             timestamp,
    
    primary key         (id),
    index               masterididx (masterid),
    index               cellididx (cellid),
    index               runclassididx (runclassid,stimspeedid),
    index               behavioridx (behavior,task),
    index               stimclassididx (stimclass)
);

CREATE TABLE gRunClass (
    id                  integer         unsigned auto_increment not null,
    name                varchar(20),
    details             varchar(255),
    task                varchar(255),
    stimclass           varchar(255),
    
    addedby             varchar(50),
    info                varchar(255),
    lastmod             timestamp,

    primary key         (id),
    index               nameidx (name)
);

CREATE TABLE gUserPrefs (
    id                  integer         unsigned auto_increment not null, 
    userid              varchar(50),
    password            varchar(255),

    lab                 varchar(255),
    lastanimal          varchar(50),
    lastwell            integer,
    lasttraining        integer         default 0,
    lastpen             varchar(50),
    dataroot            varchar(255),
    seclevel            integer         default 0,
    email               varchar(255),
    realname            varchar(255),
    
    lastallowqueuemaster  integer         default 0,
    lastmachinesort     varchar(255)    default "tComputer.load1",
    lastjobcomplete     integer         default -1,
    lastjobuser         varchar(255)    default "",

    bgcolor             varchar(255)    default "#FFFFFF",
    fgcolor             varchar(255)    default "#000000",
    linkfg              varchar(255)    default "#6666FF",
    vlinkfg             varchar(255)    default "#3333FF",
    alinkfg             varchar(255)    default "#FFFF00",
    
    avgrating           double          default 0,
    stdrating           double          default 1,
    temprat             double          default 0,
    birthday            date            default NULL,
    playlist            varchar(255),
    p1 double, p2 double, p3 double, p4 double, p5 double, 
    p6 double, p7 double, p8 double, p9 double, p10 double,
    
    primary key         (id),
    index               useridx (userid),
    index               labidx (lab)
);

CREATE TABLE sCellFile (
    id                  integer         unsigned auto_increment not null, 
    cellid              varchar(15)     not null,
    masterid            integer         not null,
    singleid            integer         not null,
    singlerawid         integer         not null,
    rawid               integer         not null,
    celldataid          integer         default 0,
    
    runclassid          integer,
    stimspeedid         integer,
    path                varchar(255),
    resplen             integer         default 0,
    repcount            integer         default 0,
    spikes              integer         default 0,
    a_state             varchar(255),
    
    respfile            varchar(255),
    respvarname         varchar(50),
    respfiletype        integer         default 1,
    nosync              integer         default 0,
    respfilefmt         varchar(50),
    respfmtcode         integer         default -1,
    
    stimfile            varchar(255),
    stimpath            varchar(255),
    stimfiletype        integer         default 1,
    stimiconside        varchar(255),
    stimfilecrf         double          default 0,
    stimwindowsize      integer         default 0,
    stimfilefmt         varchar(50),
    stimfmtcode         integer         default -1,
    
    addedby             varchar(50),
    info                varchar(255),
    lastmod             timestamp,
    
    primary key         (id),
    index               masterididx (masterid),
    index               rawididx (rawid),
    index               celldataididx (celldataid),
    index               cellididx (cellid)
);


