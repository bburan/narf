# is this a comment line?

USE cell;

DROP TABLE IF EXISTS gPenetration;
DROP TABLE IF EXISTS gCellMaster;
DROP TABLE IF EXISTS gSingleCell;
DROP TABLE IF EXISTS gSingleRaw;
DROP TABLE IF EXISTS gDataRaw;
DROP TABLE IF EXISTS gRunClass;
DROP TABLE IF EXISTS gUserPrefs;
DROP TABLE IF EXISTS gNote;
DROP TABLE IF EXISTS gNoteCat;

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
    eye                 varchar(10),
    mondist             double(6,6),
    etudeg              double(6,6),
    racknotes           text,
    probenotes          text,
    electrodenotes      text,
    impedance           integer,
    impedancenotes      text,
    stability           integer,
    stabilitynotes      text,
    crap                integer,
    training            integer         default 0,
    
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
    area                varchar(50),
    training            integer         default 0,

    depth               integer,
    umperdepth          double          default 0.5,
    findtime            varchar(50),
    polarity            varchar(10),
    handplot            text,
    comments            text,
    descentnotes        text,
    
    rfppd               double          default 0,
    rfsource            varchar(255),
    rfsize              integer         default 0,
    xoffset             integer         default 0,
    yoffset             integer         default 0,
    eyecal              varchar(255),
    quality             integer         default 0,
    crap                integer         default 0,
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
    unit                integer,

    area                varchar(50),
    handplot            text,
    rfsource            varchar(255),
    rfsize              integer         default 0,
    xoffset             integer         default 0,
    yoffset             integer         default 0,
    quality             integer         default 0,
    crap                integer         default 0,
    latency             integer,
    
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
    isolation           integer,
    
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
    training            integer         default 0,
    
    stimpath            varchar(255),
    stimfile            varchar(255),
    resppath            varchar(255),
    respfile            varchar(255),
    plexonfile          varchar(255),
    matlabfile          varchar(255),
    eyecalfile          varchar(255),
    reps                integer,
    seclength           integer,
    
    time                varchar(50),
    isolation           integer,
    snr                 double,
    syncpulse           integer,
    monitorfreq         double,
    stimconf            integer,
    healthy             integer,
    eyewin              double,
    timejuice           double,
    fixtime             integer,
    maxrate             integer,
    comments            varchar(255),
    corrtrials          integer,
    trials              integer,
    
    addedby             varchar(50),
    info                varchar(255),
    lastmod             timestamp,
    
    primary key         (id),
    index               masterididx (masterid),
    index               cellididx (cellid),
    index               runclassididx (runclassid,stimspeedid)
);

CREATE TABLE gRunClass (
    id                  integer         not null,
    name                varchar(20),
    details             varchar(255),

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
    index               useridx (userid)
);


CREATE TABLE gNote (
    id                  integer         unsigned auto_increment not null, 
    user                varchar(50),
    postdate            datetime,
    moddate             timestamp,
    moduser             varchar(50),
    modallow            integer default 0,

    title               varchar(255),
    replyto             integer,
    minseclevel         integer,
    lab                 varchar(255),
    font                varchar(50),
        
    note                longtext,
    
    primary key         (id)
);

CREATE TABLE gNoteCat (
    id                  integer         unsigned auto_increment not null, 
    noteid              integer,
    user                varchar(50),
    category            varchar(255),
    
    primary key         (id),
    index               noteidx (noteid),
    index               catidx (category)
);

