# is this a comment line?

USE cell;

DROP TABLE IF EXISTS tGlobalData;
DROP TABLE IF EXISTS tRunResult;
DROP TABLE IF EXISTS tQueue;
DROP TABLE IF EXISTS tComputer;
DROP TABLE IF EXISTS tEvent;
DROP TABLE IF EXISTS tBill;


CREATE TABLE tGlobalData (
    createdate		timestamp,
    createdby		varchar(50)	default 'svd',
    daemonclick		datetime,
    daemonhost          varchar(255)
);

INSERT tGlobalData (createdby) VALUES ('svd');


CREATE TABLE tQueue (
    id                  integer         unsigned auto_increment not null,
    
    rundataid           integer,
    progname            varchar(255),
    parmstring          text,
    user                varchar(255)   default "david",
    
    machinename         varchar(255),
    computerid          integer,
    pid                 integer        default 0,
    progress            integer        default 0,
    complete            integer        default 0,
    priority            integer        default 0,
    waitid              integer        default 0,
    
    queuedate           datetime,
    startdate           datetime,
    lastdate            timestamp,
    killnow             integer        default 0,
    mailto              varchar(255),
    mailcommand         varchar(255),
    allowqueuemaster    integer        default 0,

    primary key         (id)
);


CREATE TABLE tComputer (
    id                  integer         unsigned auto_increment not null,
    
    name                varchar(255)    NOT NULL,
    ext                 varchar(255)    NOT NULL,
    load1               double          default -1,
    load5               double          default -1,
    load15              double          default -1,
    lastdate            timestamp,
    
    location            integer         default 0,
    maxproc             integer         default 2,
    numproc             integer         default 0,
    allowqueuemaster    integer         default 0,
    owner               varchar(255),
    allowothers         integer         default 1,
    killqueueload       double          default 1.3,
    allowqueueload      double          default 0.3,
    lastoverload        integer         default 0,
    pingcount           integer         default 0,
    dead                integer         default 0,
    nocheck             integer         default 0,
    
    macaddr             varchar(255),
    os                  varchar(255),
    hardware            varchar(255),
    note                varchar(255),
    room                varchar(255),
    
    primary key         (id),
    index               nameidx (name)
);


CREATE TABLE tEvent (
    id                  integer         unsigned auto_increment not null,
    
    code                integer         default 0,
    note                varchar(255)    default 0,
    eventdate           timestamp,
    
    user                varchar(255),
    computerid          integer         default 0,
    queueid             integer         default 0,
    pid                 integer         default 0,
    
    primary key         (id),
    index               codeidx (code),
    index               queueidx (queueid)
);


