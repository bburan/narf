USE cell;

DROP TABLE IF EXISTS sCellFile;
DROP TABLE IF EXISTS sRunData;
DROP TABLE IF EXISTS sBatch;
DROP TABLE IF EXISTS sResults;

CREATE TABLE sCellFile (
    id                  integer         unsigned auto_increment not null, 
    cellid              varchar(15)     not null,
    masterid            integer         not null,
    singleid            integer         not null,
    singlerawid         integer         not null,
    rawid               integer         not null,
    celldataid          integer         default 0,
    unit                integer,
    channum             integer,
    model               integer         default 0,
    
    runclassid          integer,
    stimspeedid         integer,
    path                varchar(255),
    resplen             integer         default 0,
    repcount            integer         default 0,
    spikes              integer         default 0,
    a_state             varchar(255),
    area                varchar(255),
    
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
    stimsnr             integer         default 1000,
    
    addedby             varchar(50),
    info                varchar(255),
    lastmod             timestamp,
    
    primary key         (id),
    index               masterididx (masterid),
    index               rawididx (rawid),
    index               celldataididx (celldataid),
    index               cellididx (cellid)
);

CREATE TABLE sRunData(
    id                  integer         unsigned auto_increment not null,
    celldataid          integer         default 0 not null,
    masterid            integer,    
    singleid            integer,    
    cellid              varchar(15),
    grresid             integer	        default 0,

    batch               integer         default 0,
    complete            integer         default 0,
    rcversion           integer         default 0,
    rundate             datetime,
    respath             varchar(255),
    kernfile            varchar(255),
    resfile             varchar(255),
    archive             varchar(50),

    primary key 	(id),
    index 		celldataididx (celldataid)
);

CREATE TABLE sBatch (
    id                  integer         unsigned auto_increment not null,
    name                varchar(20),
    details             varchar(255),
    matcmd              varchar(255),
    
    runclassid          integer         default 0,
    stimspeedid         integer         default 0,
    stimfmtcode         integer         default 0,
    respfmtcode         integer         default 0,
    attstate            integer         default 0,
    
    resploadcmd         varchar(255)    default "loadresp",
    resploadparms       varchar(255),
    respfiltercmd       varchar(255)    default "respresampfull",
    respfilterparms     varchar(255),
    stimloadcmd         varchar(255)    default "loadimfile",
    stimloadparms       varchar(255),
    stimfiltercmd       varchar(255)    default "",
    stimfilterparms     varchar(255),
    stimwindowcrf       double(16,6)    default 1,
    kernfmt             varchar(255)    default "space",
    
    minlag              integer         default -8,
    maxlag              integer         default 12,
    resampcount         integer         default 20,
    resampfmt           integer         default 0,
    
    expfrac             double(16,6)    default 0.0,
    fitfrac             double(16,6)    default 0.1,
    predfrac            double(16,6)    default 0.1,
    predbatch           varchar(255),

    decorrspace         integer         default 1,
    decorrtime          integer         default 1,
    srfiltsigma         double(16,6)    default 0,
    hfiltsigma          double(16,6)    default 0,
    sffiltsigma         double(16,6)    default 0,
    sffiltthresh        double(16,6)    default 0,
    sffiltsmooth        integer         default 0,
    predsmoothsigma     double(16,6)    default 0,
    predtype            integer	        default 0,
    sfscount            integer         default 60,
    sfsstep             integer         default 1,
    nloutparm           integer         default 1,
    parmstring          text,
    
    primary key         (id)
);

CREATE TABLE sResults (
    id                  integer         unsigned auto_increment not null,
    runid               integer,
    batch               integer,
    
    matstr              longtext,

    primary key         (id)
);

