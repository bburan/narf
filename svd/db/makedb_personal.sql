DROP TABLE IF EXISTS dMusic;
DROP TABLE IF EXISTS dJuke;
DROP TABLE IF EXISTS dPlaylist;
DROP TABLE IF EXISTS dGenre;
DROP TABLE IF EXISTS dEval;
DROP TABLE IF EXISTS dComment;


CREATE TABLE dMusic (
    id                  integer         unsigned auto_increment not null,
    
    title               varchar(255),
    track               integer,
    artist              varchar(255),
    album               varchar(255),
    file                varchar(255),
    playsec             integer,
    bitrate             integer,
    genre               varchar(255),
    comment             varchar(255),
    note                varchar(255),
    user                varchar(255),
    reqcount            integer         default 0,
    id3bad              integer         default 0,
    crap                integer         default 0,
    dup                 integer         default 0,
    filemissing         integer         default 0,
    hatepoints          integer         default 0,
    source              varchar(50),
    cat1                integer         default 0,
    cat2                integer         default 0,
    modstring           varchar(255)    default NULL,
    temprating          double          default NULL,
    dateadded           datetime,
    
    lastmod             timestamp,
    
    primary key         (id),
    index               titleidx (title),
    index               artistidx (artist),
    index               albumidx (album),
    index               genreidx (genre),
    index               sourceidx (source)
);

CREATE TABLE dJuke (
    
    id                  integer         unsigned auto_increment not null,

    songid              integer,
    song                varchar(255),
    active              integer        default 0,
    command             varchar(255),
    cmduser             varchar(255),
    randsql             text,
    randuser            varchar(255),
    randexp             datetime,
    randnextid          integer        default 0,
    playlist            varchar(50),
    
    mdbsyncdate         datetime,

    createdate		timestamp,
    createdby		varchar(50)	default 'svd',
    daemonclick		datetime,
    daemonhost          varchar(255)

);

CREATE TABLE dPlaylist (
    id                  integer         unsigned auto_increment not null,
    
    songid              integer,
    note                varchar(255),
    user                varchar(255),
    
    queuedate           datetime,
    startdate           datetime,
    complete            integer        default 0,
    progress            integer        default 0,
    priority            integer        default 1,
    killnow             integer        default 0,
    playlist		varchar(50),

    primary key         (id),
    index               songidx (songid),
    index               useridx (user),
    index               startidx (startdate),
    index               queuedateidx (queuedate),
    index               completeidx (complete),
    index               listidx (playlist)
);

CREATE TABLE dGenre (
    id                  integer         unsigned auto_increment not null,
    
    songid              integer,
    genre               varchar(255),
    user                varchar(255),
    private             integer         default 0,
    
    lastmod             timestamp,

    primary key         (id),
    index               songidx (songid),
    index               prividx (user,private),
    index               genreidx (genre)
);

CREATE TABLE dEval (
    id                  integer         unsigned auto_increment not null,
    
    songid              integer,
    playid              integer,
    user                varchar(255),
    uid                 integer,
    
    rating              integer,
    note                varchar(255),
    lastmod             timestamp,
    
    primary key         (id),
    index               songidx (songid),
    index               useridx (user,rating),
    index               uididx (uid)
);

CREATE TABLE dEig (
    id                  integer         unsigned auto_increment not null,
    rr                  integer         not null,
    cc                  integer         not null,
    uid                 integer,
    u                   double,

    primary key         (id),
    index               eigidx (rr,cc)
);

CREATE TABLE dEigProj (
    songid              integer         not null,
    e1                  double,
    e2                  double,
    e3                  double,
    e4                  double,
    e5                  double,
    e6                  double,
    e7                  double,
    e8                  double,
    e9                  double,
    e10                 double,
    st1                 double default 0,
    st2                 double default 0,
    st3                 double default 0,
    st4                 double default 0,
    st5                 double default 0,
    
    primary key         (songid)
);

CREATE TABLE dComment (
    id                  integer         unsigned auto_increment not null,
    
    playid              integer,
    user                varchar(255),
    category            integer         default 0,
    
    note                text,
    lastmod             timestamp,
    
    primary key         (id),
    index               useridx (user)
);

CREATE TABLE dSelectList (
    id                  integer         unsigned auto_increment not null,
    
    songid              integer,
    user                varchar(255),
    
    queuedate           datetime,
    startdate           datetime,
    playcount           integer        default 0,
    priority            double         default 0,
    playlist		varchar(50),

    primary key         (id),
    index               songidx (songid),
    index               listidx (playlist)
);


