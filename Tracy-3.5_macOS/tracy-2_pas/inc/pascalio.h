  { Interface }

  procedure reset_(var fvar : text;
                   fname : packed array [low..high : integer] of char); external;

  procedure rewrite_(var fvar : text;
                     fname : packed array [low..high : integer] of char);
		     external;


  procedure t2init; external;

  procedure getglobv_(var globval1 : globvalrec); external;

  procedure putglobv_(var globval1 : globvalrec); external;
