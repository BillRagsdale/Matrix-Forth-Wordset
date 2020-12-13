\  WFR  Nov. 17, 2020

\ G: matrices are located by their storage address
\ H: Changes to { on data address, {{ on descriptor
\ i: 2020-11-19 revise conformal tests, low memory protection,
\    error reporting expanded. Redesigned transient{
\ J: 2020=11=23 workiing on determinant. Added first{ second{ third{
\ K: 2020-11-24 added create(* for arbitrary cell size in a matrix
\    2020-12-12 clean through }transpose.

 anew xxx   decimal   5 set-precision

((  Save these for possible future use.
warning off  \ save and recover a float via return stack
: f>r     r> rp@ b/float - rp! rp@ f! >r ;
: r>f     r> rp@ f@ b/float rp@ + rp! >r ;  warning on
))

\  address manipulation and utilties
: data>base 3 cells- ; \ data address to base address
: base>data 3 cells+ ; \ base address to data adddress
: }}Dup  6 0 do 5 pick loop ;  \ copy 'data r c b/f' parameters

\ *** Components to create a matrix ***

here CONSTANT SafeMemory  \ all matrices must reside higher

: DataValid? ( data --- )  \ is 'data' in the users address space?
      SafeMemory < abort" not in users data space" ;

500 VALUE rcLimit \ cell limit per matrix

\ following used as synonyms for groups matrices
\ typical use:  first{ second{ third{  }+
0 VALUE	 first{  0 VALUE second{  0 VALUE third{

: RCvalid?  ( r c --- ) \ negative values not allowed
   2dup * rcLimit   > abort" r*c above limit:
       0< swap 0<  or abort" no negative r & c"  ;

: StorageValid? ( b/size --- ) \ cell size = b/float?
    dup 1 < swap b/float > or  abort" cell size out of range"  ;

: StackDepth? ( n --- ) \ must be less than or equal
    depth >=  abort" insufficient parameters"  ;

: Params>Data ( r c b/f data --- ) \ place at base address
   \ parammeters retrieved for a float matrix
      4   StackDepth?   \ enough parameters
    dup   DataValid?    \ in working ram?
    over  StorageValid? \ using floats?
                        \ RCValid? checked earlier
    data>base dup >r 2 cells+ ! swap r> 2!  ;

: {{ ( data --- data r c b/f ) \ get params, do error checks
    dup DataValid?  \ check data address
    dup data>base dup 2@ swap rot 2 cells+ @  \ data r c #
    >r  2dup RCvalid? r>
    dup StorageValid?   ;

: }Zeros  ( x{ --- ) \ fill matrix with zeros
    {{ * * erase  ;

: }Bytes  ( x{ --- byte_size ) \ storage bytes used
    {{ * * nip ;

\ *** transient{ matrix as workspace ***

2variable (transient) \ its data address and status boolean

: ClearTransient ( --- ) \ block from usage
   TRUE 0 (transient) 2! ;

: >Workspace  ( r c # X{  === )
  \ emplace parameters before the data address and clear
    dup >r  Params>Data  r> }Zeros   ;

: transient{  ( --- data ) \ access to transient{ data
   \ if status FALSE then leave data address
   (transient) 2@ swap abort" transient{ not configured"  ;

: SetupTransient{  ( r c --- ) \ a float workspace
   2dup RCValid?  \ following FALSE allows transient{ use
   here 100 + base>data  dup >r (transient) !
   b/float r> >Workspace
   FALSE (transient) cell+ ! ; \ correctly structured.

: create{*  ( r c size --- ) \ create a matrix with 'size' cell size
    >r  2dup RCValid?  r>
    create   ( size here )  here  base>data dup >r >Workspace
    r> }Bytes 3 cells+ allot \ allocate & clear storage
    DOES> base>data ( returns data address) ;

: create{  ( r c  --- ) \ create a float matrix
      b/float create{*  ;

\ *** Retrieving from a matrix  ***

: }}within? ( r1 c1 b/f r2 c2   --- r1 c1 b/f r2 c2 ) \ r2 c2 in range?
    \ Note: this is non-destructive
     over  0<  over  0< or   2 pick 6 pick >=  or
     over    5 pick  >=   or  abort" index out of range "  ;

: }}  ( data r1 c1 b/f r2 c2 --- cell-address )
   \ convert descriptor plus r c to its storage address
   }}within?  2 pick *  >r   * *  r> +  swap drop  +  ;

: }Dimensions ( x{ --- r c ) \ extract r c parameters
   {{  3 roll 2drop ( r c ) ;

: 2}rcConformal? ( x{ y(  ---  ) \ r's and c' equal
  }Dimensions rot }Dimensions D<> abort" rows/col not conformal" ;

: 3}rcConformal? ( X{ Y{ Z{ --- )
     over 2}rcConformal? 2}rcConformal?  ;

: }=size?  ( X{ Y{ --- ) \ equal storage size
    }Bytes swap }Bytes <> abort" not same size" ;

 : }square? ( X{ --- ) \ r and c are equal
   }Dimensions <> abort" not square"  ;

: }dot*conformal?  ( X Y{ Z{ ---  )
  dup }square?  }dimensions drop  >r
  }dimensions  2dup min >r
  rot }dimensions swap D<>
  r> r> <> or abort" not dot* conformal" ;

: }SetupLoop  ( X{ --- c r ) \ returns col then row
    }Dimensions swap  ;

\ *** words supporting matrices ***

\ fetch, display and store from 'data r c # r c'  input
: }}F@  }}  F@ ;   : }}F? }}F@ F. ;   : }}F! }} F! ;

: }}F+!   ( data r c # r c  f: F --- )
    }} F+! ; \ add FP1 item into FP2 at addr r c # r c

: }ToTransient ( X{ Y{  == )
   \ copy X{ to Transient{ resetting its r & c
   dup }Bytes  rot data>base rot data>base  rot 3 cells+ move ;

: }Copy ( X{ Y{  == )  \ copy X{ matrix to Y{ matrix
   2dup 2}rcConformal?   }ToTransient  ;

\ **** check the above ****???


: nc}Copy ( from{ to{  --- )
  \ non-conformal copy, reversing rc in the to{ matrix
  dup data>base dup 2@ swap rot 2! }copy ;

: }List   ( X{ ---  ) \ display a matrix
   dup }SetupLoop ( c r )  cr
   0  do cr ( row)  dup  0  do  ( column)
         over {{ j i }}F?   loop   loop  2drop ;

: }Fill  ( X{  ---  ) \ fill matrix with its element numbers
   dup }SetupLoop ( c r )
   0  do  ( data c )  dup 0  do
      dup j * i 1+ + s>f over {{ j  i  }}F!  loop loop  2drop ;

\ *** Data Entry ***
\ Used in the form: 2 2 create{ A{  A{ {[ 1 2 | 3 4 ]}

-111112 CONSTANT CHECKER \ mark stack during input

: StackCheck ( n f: pi --- ) \ verify stack has not changed
   CHECKER <> FPI F= not  or abort" input values?" ;

: ]}  ( add r1 c1 #  r2 c2 checker f: pi --- )
\ conclude input, stack check, rows & columns
   StackCheck
   0 <> swap 4 pick <> or abort" missing data" \ row error
   4drop ;  \ drop all params

 warning off
: |  \ accept c1 input float values filling in a matrix row
( add r1 c1 # r2 0 checker f: pi === add r1 c1 # r2+1 0 checker f: pi )
   StackCheck
   over  5 pick = \ error if at row limit or col not zero
   over  0 <>  or  abort" r c values error"
   }}Dup \ duplicate descriptor and  r c values
   }}    \ add r c #  0 0  addr
   4 pick  0 do [compile] F# dup F!  b/float + loop \ adding float values
   drop  swap 1+ swap   \ increment row number
   CHECKER FPI  ; \ for syntax check
 warning on

: {[  ( data --- data r c # 1 0 checker f: pi)
   \ begin matrix data entry by accepting first row
    {{ 0 0  CHECKER FPI  |  ;

\ *** Creating a character matrix ***

\ characer/byte matrix
1 28 1 create{* alphabet{c

: }AlphaFill  ( x{ --- )
  to first{   'A'  ( first alpha character )
  first{ }setuploop
      0 do  dup   0  cr do   ( 'a' c )
        first{ {{ j i }} ( char c address )
        2 pick swap c! swap 1+ swap
    loop loop  2drop ;

: }alphalist  ( X{ ---  ) \ display an alapha matrix
   to first{   first{ }SetupLoop ( c r )
   0  do cr ( row)  dup  0  do  ( column)
         first{ {{ j i }} c@ emit
    loop   loop  drop ;

cr cr .( Test of character matrix. )
alphabet{c }alphafill
alphabet{c }alphalist



 \ **** Major Matrix Words *******

: }= ( A{ B{  --- boolean ) \ matrix equality
  2dup 2}rcConformal?
  to second{ to first{ False
  first{ }SetupLoop ( from to c r )
    0 do ( row)  dup 0 do ( col)
            first{ {{  j i  }}F@
           second{ {{  j i  }}F@ F= not
           swap or swap
        loop loop drop not ;

: }+  ( X{ Y{ Z{ --- )  \ X+Y>Z )  \ matrix add
  3dup 3}rcConformal?
  to third{ to second{ to first{
  first{ }SetupLoop ( from to c r )
   0 do ( row)  dup 0 do ( col)
            first{ {{  j i  }}F@
           second{ {{  j i  }}F@  F+
            third{ {{  j i  }}F!
      loop loop  drop  ;

: }-  ( X{ Y{ Z{ --- )  \ x-y>Z )  \ matrix differnce
 3dup 3}rcConformal?
 to third{ to second{ to first{
  first{  }SetupLoop ( from to c r )
    0 do ( row)  dup 0 do ( col)
            first{ {{  j i  }}F@
           second{ {{  j i  }}F@  F-
            third{ {{  j i  }}F!
      loop loop  drop  ;


: }*  ( X{ Y{ Z{ --- )  \ x*y>Z ) \ matrix product
   3dup 3}rcConformal?
   to third{ to second{ to first{
   first{ }SetupLoop ( from to c r )
    0 do ( row)  dup 0 do ( col)
            first{ {{  j i  }}F@
           second{ {{  j i  }}F@  F*
            third{ {{  j i  }}F!
      loop loop  drop  ;

: }dot*    (  x{ y{ z{ --- )    \ matrix dot product
   3dup }dot*conformal?
   to third{ to second{ to first{
   second{ }setuploop swap ( c r)
   2dup  min dup dup * ( cols of product & its cell count )
   0 do ( over cells in product )
  over 0 do ( over columns of first, rows of second matrix)
   ( r c divisor )
   j over /mod  (  r c divisor remainder quotient )
      first{ {{  4 pick      i ( = quot i   ) }}F@
     second{ {{  i      6 pick ( = i rem    ) }}F@  F*
      third{ {{  4 pick 6 pick ( = quot rem ) }}F+!
      2drop ( the rem quot ) loop loop
     3drop  ; \ clean up r c and divisor

cr cr .( Test of }}dot* )

2 2 create{ A{  A{ {[ 1 2 | 3 4 ]}
2 2 create{ B{  B{ {[ 5 6 | 7 8 ]}
2 2 create{ C{
A{ B{ C{ }dot*  C{ }list

\ 19.000 22.000
\ 43.000 50.000

3 3 create{ G{ G{ {[ 1 2 3 | 4 5 6 | 7 8 9 ]}
3 3 create{ H{ H{ {[ 1 2 1 | 2 4 6 | 7 2 5 ]}
3 3 create{ I{
 G{ H{ I{ }dot*  I{ }list

\ 26.000 16.000 28.000
\ 56.000 40.000 64.000
\ 86.000 64.000 100.00  ok

2 3 create{ K{        K{ {[ 1 2 3 | 4 5 6       ]}
3 2 create{ L{        L{ {[ 1 2   | 2 4   | 7 2 ]}
2 2 create{ M{
      K{ L{ M{ }dot*  M{ }list

: }Transpose   ( X{ --- ) \ via transient{
   to first{    first{  }dimensions swap setuptransient{
      first{  }setuploop
   0 do ( row)  dup 0 do ( to c from )
         first{ {{  j i  }}F@   \ from
     transient{ {{  i j  }}F!   \ to
     loop loop  transient{ first{ nc}copy ;

cr cr .( Test of transpose)
 2 4  create{ J{
           J{  }fill
           J{  }list
           J{  }transpose
   transient{  }list
	   j{  }list

\ ***** Good to here *****

