\  WFR  Jan. 30, 2021

\ A: 2021-01-30 A full rebuild dropping b/f after {{
\ B: 2021-02-01 Added back b/f after {{. OK data input tools.
\ C: 2021-02-04 Fixed transient{, All OK to row 380 }det.
\ D: 2021-02-05 Eliminated use of first{, second{ & third{
\ E: 2021-02-07 Rationalized structure. transient{ now re-entrant
\ F: 2021-02-19 Numeric input supports compillation.
\ G: 2021-02-22 Revised syttax checking. {{ & }}
\ H: 2021-03-01 Updated testing routines.
\ I: 2021-06-20 Refactored for sub-matrices
\ J: 2021-06-04 Adding column summation and Npick.

((
Should be ANS compliant except for NEXTWORD in {[int & {[comp.
NEXTWORD gets the next word to POCKET and reloads if necessary.
All numeric data is floating point. 120e or 12e1 for example.
Execution of aMatrix{ yields its initial data adddress.
Syntax summary:
matrix{ [optional r, c] }operation.  Operation on full or partial matrix.
matrix{ >rows< >cells<  }operation.  Specifies all rows and columns.
matrix{ {{ r c }} => address. Yield the address of the FP cell matrix{[r,c].
))

timer-reset \ Fun to watch elapsed time to compile and test.

anew MatrixSuite   decimal   5 set-precision

\ Stack access support
: Npick  create c, does> c@ pick ;
: Nroll  create c, does> c@ roll ;

 2 Npick  2pick   3 Npick  3pick   4 Npick   4pick   5 Npick 5pick
 6 Npick  6pick   7 Npick  7pick   8 Npick   8pick   9 Npick 9pick
10 Npick 10pick  11 Npick 11pick  12 Npick  12pick  13 Npick 13pick
14 Npick 14pick

 2 Nroll  2roll   3 Nroll  3roll   4 Nroll   4roll   5 Nroll 5roll
 6 Nroll  6roll   7 Nroll  7roll   8 Nroll   8roll   9 Nroll 9roll
10 Nroll 10roll  11 Nroll 11roll  12 Nroll  12roll  13 Nroll 13roll

1  VALUE TestLimit \ Execute tests numbered N and greater.
\ Set above parameter number (and above) for which tests will be done.
\ Set to zero to execute NO tests.

sys-warning-off
: *IF  ( n --- ) \ skip test to [then] if number is > ExecuteTests.
    TestLimit 0>  swap 1+  TestLimit > and dup [compile] [IF]
    if cr cr then  ;
\ Typical use: 3 *if cr cr .( 3  Test of transient{ and its support )
\   <some tests here>   [THEN]
sys-warning-on

((  Save these for possible future use.
warning off  \ save and recover a float via return stack
: f>r     r> rp@ b/float - rp! rp@ f! >r ;
: r>f     r> rp@ f@ b/float rp@ + rp! >r ;  warning on
))

\  *** Matrix address manipulation *** \

: data>base 3 cells- ; \ matrix data address to base address
: base>data 3 cells+ ; \ matrix base address to data adddress

\ *** Basic error checks  *** \

here CONSTANT SafeMemory  \ all matrices must reside higher

500 VALUE RClimit \ cell limit per matrix

: StackDepth? ( n --- ) \ must beless than or equal
    depth >=  abort" insufficient parameters"  ;

: SizeValid?  ( r1 c1 --- ) \ Number of cells within in limit:
    * RClimit > abort" Number of cells above limit." ;

: rcValid?    ( r1 c1 --- ) \ validate storage size, sign
   2dup SizeValid?  0< swap 0<  or abort" no negative r or c"  ;

: rcWithin? ( x{ # R C r1 c1 --- x{ # R C r1 c1 )    \ used within }} Note: this is non-destructive
    2dup rcValid?   over  4pick >=  over 4pick >= or
    abort" index out of range " ; \ error if r>=R or c>=C

\ *** Creation of matrices *** \

: create{*  ( size R C --- ) \ create a matrix with 'size' cell size
   3 Stackdepth?   2dup SizeValid?
   create  rot 1 max >r  ( have R C size size )
           r@ b/float > abort" cell size too large"
           r@ , 2dup , ,  r> * *
           here swap dup allot erase
   does>   base>data ;

: create{  ( R C  --- ) \ create a float matrix from its rows and columns.
     b/float rot rot create{*  ; \ b/float is 8 bytes

1 *if  .( 1  Testing StackDepth? ) 4 3 2 1   4 stackdepth? 4drop
    .( Stackdepth? is ok)
cr .( Testing SizeValid? )   45 10 SizeValid?
   .(     SizeValid?  is ok )
cr .( Testing rcValid? )     45 10 rcValid?
   .(       rcValid?    is ok)
cr .( Testing rcWithin? ) 8 5 6 4 5 rcwithin?  3drop 2drop
   .(      rcWiithin?  is ok)
cr cr .( Test of 4 3 create{  )
cr  .( see:       08 00 00 00 03 00 00 00  04 00 00 00 )
4 3 create{ A4x3{  A4x3{ data>base 3 cells dump  forget A4x3{
[then]

\ *** Access within matrices *** \

: {{ ( x{ --- addr size R C ) \ recover all parameters
\ it is important to keep the three checks independent.
  dup SafeMemory <                   abort" matrix syntax error"
  dup data>base @ dup b/float >      abort" matrix syntax error"
  over -2 cells+ 2@ 2dup * RClimit > abort" matrix syntax error" ;

: }}  ( addr size R C r1 c1 --- cell-address )
   \ convert descriptor to its storage address
   rcWithin? 4 pick * >r  * nip * + r> +
   dup SafeMemory < abort" parameter syntax error" ;

2 *if  .( 2  \ *** Testing {{ and  }} *** )
cr cr .( Test of A4x3{ {{ )   4 3 create{ A4x3{
A4x3{ {{ .( See  [4] xxxxxxx 8 4 3 ) cr 22 spaces .s 4drop
cr .( Test of 777e in A4x3{ {{ }} ) 777e A4x3{ 2 cells+ F!  A4x3{ {{ 0 1 }}
F@ f>s 777 =  [if]  .( got 777, correct ) [else] .( got zero, wrong) [then]
forget A4x3{
[then]

\ *** Support for transient{ work space. *** \
\ transient{ is a temporary workspace used for sub-matrix manipulation.
\ Transient workspace is 1050 bytes ahead of HERE in ram.
\ New definitions should not be added with an opem transient{.
\ You may have sequential opentransient{s. That is, they are reentrant.
\ Each opentransient{ must have a paired closetransient{.
\ They must be nested and closed in the reverse order of their creation.

0 VALUE (transient{) \ gives the data address of temporary matrix

1050 CONSTANT T{-offset \ offset above here of first transient{ allocation

: transient{ ( --- addr )  \ give the data-address for transient{
   (transient{) dup 0= abort" transient{ not initialized" ;

: OpenTransient{  ( R C --- ) \ add a transient F matrix above here.
   2 stackdepth?   2dup SizeValid?
   (transient{) dup 0=           ( R C 0 | data flag )
   if here T{-offset +            ( R C 0   new    )
       else dup {{ * * + 4 cells+ ( R C old new    )  then
   dup >r TO (transient{)           ( R C old >> r: new )
   r@ -4 cells+ !  b/float r@ -3 cells+ !  r> -2 cells+ 2!
   transient{ {{ * * erase ;

: CloseTransient{  ( --- ) \ return to prior transient{ allocation
   (transient{) dup 0= abort" transient{ not open"
    -4 cells+ @ to (transient{)  ;

3 *if  .( 3  Test of transient{ and its support )
3 4 opentransient{  transient{ {{ 1 1  }} F@ f>s closetransient{
cr [if]  .( got non-zero, wrong ) [else] .( got zero, correct) [then] [then]

\ *** Internal support words *** \

: }CellSize  ( x{ --- b/cell) \ get bytes per cell
    {{ 2drop nip ;

: }Cells      ( x{ --- #cells ) \ Matrix size in cells
    {{ *  nip nip ;

: }Bytes      ( x{ --- #bytes ) \ storage bytes used
    {{ * * nip  ;

: }Dimensions ( x{ --- R C ) \ extract the number of rows and columns
    {{ 2nip ;

: }SetupLoop  ( x{ --- C R ) \ used for nested DO-LOOPs.
    }Dimensions swap ;

: }square? ( x{ --- ) \ r and c are equal
   }Dimensions <> abort" not square"  ;

: 2}rcConformal? ( x{ y(  ---  ) \ both r's & c's equal
  }Dimensions rot }Dimensions D<> abort" rows/col not conformal" ;

: 3}rcConformal? ( x{ y{ z{ --- ) \ equal r's & c's for all three.
     over 2}rcConformal? 2}rcConformal?  ;

: 2RangeValid? ( r1 r2 r3 r4 --- )  \ used by }SubCopy
\  r1-r2 must equal r3-r4
    - >r - r> <> abort" row or column ranges differ" ;

: 3RangeValid? ( r1 r2 r3 r4 r5 r5 r6 --- ) \ used by (}SubX)
    2over 2RangeValid? 2RangeValid? ;

: >rows< ( x{ --- x{ r1 r2  ) \ note: non-destructive
   \ add all rows to matrix address being expanded
    dup }Dimensions     drop 1- 0 swap ;

: >cols< ( x{ r1 r2 --- x{ r1 r2 c1 c2 )
     \ add all columns to matrix address being expanded
    2pick }Dimensions   nip 1- 0 swap ;

4 *if  .( 4  Test of matrix parameters )
5 5 create{ A5x5{ 2 3 create{ A2x3{
cr .( Test of }square? ) A5x5{ }square? .(       Is ok )
cr .( Test of 3}rcConformal? ) A5x5{ dup dup 3}rcConformal? .( is OK)
cr .( Test of >rows< >cols<  See xxxxxxx 0 1 0 2 )
cr 23 spaces A2x3{ >rows< >cols< .s 3drop 2drop
forget A5x5{  [then]

\ *** Fetch, display and store from 'data f R C r c' input *** \

: }}@  }}  F@ ;   : }}? }}@ F. ;  : }}! }} F! ;

: }}+!   ( data f R C r c  F: f1 --- )
    }} F+! ; \ add f1 float into f2 currently at addr r c

5 *if  .( 5  Test of matrix access ) 2 3 create{ A2x3{
A2x3{ {{ 1 1 777e }}!  A2x3{ {{ 1 1 111e }}+! A2x3{ {{ 1 1 }}@
f>s 888 =  [if]  .( got 888. correct ) [else] .( wrong) [then]
forget A2x3{ [then]

\ *** Matrix based application words *** \

: }SubList cr  ( x{ r1 r2 c1 c2 --- ) \ Formatted list of a matrix.
    2 roll 1+             ( x{ r1 c1 c2 r2+1    )
    3 roll                ( x{    c1 c2 r2+1 r1 )
   do  dup 1+ 2 pick  do  ( x{    c1 c2         )
       2 pick {{ j i }}?  loop cr loop  3drop ;

: }list ( a{ --- ) \ list full matrix
  >rows< >cols<  }SubList  ;

: }rlist ( a{ r1 --- ) \ list a row
( want a{ r1 r1 c1 c2 )
   dup >cols< }SubList ;

cr cr .( Example Three, cell list 2nd row 2nd column )

: }CellList  ( a{ r1 c1  --- )
   swap dup rot dup }SubList ;

: }Zeros  ( x{ --- ) \ fill matrix with zeros
    {{ * *  erase  ;

: }Fill ( x{ --- ) \ Fill cells with their row & column numbers.
   dup  }SetupLoop 0 do  dup 0 do
   over {{ j 10 * i + s>f j i }}!  loop loop 2drop ;

6 *if   .( 6  Test of }zeros }list }fill ) 3 3 create{ A3x3{
A3x3{ }zeros A3x3{ }list
A3x3{ }Fill A3x3{ }list
A3x3{ 1 }rlist  A3x3{ 1 1 }CellList A3x3{ 1 2 1 2 }sublist
forget A3x3{  [then]

\ *** Data Entry {[ | ]} ***

: StackCheck ( n f: pi --- ) \ verify stack has not changed
    Fpi F= not  abort" data structure" ;

: ]}_  ( x{    F: Fpi n n n --- )
 \ conclude input converting stack value to floats into x{
   dup }dimensions   ( x{ R C   F: pi n . . . n )
   0    rot 1-   do  ( x{ C    over    rows R-1 up to 0 )
   0   over 1-   do  ( x{ C    over columns C-1 up to 0 )
       over {{ j i }}!
   -1 +loop  StackCheck  -1 +loop  2drop ;
sys-warning-off
: {[int ( x(  ---  ) \ interpreting, accept FP data until ]}
   Fpi     \ need leading StackCheck
   BEGIN        \ parsing word by word
      BL NEXTWORD 0= abort" need input text" count >float
           if ( leave FP on its stack ) false else
      s" |"  pocket count compare
             0= if Fpi ( place Fpi)     false else
      s" ]}" pocket count compare 0=
               if  ]}_                   true  else
             1 abort" data format? "     then  then  then
    UNTIL    ;

: {[comp    ( --- )  \ in : definition, convert and compile into FP
   compile Fpi  ( lead with StackCheck value )
   BEGIN        \ parsing word by word
      BL NEXTWORD 0= abort" need input text" count >float
            if [compile] Fliteral    false else
      s" |"  pocket count compare 0=
            if compile Fpi            false else
      s" ]}" pocket count compare 0=
            if compile  ]}_           true  else
             1 abort" data format? "  then  then  then
    UNTIL    ;
sys-warning-on
: {[  ( --- ) \ process remaing text until ]}
   state @     \ check STATE
   if {[comp  else ( with x{ ) {[int then ; immediate

7 *if  .( 7  Test of {[ | ]} executing and compiling )
4 4 create{ A4x4{
A4x4{ {[ 1  2  3  4 |  5  6  7  8  |
         9 10 11 12 | 13 14 15 16  ]}
A4x4{ }list A4x4{ }zeros
: test{[  {[ 1  2  3  4 |  5  6  7  8 |
             9 10 11 12 | 13 14 15 16 ]} ;
A4x4{ test{[   A4x4{ }list  forget A4x4{  [then]

\ *** Application words for region access within matrices *** \

: }Sub#Fill ( x{ r1 r2 c1 c2  F: f1 --- ) \ load float constant into area
\   4 pick }Validate?     ( done in {{      )
    2roll 1+  3roll     ( x{    c1 c2 r2+1 r1  )
   do  dup 1+ 2pick      ( x{    c1 c2 c2+1 c1  )
      do  ( x{ c1 c2 ) 2pick {{ j i  Fdup  }}!  loop loop
   3drop Fdrop ;

: }r#Fill ( x{ r1  F: F1 --- ) \ fill a row with a number
     dup >cols< }Sub#Fill ;

: }c#Fill  ( x{ c1 F: F11 --- ) \ fill a column with a number
   >r >rows< r> dup }Sub#Fill ;

: }#Fill  ( x{ F: f1 --- ) \ load full matrix with a number
    >rows< >cols< }Sub#Fill ;

: }SubNegate ( x{ r1 r2 c1 c2  --- ) \ negate a sub-matrix
\   4 pick }Validate?     ( done in {{      )
    2roll 1+             ( x{ r1 c1 c2 r2+1    )
    3 roll                ( x{    c1 c2 r2+1 r1 )
   do  dup 1+ 2pick  do  ( x{    c1 c2         )
       2pick {{ j i  }}  ( x{ c1 c2 cell-address )
       dup F@ Fnegate F!
      loop loop  3drop  ;

: }rNegate ( x{ r1 --- ) \ negate row of a matrix
      dup >cols< }SubNegate ;

: }cNegate ( x{ c1 --- ) \ negate columns of a
    >r >rows< r> dup  }SubNegate ;

: }negate ( x{ --- ) \ negate matrix contents
    >rows< >cols< }SubNegate ;

: }SubRandom ( x{ r1 r2 c1 c2  F: f1 --- )
 \ load random floats between zero and f1 into area
    2roll 1+   3 roll      ( x{    c1 c2 r2+1 r1 F: f1 )
   do  dup 1+ 2pick  do    ( x{    c1 c2         F: f1 )
       2pick {{ j i  Fdup f>s random s>f }}!
      loop loop  3drop Fdrop ;

: }Random ( x{  F: f1 --- ) \ fill matrix with random floats 0..f1
      >rows< >cols< }SubRandom ;

8 *if  .( 8  Test of matrix sub-areas. )
3 3 create{ A3x3{
A3x3{ }zeros A3x3{ 0 1 0 1 77e }Sub#Fill A3x3{ }list
A3x3{ }zeros A3x3{ 1 888e   }r#Fill      A3x3{ }list
A3x3{ }zeros A3x3{ 1 888e   }c#Fill      A3x3{ }list
A3x3{ }zeros A3x3{ 555e     }#fill       A3x3{ }list
             A3x3{  1 2 1 2 }SubNegate   A3x3{ }list  A3x3{ }zeros
        A3x3{ 1 2 1 2 1000e }SubRandom   A3x3{ }list
                A3x3{ 1000e }Random      A3x3{ }list  forget A3x3{ [then]

\ *** Copy submatrix, rows, columns and full matrices ***

: }RawCopy ( X{ Y{  == ) \ copy X to Y, possibly re-dimensioning
\ }RawCopy is very risky as limited error checks are made and it
\   directly writes into memory.
   over }cells      ( x{ y{ #xcells)
   over }cells      ( x{ y{ #xcells #ycells )
   <> abort" Size mismatch"
   over  data>base  ( x{ y{ Xaddr )
   swap  data>base  ( x{    Xaddr Yaddr )
   2roll }bytes    (       Xaddr Yaddr Xbytes )
   3 cells+  move ;

: }SubCopy  ( x{ r1 r2 c1 c2  y{ r3 r4 c3 c4 ---  )
 \ copy over x{ to y{ r1 r2  c1 c2 submatrix
 \ matrices can differ in dimensions!
   8pick 8pick 5pick 5pick  2RangeValid?
   6pick 6pick 3pick 3pick  2RangeValid?
   9pick }dimensions opentransient{
   9pick transient{ }RawCopy  \ transient{ is now the source matrix
                       ( x{ r1 r2 c1 c2  y{ r3 r4 c3 c4          )
   drop  nip           ( x{ r1 r2 c1 c2  y{ r3    c3             )
   3roll              ( x{ r1    c1 c2  y{ r3    c3 c2          )
   4pick - 1+         ( x{ r1    c1     y{ r3    c3  c2-c1+1    )
   5roll              ( x{ r1    c1     y{ r3    c3  c2-c1+1  r2  )
   6pick - 1+         ( x} r1    c1     y{ r3    c3  c-lim    r2-r1+1 )
                       ( x{ r1    c1     y{ r3    c3  c_lim r_lim  )
         0 do          ( x{ r1    c1     y{ r3    c3  c_lim   )
     dup 0 do          ( x{ r1    c1     y{ r3    c3  c_lim   )
     transient{        ( x{ r1 c1  y{ r3 c3  c_lim  trans{    )
             {{        ( x{ r1 c1  y{ r3 c3  c_lim  addr # R C)
                9pick j +  9pick i +
                       ( x{ r1 c1  y{ r3 c3 c_lim addr # R C r1+j c1+i )
  }}@      \ cell trans{ r1 c1 fetch
  3pick {{            ( x{ r1 c1  y{ r3 c3  c_lim  addr # R C)
  6pick j +  6pick i + ( x{ r1 c1  y{ r3 c3  c_lim  addr # R C r3+j c3+i )
  }}!      \ cell     Y{ r3 c3 store
  loop loop 4drop 3drop  closetransient{  ;

: }rCopy ( x{ y{ r1  --- )  \ a single row copy
     dup >cols< 4 roll ( x  r r c c y )
     3pick dup >cols<  }SubCopy ;

: }cCopy ( x{ y{ c1  --- )  \ a sigle column copy
  >r >rows< r> dup  ( x y r r c c )
  4roll >rows<  ( x r r c c y r r )
  3pick dup }SubCopy ;

: }rrCopy ( x{ r1 y{ r3  --- )  \ copy between rows
   >r >r dup >cols< r> r> dup >cols< }SubCopy ;

: }ccCopy ( x{ c1 y{ c4  --- )  \ copy between columns
   >r >r >r >rows< r> dup   r> >rows< r> dup }SubCopy ;

: }CellCopy ( x{ y{ r c --- ) \ copy a specific row/col cell
   2dup rot ( x{ y{ r r c c )
   2dup 2>r 2over 2>r 4roll 2r> 2r> }SubCopy ;

: }Copy   ( x{ y{ --- ) \ full matrix copy
   >r >rows< >cols< r> >rows< >cols< }SubCopy ;

9 *if  .( 9  Test of matrix copying words. )
5 5 create{ A5x5{ 5 5 create{ B5x5{ 6 6 create{ B6x6{
cr cr .( Test of }subcopy, }rcopy, }ccopy, }rrcopy, }}ccCopy, }CellCopy )
A5x5{ }fill A5x5{ 0 3 0 3  B6x6{ 2 5 2 5 }subcopy  B6x6{ }list
            B5x5{ }zeros A5x5{   B5x5{ 1 }rcopy    B5x5{ }list
            B5x5{ }zeros A5x5{   B5x5{ 1 }ccopy    B5x5{ }list
            B5x5{ }zeros A5x5{ 4 B5x5{ 1 }rrcopy   B5x5{ }list
            B5x5{ }zeros A5x5{ 4 B5x5{ 1 }cccopy   B5x5{ }list
            B5x5{ }zeros A5x5{   B5x5{   }Copy     B5x5{ }list
            B5x5{ }zeros A5x5{   B5x5{ 2 3 }Cellcopy B5x5{ }list
forget A5x5{   cr cr .( Stack check ) .s
[then]

1 [IF] \ *** Optional: Creating a character byte matrix ***

1 1 60 create{* alphabet{c

: }AlphaFill  ( x{ --- ) \ fill alpha matrix with  alphabet
   dup }SetupLoop 'A' swap    ( x{ c 'a' r )
   0 do   over  0  cr do      ( x{ c 'a' )
        dup 3pick {{ j i }}  ( x{ c 'a' 'a' addr )
        c! 1+
     loop loop  3drop ;

: }alphalist  ( X{ ---  ) \ display an alapha matrix
   dup }SetupLoop ( c r )  cr
   0  do cr ( row)  dup  0  do  ( column)
         over {{ j i }} c@ emit
    loop  loop  2drop  ;

10 *if  .( 10  Test of character matrix.)
 alphabet{c  }alphafill  alphabet{c  }alphalist
cr cr  .( Stack check: )  .s
[then]

\ *** Arithmetic by rows and columns ***

((  Test X range valid, test Y range valid, open transient matrix
    copy Z matrix to transient, form row range, form column range, copy X cell, copy Y cell
	execute operation, locate transient cell address, store result in transient cell
	loop over columns, loop over rows, copy transient matrix back to Z matrix
	close transient, clean up 11 itemns on stack. ))

: (}SubX) \ do math over submatrices
   ( x{ r1 r2 c1 c2 y{ r3 r4 c3 c4  z{ r5 r6 c5 c6 exec ---  )
 \ arithmetic over r1 r2 c1 c2 sub-matrix
  14pick 14pick 11pick 11pick 8pick 8pick 3RangeValid?
  12pick 12pick  9pick  9pick 6pick 6pick 3RangeValid?
   5pick dup }dimensions opentransient{
              transient{ }RawCopy  \ copy z{ to transient{
  ( x{ r1 r2 c1 c2  y{   r3  r4 c3 c4  z{  r5 r6 c5 c6 exec )
  13roll 12roll 2drop   8roll 7roll 2drop
    swap  2pick - 1+      3roll 4pick - 1+
    0 do dup 0 do ( by column, then row )
               ( x{ r1 c1 y{ r3 c3 z{ r5 c5 exec c_loop )
   10pick  {{ ( x{ r1 c1 y{ r3 c3 z{ r5 c5 exec c_loop addr # R C )
              13pick j +  13pick i + }}@
   7pick  {{ 10pick j +  10pick i + }}@  over execute
transient{ {{  7pick j +   7pick i + }}!
   loop loop  transient{ 5pick }copy   closetransient{
   4drop 4drop 3drop  ;

: }Sub+   \ add areas within matrices
  ( x{ r1 r2 c1 c2 y{ r3 r4 c3 c4 z{ r5 r6 c5 c6  --- )
        [']  F+   (}SubX)  ;

: }Sub-   \ subtract areas with matrices
  ( x{ r1 r2 c1 c2 y{ r3 r4 c3 c4 z{ r5 r6 c5 c6  --- )
        [']  F-   (}SubX)  ;

: }Sub.*  \ multiply areas within matrices
  ( x{ r1 r2 c1 c2 y{ r3 r4 c3 c4 z{ r5 r6 c5 c6  --- )
        [']  F*   (}SubX)  ;

: (F/)  ( F: f1 f2 --- f3 ) \ floating devision with zero test
    Fdup F0.0 F= abort" division by zero" F/ ;

: }Sub./  ( F: f1 f2 --- f3 ) \ sub-region floating div
  ( x{ r1 r2 c1 c2 y{ r3 r4 c3 c4 z{ r5 r6 c5 c6  --- )
       ['] (F/)  (}SubX) ;

11 *if  .( 11  Testing arithmentic by row and column)
5 5 create{ A5x5{   5 5 create{ B5x5{   5 5 create{ C5x5{
cr cr .( Testing (}SubX)
A5x5{ }fill  B5x5{ }fill  C5x5{ }zeros
A5x5{ 0 2 0 3  B5x5{ 0 2 0 3  C5x5{ 0 2 0 3 ' F+ (}SubX)
C5x5{ }list

cr cr .( Testing }Sub+, }Sub-, }Sub* amd }Sub/)
A5x5{ }fill   B5x5{ }fill    C5x5{ }zeros
A5x5{ 1 2 0 3  B5x5{ 2 3 0 3  C5x5{ 3 4 0 3 }Sub+
C5x5{ }list

A5x5{ }zeros   B5x5{ }fill    C5x5{ }zeros
A5x5{ 1 2 0 3  B5x5{ 2 3 0 3  C5x5{ 3 4 0 3 }Sub-
C5x5{ }list

A5x5{ }fill    B5x5{ }fill    C5x5{ }zeros
A5x5{ 1 2 0 3  B5x5{ 2 3 0 3  C5x5{ 3 4 0 3 }Sub.*
C5x5{ }list

B5x5{ }fill  B5x5{ }fill C5x5{ }zeros
a5x5{ 1 1 0 4 b5x5{ 2 2 0 4 C5x5{ 3 3 0 4 }Sub./
c5x5{ }list
forget A5x5{  cr cr .( Stack check: )  .s
[then]

\ ***** Full arithmetic by cell ******

: 3}expand ( x{ y{ z{ --- added rows and columns )
  >r >r >rows< >cols< r> >rows< >cols< r> >rows< >cols< ;

: }+  ( x{ y{ z{ --- ) \ matrix addition  x{ + y{ = z{
    3}expand }Sub+ ;

: }-  ( x{ y{ z{ --- ) \ matrix subtraction   x{ - y{ = z{
    3}expand  }Sub- ;

( Note: the following are the element by element multiply and divide. )

: }.*  ( x{ y{ z{ --- ) \ element multiplication  x{ .* y{ = z{
    3}expand }Sub.* ;

: }./  ( x{ y{ z{ --- ) \ element division  x{ ./ y{ = z{
    3}expand }Sub./ ;

12 *if  .( 12  Testing }+, }-, }* and }/ )
5 5 create{ A5x5{   5 5 create{ B5x5{   5 5 create{ C5x5{
A5x5{ }fill   B5x5{ }fill    C5x5{ }zeros
A5x5{ B5x5{ C5x5{  }+  C5x5{ }list

A5x5{ }zeros   B5x5{ }fill    C5x5{ }zeros
A5x5{  B5x5{ C5x5{ }-  C5x5{ }list

A5x5{ }fill    B5x5{ }fill    C5x5{ }zeros
A5x5{  B5x5{ C5x5{ }.*  C5x5{ }list

A5x5{ }fill    B5x5{ 10e0 }#fill    C5x5{ }zeros
A5x5{  B5x5{ C5x5{ }./  C5x5{ }list
forget A5x5{   cr cr .( Stack check: )  .s
[then]

cr cr .( Expanding matrix math words. )

: (}r)  ( a{ b{ c{ r1 --- ) \  Internal for specfic row sum.
  >r  rot r@ dup >cols<   ( b{ c{ a{ r1 r1 c1 c2 )
  6 roll r@ dup >cols<   ( c{ a{ r1 r1 c1 c2 b{ r1 r1 c1 c2 )
  10roll r> dup >cols< ; ( a{ r r c c b{ r r c c c{ r r c c )

: }r+  ( a{ b{ c{ r1 --- ) \  sum specified row across three matrices.
  (}r)  ['] F+ (}SubX) ;

: }r-  ( a{ b{ c{ r1 --- ) \  product specified row across three matrices.
  (}r)  ['] F- (}SubX) ;

: }r.*  ( a{ b{ c{ r1 --- ) \  sum specified row across three matrices.
  (}r)  ['] F* (}SubX) ;

: }r./  ( a{ b{ c{ r1 --- ) \  sum specified row across three matrices.
  (}r)  ['] (F/)  (}SubX) ;

\ *** And similarly for minus, times and divide.
13 *if  .( 13  Test of }r+, }r-, }r*, }r/  ) cr
3 3 create{ a{  a{ 7e }#fill a{ }list
a{ a{ a{ 1 }r+  a{ }list a{ 7e }#fill
a{ a{ a{ 1 }r-  a{ }list a{ 7e }#fill
a{ a{ a{ 1 }r.*  a{ }list
a{ a{ a{ 1 }r./  a{ }list  forget a{   [THEN]

: (}rrr)  ( a{ r1 b{ r2 c{ r3 --- 15v ) \  Expand rrr to 15 values.
                2>r  2swap  ( b{ r2 a{ r1 )
                dup >cols<  ( b{ r2 a{ r1 r1 c1 c2 )
  6 roll 6 roll dup >cols<  ( a{ r1 r1 c1 c2 b{ r2 r2 c1 c2 )
            2r> dup >cols<  ;

: }rrr+  ( a{ r1 b{ r2 c{ r3 --- ) \  sum cells in the specified rows .
    (}rrr)   ['] F+  (}SubX) ;

14 *if  .( 14  .( Check of }rrr+, }rrr-, }rrr*, }rrr/  ) cr
3 3 create{ a{  3 3 create{ c{  a{ 7e }#fill
a{ 0  a{ 1 c{ 2 }rrr+  c{ }list c{ }zeros
forget a{   [THEN]

: (}c)    ( a{ b{ c{  c1 --- ) \ Exapd }c to 15 parameters
  >r rot  >rows< r@ dup    ( b{ c{ a{ r1 r2 c1 c1 )
  6  roll >rows< r@ dup    ( c{ a{ r1 r2 c1 c1 b{ r1 r2 c1 c1 )
  10 roll >rows< r> dup ;  ( a{ r r c c b{ r r c c c{ r r c c )

: }c+    ( a{ b{ c{  c1 --- ) \ Sum cells in the common columns
    (}c)  ['] F+ (}SubX) ;

15 *if  .( 15 Check of }c+, }c-, }c*, }c/  ) cr
3 3 create{ a{  3 3 create{ c{   a{ 7e }#fill
a{ a{ c{ 1 }c+  c{ }list c{ }zeros
forget a{  [THEN]

: (}ccc)  ( a{ c1 b{ c2 c{ c3 --- ) \  Expand }cccx parameters.
  2>r  2swap                     ( b{ c2 a{ c1 )
  >r >rows< r> dup               ( b{ c2 a{ r1 r2 c1 c1 )
  6roll 6roll >r >rows< r> dup ( a{ r1 r2 c1 c1 b{ r1 r2 c2 c2 )
  2r> >r >rows< r> dup ; ( a{ r1 r1 c1 c2 b{ r2 r2 c1 c2 c{ r1 r2 c3 c3 )

: }ccc+  ( a{ c1 b{ c2 c{ c3 --- ) \  sum specified columns.
  (}ccc)  [']  F+  (}SubX) ;


16 *if  .(  16 Check of }ccc+, }ccc-, }ccc*, }ccc/  ) cr
3 3 create{ a{  3 3 create{ c{  a{ 7e }#fill
a{ 0 a{ 1 c{ 2 }ccc+  c{ }list c{ }zeros
forget a{  [then]

: }com+ ( a{ b{ c{ r1 r2 c1 c2 --- ) \ sums cells over the common area.
  2over 2over 2>r 2>r ( a{ b{ c{ r1 r2 c1 c2 )
  2over 2over 2>r 2>r ( a{ b{ c{ r1 r2 c1 c2 )
  5roll 2r> 2r>  ( a{ c{ r1 r2 c1 c1 b{ r1 r2 c1 c2 )
  9roll 2r> 2r>  }sub+  ;

: }rc+  ( a{ b{ c{ r1 c1 --- ) \  sum specified a{ b{ cells to c{
    swap dup rot dup  }com+ ;


17 *if  .( 17  Testing }r+, }rrr+, }c+, }ccc+  }rc+,  }com+, }sub+ )
cr 3 3 create{ a{   3 3 create{ b{  3 3 create{ c{
a{ }fill  b{ }fill  c{ }zeros
a{ }list   \  b{ }list c{ }list
a{ b{ c{ 1       }r+     c{ }list  c{ }zeros
a{ 0  b{ 1 c{ 2  }rrr+   c{ }list  c{ }zeros
a{ b{ c{ 1       }c+     c{ }list  c{ }zeros
a{ 0 b{ 1 c{ 2   }ccc+   c{ }list  c{ }zeros
a{ b{ c{ 1 1     }rc+    c{ }list  c{ }zeros
a{ b{ c{ 1 2 1 2 }com+   c{ }list  c{ }zeros
a{ 0 1 0 1 b{ 0 1 1 2 c{ 1 2 1 2 }sub+ c{ }list  c{ }zeros
\ forget a{
[then]

cr cr .( Summation of columns a{ int o b{ )

: }(SubCsum) ( a{ r1 r2 c1 c2 b{ r3 c3 ) \  Sub Column Summation
\ sum columns of a{ sub-matrix into row of c{ with r3 c3 top left
    5pick 7pick <  4pick 6pick < or abort" r/c range error"
    2pick dup }dimensions opentransient{
              transient{ }RawCopy  \ copy b{ to transient{
    5roll  6pick   \ a{ r1 c1 c2 b{ r2 c3 r2 r1  row span
    - 1+           \ a{ r1 c1 c2 b{ r2 c3 r2-r1+1
    4roll  5pick   \ a{ r1 c1  b{ r2 c3 r2-r1+1 c2 c1
    - 1+   0       \ a{ r1 c1  b{ r2 c3 r2-r1+1 c2-c1+1 0
   do              \ a{ r1 c1  b{ r2 c3 r2-r1+1
    F0.0   dup 0   \ a{ r1 c1 b{ r3 c3 r2-r1+1   F: sum
   do              \ a{ r1 c1 b{ r3 c3 r2-r1+1
     6pick {{              \ a{ r1 c1 b{ r3 c3 r2-r1+1  add r R C
     9pick i +  9pick j +  \ a{ r1 c1 b{ r3 c3 r2-r1+1  add r R C r1+i c1+j
     }}@   F+              \ a{ r1 c1 b{ r3 c3 r2-r1+1
   loop  \ down one row
     transient{ {{ 6pick 6pick   \ a{ r1 c1 b{ r3 c3 r2-r1+1 add f R C r3 c3
     i + }}!                     \ a{ r1 c1 b{ r3 c3 r2-r1+1
   loop
     transient{ 4pick }RawCopy   closetransient{
     4drop 3drop ;

: }ccSum  ( a{ c1 c2 b{ r3 ) \ sum columns c1-c2 into r3 of b{
   2>r 2>r >rows< 2r> 2r>    \ a{ r1 r2 c1 c2 b{ r3
   3pick }(SubCsum) ;

: }cSum  ( a{ c1 b{ r3 ) \ sum column c1 into row r3 of b{
   2>r dup 2r> }ccSum ;

18 *if  .( 18 Check of }{subCsum. }cSum, ]ccSum, }cccSum )

4 4 create{ 4x4A{   4X4a{ }fill
2 4 create{ 2x4B{   2x4B{ }zeros
4x4A{ }list
4x4A{ 0 3 0 3 2x4B{ 1 0 }(subCsum)
2x4B{ }list   2x4B{ }zeros
4x4A{ 0 3 0 1 2x4B{ 1 2 }(subCsum)
2x4B{ }list   2x4B{ }zeros
4x4A{ 0 2 2x4B{ 0 }ccSum
2x4B{  }list   2x4B{ }zeros
4x4A{ 2 2x4B{ 1 }cSum
2x4B{  }list    forget 4x4A{
[then]

\   General form across columns
\ a{ c1 c2 b{ r1  }ccN     \ over columns, count of each column's cells
\ a{ c1 c2 b{ r1  }ccMean  \ over columns, mean of each column
\ a{ c1 c2 b{ r1  }ccSum   \ over columns, total of each column
\ a{ c1 c2 b{ r1  }cc^2Sum \ over columns, sum of squares of each cell

: }ccN  ( a{ c1 c2 b{ r3  --- )   \ over the column range, the number of rows
   2pick 4pick <  abort" columns error"
   4pick }dimensions drop S>F  \ row count
   2roll 1+ 3roll  \ a{ b{ r3  c2+1 c1
   do  over {{  4pick  i Fdup }}!  loop
   3drop Fdrop ;

: }cN  ( a{ c1 b{ r3  --- )   \ the column's number of rows
  2>r dup 2r>  }ccN ;

: }Sub=   \ sub-matrix equality
 ( x{ r1 r2 c1 c2  y{ r3 r4 c3 c4  --- boolean )
   8pick 8pick 5pick 5pick  2RangeValid?
   6pick 6pick 3pick 3pick  2RangeValid?
                       ( x{ r1 r2 c1 c2  y{ r3 r4 c3 c4          )
   drop  nip           ( x{ r1 r2 c1 c2  y{ r3    c3             )
   3roll              ( x{ r1    c1     y{ r3    c3 c2          )
   4pick - 1+         ( x{ r1    c1     y{ r3    c3  c2-c1+1 )
   true \ for tracking result flag
   6roll              ( x{ r1    c1     y{ r3    c3  c2-c1+1 flag r2      )
   7pick - 1+         ( x} r1    c1     y{ r3    c3  c-lim   flag r2-r1+1 )
                                              ( x{ r2 c1 y{ r3 c3 c_lim flag r_lim   )
         0 do                                 ( x{ r1 c1 y{ r3 c3 c_lim flag         )
     over 0 do                                ( x{ r1 c1 y{ r3 c3 c_lim flag   )
     7pick {{ 10pick j + 10pick i + }}@    ( x{ r1 c1 y{ r2 c2 c_lim flag   F: F@ )
     4pick {{  7pick j +  7pick i + }}@ F= ( x{ r1 c1 y r2 c2 c_lim flag flag )
       and                                    ( x{ r1 c1 y{ r2 c2 c_lim flag )
     loop loop >r 4drop 3drop r> ;

\  *** Maybe change syntax to x{ y{ r1 }r=  ??

: }r=   ( x{ y{ r1  --- boolean ) \ equality over a row
      dup >cols<        2dup 2>r 2over 2>r  4roll 2r> 2r> }Sub=  ;

: }c=   ( x{ y{ c1  --- boolean ) \ full matrix equality
      >r >rows< r> dup  2dup 2>r 2over 2>r  4roll 2r> 2r>  }Sub=  ;

: }=   ( x{ y{  --- boolean ) \ full matrix equality
      >r >rows< >cols< r> >rows< >cols< }Sub=  ;

19 *if  .( 19  Test of matrix and sub-matrix equality )
5 5 create{ A5x5{  5 5 create{ B5x5{
A5x5{ }fill  B5x5{ }fill
A5x5{ {{ 0 0 }} 123e0 F!  cr  A5x5{ }list  cr
A5x5{ 0 3 0 3 B5x5{ 0 3 0 3 }Sub=   \ is false
cr  .( Test of }Sub= )  [if] .( got true ) [else] .( got false, correct) [then]

A5x5{ 1 4 1 4 B5x5{ 1 4 1 4 }Sub=   \ is is true
cr  .( Test of }Sub= )  [if] .( got true,  correct ) [else] .( got false) [then]

A5x5{ B5x5{ }=                      \ is is false
cr .( Test of }= ) [if]  .( got true ) [else] .(    got false, correct) [then]

A5x5{ }zeros B5x5{ }zeros B5x5{ {{ 2 2 7e }}!  B5x5{ }list
A5x5{ B5x5{ 2 }r=
cr .( Test of }r= ) [if] .( got true, wrong) [else] .( got false correct) [then]

A5x5{  B5x5{ 2 }c=
cr .( Test of }c= ) [if] .( got true, wrong) [else] .( got false correct) [then]
forget A5x5{
[then]

cr cr .elapsed

cr cr .( This is end of operational code.
\s

 ********* UNDER DEVELOPMENT, HERE ONWARD **********

((  wamt mean, stdev, variance,
row 0 is row count  N
row 1 is column sum  sum(column)
row 2 is column mean sum(column)/N
transient column 0 is corrected means   x - mean(x)
transient column 1 is square of corrected means
row 3 is sum of corrected means / n  is variance

))

cr cr .( Testing Descriptive Statistics
4 4 create{ 4x4A{   4x4A{ 3.e0 }#fill   4x4A{ }list
4 4 create{ 4x4B{
4x4A{ 1 2 4x4B{ 2  }ccN   4x4B{ }list 4x4B{ }zeros
4x4A{ 3   4x4B{ 1  }cN    4x4B{ }list 4x4B{ }zeros


: }*conformal?                    ( X{ Y{ Z{ ---  ) \
   2roll }dimensions                (    y{ z{ xR xC )
   3pick }dimensions                (    y{ z{ xR xC yR yC )
   swap d= not                       (    Y{ Z{ flag )
   2roll }dimensions nip            (       z{ flag CofY )
   2pick }dimensions nip            (       z{ flag CofY CofZ )
   = not or                          (       z{ flag )
   swap }dimensions <> or  abort" Not dot* conformal" ;

( Note: }* is the linear agebra matrix multiply.  For element by element, use }.*  )

: }*    (  x{ y{ z{ --- )    \ matrix product
   3dup }*conformal?  dup }zeros
   over }setuploop swap     ( x{ y{ z{ r c )
   2dup  min dup dup *      ( x{ y{ z{ r c min z{cellcount )
   ( cols of product & its cell count )
       0 do     ( x{ y{ z{ r min )  ( over cells in product )
  over 0 do     ( x{ y{ z{ r min )
      ( over columns of x{, rows of y{  ( r c divisor )
      j over /mod  ( x{ y{ z{ r c divisor remainder quotient )
      7pick {{  4pick      i ( = quot i   ) }}@
      6pick {{  i     6pick ( = i rem    ) }}@  F*
      5pick {{  4pick 6pick ( = quot rem ) }}+!
      2drop ( the rem quot ) loop loop
      3drop 3drop  ; \ clean up r c and divisor

14 *IF  cr cr .( 14 Test of }*conformal? )
4 3 create{ A4x3{  3 4 create{ B3x4{  4 4 create{ C4x4{
2 2 create{ A2x2{  2 2 create{ B2x2{  2 2 create{ C2x2{
3 3 create{ A3x3{  3 3 create{ B3x3{  3 3 create{ C3x3{
2 3 create{ A2x3{  3 2 create{ B3x2{
A4x3{ B3x4{ C4x4{ }*conformal?
cr .(            }*conformal? is OK )

cr cr .( Test #1 of }* )  \ tests OK
A2x2{ {[ 1 2 | 3 4 ]}
B2x2{ {[ 5 6 | 7 8 ]}
A2x2{ B2x2{ C2x2{ }*  C2x2{ }list    \ tests OK

\ 19.000 22.000
\ 43.000 50.000

cr cr .( Test #2 of }* )
A3x3{ {[ 1 2 3 | 4 5 6 | 7 8 9 ]}
B3x3{ {[ 1 2 1 | 2 4 6 | 7 2 5 ]}
A3x3{ B3x3{ C3x3{ }*  C3x3{ }list

\ 26.000 16.000 28.000
\ 56.000 40.000 64.000
\ 86.000 64.000 100.00  ok

cr cr .( Test #3 of }* )

A2x3{ {[ 1 2 3 | 4 5 6       ]}
B3x2{ {[ 1 2   | 3 4   | 5 6 ]}
A2x3{ B3x2{ C2x2{ }*  C2x2{ }list
\ 7.0000 10.000
\ 19.000 28.000
forget A4x3{ ( And all matrices after.)    [then]

\ *** Major applications words *** \

: }Transpose   ( X{ --- ) \ via transient{
      dup }Dimensions swap OpenTransient{
      dup }setuploop       ( x{ c r )
   0 do ( row)  dup 0 do   ( x{ c )
              over  {{  j i  }}@   \ from
         transient{ {{  i j  }}!   \ to
     loop loop drop ( c )
     transient{ swap ( tr{ x{ ) }RawCopy CloseTransient{  ;

15 *if  cr cr .( 15  Test of transpose)
4 3 create{ A4x3{ A4x3{ }fill  A4x3{ }list
A4x3{ }transpose A4x3{ }list  forget A4x3{     [then]

\ *** Determinant ***

: }replicate2   ( x{ ---  ) \ replicate two x{ in transient
   dup  }setuploop  ( c r )
   0 DO dup 0 DO
         over  {{ j i  }}@   Fdup
    transient{ {{ j i  }}!
    transient{ {{ j i 6pick + }}!
     loop loop  2drop  ;

: }det ( x{ --- f: determinate)
       dup }Square?       \ just checking
       dup  }Dimensions ( r c ) dup + ( r c*2 ) OpenTransient{
       dup  }replicate2   \ duplicate copy at transient{
    0e0  1e0  1e0         \ final value, partial sum, partial difference
       dup }SetupLoop     \  cols, rows of original
   0 DO  dup  0 DO   ( x{ c )
       transient{ {{ i  j i +              }}@  Fswap  F*  Fswap  \ partial sum
       transient{ {{ i  5pick  1- i - j +  }}@  F*    \ parital difference
    loop F- F+  1e0 1e0   \ final, partial sum, partial difference
    loop F2drop 2drop ( c value )
    CloseTransient{  ;

16 *if  cr cr .( 16  Test of }det)  3 3 create{ A3x3{   3 3 create{ B3x3{
A3x3{ {[ 2 -5 3 |  0 7 -2 | -1 4 1 ]}
\ det = 41
B3x3{ {[ -3 9 5 | -4 0 1 |  6 3 0 ]}
\ det = 3
A3x3{ }det cr cr f.  .( expect 41.000)
B3x3{ }det cr cr f.  .( expect 3.0000)  forget A3x3{
[then]

: }rExchange ( x{ r1 r2 --- )
  \ in matrix x{ exchange rows r1 & r2 all columns )
  \ copy x{ to transient{, then copy rows back.
   2pick dup }dimensions opentransient{ transient{  }copy
   transient{  2pick  4pick 3pick  }rCopy
   transient{   swap  3roll 3roll  }rCopy  closetransient{ ;

17 *if  cr cr .( 17  Test of }Rexchange )  5 5 create{ A5x5{
A5x5{ }fill  A5x5{ 0 4  }Rexchange  A5x5{ }list forget A5x5{
cr cr .( Stack check  ) .s
[then]

\ *********** Eschelon for linear algebra ************

cr cr .(  WORK HERE FOR A FUNCTIONAL WORD SET )

: }rc#= ( x{ r c F: f1 --- boolean )
  \ tests the r c cell for equality to f1
    >r >r {{ r> r> }}@ F= ;

: }RxCreduce ( x{ r1 c1  --- )
 \ divide row r by value in column c. Sets [r2 c2] to 1.
  2pick {{ ( x{ r1 c1 addr # r2 c2 )
  5pick 5pick ( x{ r1 c1 addr # r2 c2 r1 c1 ) }}@ ( divisor)
  2pick }dimensions nip 0 do ( {x r1 c1 )
  2pick {{ ( x{ r1 c1 add # R C )
            5pick i }} dup F@ Fover (F/) F!
      loop 3drop Fdrop ;

: }Reduce   ( x{ r1 c1 --- )
   \ Reduce x{ over row1 through R at column1
   2pick }dimensions drop ( x{ r1 c1 R )
   2roll            ( x1 c2 R r1 ) do \ over rows
   over  i  2pick   ( x{ c2 x{ i col )  }RxcReduce
      loop  2drop ;

: }R-Rreplace ( x{ r1 r2 --- )
  \ row r1 - row r2 replaces row r2
  2pick }dimensions nip 0  do \ over all columns
  2pick {{ ( x{ r1 r2 addr # R C ) 5pick i }}@
  2pick {{                         4pick i }} dup F@
    Fswap F- F!  \ store in r2 saved address
    loop  3drop  ;

: }Replace ( x{ r1 --- )
\  row replacement at r1  down to all lower rows )
   over }dimensions drop  ( x{ r1  R )
   over - over do
   2dup i 1+ }R-Rreplace loop  2drop ;

18 *IF  cr cr .( 18  Test of linear equation word set)
3 4 create{ A3x4{  5 5 create{ A5x5{
: reload A3x4{ {[ 2 -3 2 14 | 1 2 3 6 | 3 1 -1 -2 ]} ;

cr cr .( Test of }RxCreduce [1 0] )
reload   A3x4{ }list
A3x4{ 0 0 }RxCreduce          A3x4{ }list

cr cr .( Test of }Reduce)
reload   A3x4{ 0 0 }Reduce    A3x4{ }list

cr cr .( Test of }Replace )
reload A3x4{ 0 0 }R-Rreplace  A3x4{ }list

reload A3x4{ 1 0 }reduce  A3x4{ 1 }Replace A3x4{ }list cr
A5x5{ }fill  A5x5{ 1 0 10.0e }rc#=
cr .( Test of }rc#= ) [if]  .( got true, correct ) [else] .( got false wrong) [then]
cr cr .( Stack check  ) .s  forget A3x4{
[then]

cr cr .(  HAS THE WORD SET BEEN REFINED? )


cr cr .elapsed ( 64 msec to compile with no testing.)


