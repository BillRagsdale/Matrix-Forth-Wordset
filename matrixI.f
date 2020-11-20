\  WFR  Nov. 17, 2020
\ Added number size to matrix definition.
\ Restructures the matrix parameter order.
\ Added matrix mulitiply
\ Major renaming using {{ and }}
\ Added }}dM*
\ G: matrices are located by their storage address
\ H: matrices now yield their address,  changes to { on data address, {{ on descriptor
\ i: revise conformal tests 2020-10-19

 anew xxx   decimal   5 set-precision

warning off  \ save and recover a float via return stack
: f>r     r> rp@ b/float - rp! rp@ f! >r ;
: r>f     r> rp@ f@ b/float rp@ + rp! >r ;  warning on

\ descriptor address manipulation
: data>base 3 cells- ; \ data address to base address
: base>data 3 cells+ ; \ base address to data adddress

\ : }Drop   4 0 do   drop loop ;  \ drop descriptor
\ : }Dup    4 0 do 3 pick loop ;  \ copy descriptor
\ : }Over   4 0 do 7 pick loop ;  \ copy descriptor
\ : }Swap   4 0 do 7 roll loop ;  \ swap descriptor

(( consider to make clear at loop end  data1 data2 c
: c}Drop   2drop ; \ after a loop, drop c and data1
: c2}Drop  3drop ; \ after a loop, drop c and 2 data
: c3}Drop  4drop ; \ after a loop, drop c and 3 data
))

: }}Dup  6 0 do 5 pick loop ;  \ copy descriptor & r c parameters

: scratch here 100 + ; ( === data)  \ workspace for transient{
\  used in the form  r c b/f scratch transient{

\ *** The most basic matrix words ***

: }Datavalid? ( data --- )  \ is in the users address space
      4495384 < abort" not in users data space" ;

: }GetParams ( data --- data r c # )
   dup }Datavalid? dup data>base dup 2@ swap rot 2 cells+ @  ;

: {{  ( data --- data r c # ) \ get matirx' parameters
   }GetParams dup b/float <>  abort" invalid descriptor"  ;

: Params>Data ( r c b/f data --- ) \ place at base address
   \ storage for floats and must be in ram
    over b/float <> over 4495384 < or abort" input incorrect"
    data>base dup >r 2 cells+ ! swap r> 2!  ;

: }Zeros  ( data --- ) \ fill matrix data with zeros
    {{ * * erase  ;

: transient{ ( r c # data === )
\ emplace parameters at base address and clear
\ used in the form: r c b/float scratch transient{
    dup >r  Params>Data  r> }Zeros   ;

: }Bytes ( data --- byte_size ) \ storage bytes used
    {{ * * nip ;

: create{  ( r c b/f --- ) \ full set required
    create  here  base>data transient{
    here  base>data }Bytes 3 cells+ allot \ allocate & clear storage
    DOES> base>data ( returns data address) ;

: create{F ( r c --- ) \ create an Float matrix
    B/FLOAT create{ ;

\ *** Error check utilities. All end in?  ***

: }Dimensions ( data --- r c ) \ show parameters
   {{  3 roll 2drop ( r c ) ;

: }}within? ( r1 c1 b/f r2 c2   --- r1 c1 b/f r2 c2 ) \ r2 c2 in range?
    \ Note: this is non-destructive
     over  0<  over  0< or   2 pick 6 pick >=  or
     over    5 pick  >=   or  abort" index out of range "  ;

: 2}rcConformal? ( data1 data2 --- ) \ r's and c' equal
    >r }Dimensions r> }Dimensions ( r c r c )
    rot <> >r <> r> or  abort" rows/cols not conformal"  ;

: 3}rcConformal? ( X{ Y{ Z{ --- )
     over 2}rcConformal? 2}rcConformal?  ;

: }=size?  ( x{{ y{{ --- ) \ equal storage size
    }Bytes swap }Bytes <> abort" not same size" ;

 : }square? ( data --- ) \ r and c are equal
   }Dimensions <> abort" not square"  ;

: }dot*conformal?  ( x{ y{ z{ ---  )
  dup }square?  }dimensions drop  >r
  }dimensions  2dup min >r
  rot }dimensions swap D<>
  r> r> <> or abort" not dot* conformal" ;

: }transient? ( data --- )
    here < abort" transient location error?" ;

: }SetupLoop  ( data --- c r )
\ will crash of dats is bad
    }Dimensions swap   ; \ from data address, setup loop C R

: n}SetupLoop  ( dataN data1 N --- N-thC N-thR )
 \ select N-th data address and setup loop C R
  pick  }SetupLoop ( C R )  ;

: dup{{  ( data0 --- data0   data0 r c bf )
        dup {{ ;    \ maintain data address and add descriptor

: over{{  ( data1 data0 --- data1 data0 data1 r c bf )
        over {{ ;    \ maintain data address and add descriptor

: pick{{  ( data1 data0 N --- data1 data0   dataN r c bf )
        pick {{ ;  \ maintain addressss, adding N-th descriptor

: }}  ( data r1 c1 b/f r2 c2 --- cell-address )
   \ convert descriptor plus r c to its storage address
   }}within?  2 pick *  >r   * *  r> +  swap drop  +  ;

\ fetch, display and store from {{ r c input
: }}F@  }}  F@ ;   : }}F? }}F@ F. ;   : }}F! }} F! ;

: }}F+!   ( data r c # r c  f: F --- )
    }} F+! ; \ add FP1 item into FP2 at addr r c # r c

: }Copy ( data2 data1  == )  \ copy data2 matrix to data1 matrix
   2dup }=size?   dup }Bytes
   rot data>base rot data>base  rot 3 cells+ move ;

: }List   ( data ---  ) \ display a matrix
   dup }SetupLoop ( c r )  cr
   0  do cr ( row)  dup  0  do  ( column)
         over{{ j i }}F?   loop   loop  2drop ;

: }Fill  ( data ) \ fill matrix with its element numbers
   dup }SetupLoop ( c r )
   0  do  ( data c )  dup 0  do
      dup j * i 1+ + s>f over{{ j  i  }}F!  loop loop  2drop ;

\ *** Data Entry ***

-111112 CONSTANT CHECKER \ mark stack during input

: StackCheck ( n f: pi --- ) \ verify stack has not changed
   CHECKER <> FPI F= not  or abort" values?" ;

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
   4 pick  0 do [compile] F# dup F!  b/float + loop
   drop  swap 1+ swap   \ increment row number
   CHECKER FPI  ;
warning on
: {[  ( data --- data r c # 1 0 checker f: pi)
   \ begin matrix data entry by accepting first row
    {{ 0 0  CHECKER FPI  |  ;

\ **** Major Matrix Words *******

: }= ( data1 data2 --- boolean ) \ compare matrices
  2dup 2}rcConformal?  False
  over }SetupLoop ( from to c r )
    0 do ( row)  dup 0 do ( col)
           3 pick{{  j i  }}F@
           2 pick{{  j i  }}F@
           swap F=  not or swap
        loop loop 2swap 3drop  ;

: }+  ( X{ Y{ Z{ --- )  \ X+Y>Z
 3dup 3}rcConformal?
  dup  }SetupLoop ( from to c r )
    0 do ( row)  dup 0 do ( col)
           3 pick{{  j i  }}F@
           2 pick{{  j i  }}F@  F+
             over{{  j i  }}F!
      loop loop  3drop  ;

: }-  ( X{ Y{ Z{ --- )  \ x-y>Z
 3dup 3}rcConformal?
 dup }SetupLoop ( from to c r )
    0 do ( row)  dup 0 do ( col)
           3 pick{{  j i  }}F@
           2 pick{{  j i  }}F@  F-
             over{{  j i  }}F!
      loop loop  3drop  ;

: }*  ( X{ Y{ Z{ --- )  \ x*y>Z
 3dup 3}rcConformal?
 dup }SetupLoop ( from to c r )
    0 do ( row)  dup 0 do ( col)
           3 pick{{  j i  }}F@
           2 pick{{  j i  }}F@  F*
             over{{  j i  }}F!
      loop loop  3drop  ;

: }dot*    (  x{ y{ z{ --- )  \ matrix dot product
   3dup }dot*conformal?
   2 n}setuploop swap ( c r)
   2dup  min dup dup * ( cols of product & its cell count )
   0 do ( over cells in product )
over 0 do ( over columns of lead, rows of second matrix)
   ( Xadd Yadd Zadd r c divisor )
   j over /mod  ( X Y Z r c divisor remainder quotient )
   7 pick{{  4 pick      i ( = quot i   ) }}F@
   6 pick{{  i      6 pick ( = i rem    ) }}F@  F*
   5 pick{{  4 pick 6 pick ( = quot rem ) }}F+!
   2drop ( rem quot ) loop loop
   3drop 3drop  ; \ clean up

cr cr .( Test of }}dot* )

2 2 create{F A{  A{ {[ 1 2 | 3 4 ]}
2 2 create{F B{  B{ {[ 5 6 | 7 8 ]}
2 2 create{F C{
a{ b{ c{ }dot*
c{ }list

\ 19.000 22.000
\ 43.000 50.000

3 3 create{F G{ G{ {[ 1 2 3 | 4 5 6 | 7 8 9 ]}
3 3 create{F H{ H{ {[ 1 2 1 | 2 4 6 | 7 2 5 ]}
3 3 create{F I{
 G{ H{ I{ }dot*  I{ }list

\ 26.000 16.000 28.000
\ 56.000 40.000 64.000
\ 86.000 64.000 100.00  ok

: }Transpose   ( data --- ) \ via scratch
   dup  {{ >r swap r> scratch transient{ drop
   scratch  over }SetupLoop
   0 do ( row)  dup 0 do ( to c from )
        2 pick{{  j i  }}F@   \ from
          over{{  i j  }}F!   \ to
    loop loop  drop swap }Copy ; \ copy scratch to original
cr .( Test of transpose)
g{ }list g{ }transpose
g{ }list

3 3 create{F  J{    J{ {[ 2 -5 3 | 0 7 -2 | -1 4 1 ]}

: }det   ( addr r c #  --- f: determinate )
   J{{ ( testing ) \ F0.0 F1.0 F1.0  \ final-value diffs sums
                        \ create proto at scratch
   }}over }}over }}Copy  \ duplicate matrix at scratch
   scratch cell+  dup @ 2 * swap !  \ double the column #
   }}dup * * over over + swap ( from to quantity ) move

 \   over 0 DO  \  over columns
 \ 2 pick 0 DO  \  over rows
 \   }}Dup  i j i +          }}F@ F*              \ addumulate sum
 \   }}Dup  i i 2 pick + i - }}F@ FSWAP F* FSWAP  \ accumulate diff
 \  loop
        \ rotate third item, add second, swap negate add
 \  FROT  F+ FSWAP FNEGATE F+ F1.0 F1.0 \ f1.0 f1.0 as place holders
 \ loop
 \  }}drop  fdrop fdrop

 ;   \ have determinate on F stack



