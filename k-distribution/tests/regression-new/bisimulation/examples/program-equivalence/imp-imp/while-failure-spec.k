module SYNC-POINT-CANDIDATES

imports MAIN

syntax Id ::= "n" [token]

rule
<s> ~ a => ^ b ... </s>
<ak> while (_) _ </ak>
[trusted]

rule
<s> ~ a => ^ b ... </s>
<ak> .K </ak>
[trusted]

rule
<s> ~ b => ^ a ... </s>
<bk> while (_) _ </bk>
<bstate> n |-> N:Int _ </bstate>
requires N modInt 2 ==Int 0
[trusted]

rule
<s> ~ b => ^ a ... </s>
<bk> .K </bk>
[trusted]

endmodule

module WHILE-FAILURE-SPEC

imports MAIN
imports SYNC-POINT-CANDIDATES

// wrong loop unrolling (`<` should have been used instead of `<=`)
rule
<k> check => . </k>
<s> ^ a ... </s>
<aT>
  <ak>
    while (n <= 10) {
      n = n + 2;
    }
  =>
    _
  </ak>
  <astate> n |-> N1:Int _ => _ </astate>
</aT>
<bT>
  <bk>
    while (n <= 10) {
      n = n + 1;
    }
  =>
    _
  </bk>
  <bstate> n |-> N2:Int _ => _ </bstate>
</bT>
requires N1 ==Int N2
 andBool N1 modInt 2 ==Int 0

rule
<k> check => . </k>
<s> ^ a ... </s>
<aT>
  <ak> .K </ak>
  <astate> n |-> N1:Int _ </astate>
</aT>
<bT>
  <bk> .K </bk>
  <bstate> n |-> N2:Int _ </bstate>
</bT>
requires N1 ==Int N2
[trusted]

endmodule
