% wilsonone.mp -- do not edit, generated automatically by main.tex
input feynmp
require_RCS_revision "1.30";
beginchar(1, 150*1.0pt#, 100*1.0pt#, 0);
"feynmf: 1";
LaTeX_unitlength:=1.0pt;
subgraph (0, 0, w, h);
vinit;
pickup pencircle scaled thin;
vleft (__i1);
vright (__o1, __o2);
vconnect ("fermion", __i1, __v1);
vconnect ("fermion", __i1, __v2);
vconnect ("plain", __v1, __o1);
vconnect ("plain", __v2, __o2);
vlabel ("$\lambda _1 n_1$", __v2);
vlabel ("$\lambda _2 n_2$", __v1);
vfreeze;
vconnect ("gluon", __v1, __v2);
vfreeze;
vdraw;
endsubgraph;
endchar;
beginchar(2, 150*1.0pt#, 100*1.0pt#, 0);
"feynmf: 2";
LaTeX_unitlength:=1.0pt;
subgraph (0, 0, w, h);
vinit;
pickup pencircle scaled thin;
vleft (__i1);
vright (__o1, __o2);
vconnect ("fermion", __i1, __o1);
vconnect ("plain,label=$\lambda _{1}n$,l.side=left", __i1, __v2);
vconnect ("fermion", __v2, __v3);
vconnect ("plain, label=$\lambda _2 n$,l.side=left", __v3, __o2);
vfreeze;
vconnect ("gluon,left,tension=0", __v2, __v3);
vfreeze;
vdraw;
endsubgraph;
endchar;
% the end.
end.
endinput;
