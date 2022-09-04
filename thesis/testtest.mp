% testtest.mp -- do not edit, generated automatically by main.tex
input feynmp
require_RCS_revision "1.30";
beginchar(1, 40*1.0pt#, 40*1.0pt#, 0);
"feynmf: 1";
LaTeX_unitlength:=1.0pt;
subgraph (0, 0, w, h);
vinit;
pickup pencircle scaled thin;
vleft (__i1);
vright (__o1);
vconnect ("plain,label=$k\rightarrow $,l.side=left", __i1, __o1);
vfreeze;
vdraw;
endsubgraph;
endchar;
beginchar(2, 40*1.0pt#, 40*1.0pt#, 0);
"feynmf: 2";
LaTeX_unitlength:=1.0pt;
subgraph (0, 0, w, h);
vinit;
pickup pencircle scaled thin;
vleft (__i1, __i2);
vright (__o1, __o2);
vconnect ("plain", __i1, __v, __o1);
vconnect ("plain", __i2, __v, __o2);
vfreeze;
vdraw;
endsubgraph;
endchar;
beginchar(3, 40*1.0pt#, 40*1.0pt#, 0);
"feynmf: 3";
LaTeX_unitlength:=1.0pt;
subgraph (0, 0, w, h);
vinit;
pickup pencircle scaled thin;
vleft (__i1);
vright (__o1);
vdot (__o1);
vconnect ("plain,label=$k\rightarrow $,l.side=left", __i1, __o1);
vfreeze;
vdraw;
endsubgraph;
endchar;
% the end.
end.
endinput;
