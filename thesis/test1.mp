% test1.mp -- do not edit, generated automatically by main.tex
input feynmp
require_RCS_revision "1.30";
beginchar(1, 65*1.0pt#, 50*1.0pt#, 0);
"feynmf: 1";
LaTeX_unitlength:=1.0pt;
subgraph (0, 0, w, h);
vinit;
pickup pencircle scaled thin;
vleft (__i1, __i2);
vright (__o1);
vconnect ("fermion,tension=1/3", __i1, __v1);
vconnect ("plain", __v1, __v2);
vconnect ("fermion", __v2, __v3);
vconnect ("plain", __v3, __i2);
vconnect ("photon", __v1, __o1);
vfreeze;
vdraw;
endsubgraph;
endchar;
% the end.
end.
endinput;
