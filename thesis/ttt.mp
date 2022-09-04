% ttt.mp -- do not edit, generated automatically by main.tex
input feynmp
require_RCS_revision "1.30";
beginchar(1, 180*1.0pt#, 100*1.0pt#, 0);
"feynmf: 1";
LaTeX_unitlength:=1.0pt;
subgraph (0, 0, w, h);
vinit;
pickup pencircle scaled thin;
vleft (__i);
vblob ((.20w), __i);
vright (__o);
vconnect ("fermion,label=$p+k\rightarrow $, l.side=left", __i, __v3);
vconnect ("fermion,label=$p\rightarrow $, l.side=left", __v3, __o);
vfreeze;
vright (__v3, __o3);
vforce (((0.5w,0.5h)), __v3);
vforce (((0.5w,0h)), __o3);
vdot (__v3);
vconnect ("gluon, label=$k\downarrow $, l.side=left", __v3, __o3);
vfreeze;
vdraw;
endsubgraph;
endchar;
% the end.
end.
endinput;
