README for DarkSUSY documentation.

The documentation to DarkSUSY is distributed together with the code and
organized as follows. The main documentation is a latex file created from the
documentation of the various parts of DarkSUSY. This latex file is then used
as a source both for the pdf manual and for an html version (using latex2html
to generate the html output).

The latex file is created by a script 'scr/headers2tex.pl' and generates the
master tex file Manual.tex in docs/ (in the DarkSUSY root). This script
performs the following tasks:

* It generates a template file, adding the definitions and macros in 
  src/docs/headers

* It then adds front page material, which are all files src/docs/H*.tex
  (added in alphabetical order)

* It then generates the table of contents

* It then adds introductory chapters, i.e. all files src/docs/I*.tex

* The main task is then to add documentation for the various subdirectories
  of DarkSUSY. For each subdirectory in src (really those defined in
  headers2tex.pl), 
  - it first adds the files *.tex (in alphabetical order)
  - it then reades the headers of the Fortran files (everything above
    the function / subroutine declaration) and add these.
    If the keywords 'BeginTex' and 'EndTex' occurs, the text in between
    (apart from Fortran comment tags in the first column(s)) is 
    added as latex code. Other text in the headers is added as verbatim
    text

* After all subdirectories are gone through, appendices are added from
  src/docs/A*.tex (currently none)

* The end material (acknowledgements, references) are then added from
  src/docs/E*.tex

Eventual figues should be refered to as fig/<filename> and be put in
src/docs/fig/. When headers2tex.pl has completed, a file Manual.tex is created
in docs/ together with all the figures in docs/figs/. This file can be run
through latex (or latex2html) to create the pdf (or html) manual.

Questions?

Ask Joakim Edsjo, edsjo@physto.se


