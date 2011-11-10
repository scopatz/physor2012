#! /bin/bash
NAME="scopatz_physor2012"
latex ${NAME}.tex 
bibtex ${NAME}.aux 
bibtex ${NAME}.aux 
latex ${NAME}.tex 
latex ${NAME}.tex 
dvipdf ${NAME}.dvi

