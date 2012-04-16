#!/bin/bash
PRESNAME="scopatz_physor2012_pres"

rst2pdf ${PRESNAME}.rst -b1 -s slides.style -o ${PRESNAME}.pdf --fit-background-mode=center
