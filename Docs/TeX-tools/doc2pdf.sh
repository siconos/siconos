#!/bin/sh 

# -----------------------------------------------------------------------------
# doc2pdf v1.3
# SICONOS WP2 - Projet BIPOP - INRIA
# JB CHARLETY - 16/11/2004
#
# Script which creates directly pdf document from LaTeX
#
# -----------------------------------------------------------------------------

Version=1.3

# -----------------------------------------------------------------------------
# PATHS
# -----------------------------------------------------------------------------

# Path to script glosstex
GlosstexPath=/Users/dubois/Tex
TeXtoolsPath='/Users/dubois/WORK/CODING/SICONOS/TeX-tools'

if !(test -e $GlosstexPath/glosstex.csh) then
    echo "doc2pdf --- Please edit doc2pdf.sh and set variable GlosstexPath"
    exit
fi

if !(test -e $TeXtoolsPath/acronym.gdf) then
    echo "doc2pdf --- Please edit doc2pdf.sh and set variable TeXtoolsPath"
    exit
fi


# -----------------------------------------------------------------------------
# COMMANDS
# -----------------------------------------------------------------------------

#echo parameters are $@

while getopts c:d:z:v'-version''-help'l: option
  do
#  echo option= $option
  case $option in
      c)   	
	  $GlosstexPath/glosstex.csh $2 $TeXtoolsPath/acronym.gdf $TeXtoolsPath/notation.gdf;	
	  dvips $2.dvi -o $2.ps;	
	  ps2pdf -sPAPERSIZE=a4 $2.ps;
	  echo ;
	  echo "doc2pdf --- pdf created"
	  echo ;
	  ;;
      
      d)
	  if test -e ./$2.tex
	      then
	      echo "doc2pdf --- TeX file exists"
	      for i in $2.*
		do
		if test "$i" != "$2.tex"
		    then
		    rm $i
		fi
	      done
	      echo "doc2pdf --- repertory cleaned"
	  else
	      echo "doc2pdf --- file does not exist"
	  fi
	  ;;
      
      z)
	  if test -d ../$2
	      then
	      echo "doc2pdf --- repertory exists"
	      cd ..
	      tar -cvzf $2.tar.gz $2
	      cd $2
	      echo "doc2pdf --- archive created"
	  else
	      echo "doc2pdf --- repertory does not exist"
	  fi
	  ;;	
      
      l)
	   if test -e ./$2.tex
	      then
	      echo "doc2pdf --- TeX file exists"
	      $GlosstexPath/glosstex.csh $2;
	      echo "doc2pdf --- dvi created"
	  else
	      echo "doc2pdf --- file does not exist"
	  fi
	  ;;

	  v)
	  echo ;
	  echo doc2pdf $Version --- JB CHARLETY - November 2004 - INRIA
	  echo ;
	  ;; 	
      
      h)
	  echo ;
	  echo "doc2pdf $Version --- Options"
	  echo "-c DocName	: compiles DocName.tex and creates DocName.pdf"
	  echo "-d DocName	: cleans DocName.pdf, and TeX intermediate files"
	  echo "-d DocName	: compiles DocName.tex and creates DocName.dvi"
	  echo "-z DirName	: creates an archive DirName.tar.gr of the directory given in parameter"
	  echo "-version	: gives version number"
	  echo "-help 		: gives parameters options"
	  echo ;
	  ;; 
  esac
done

# -----------------------------------------------------------------------------
# CHANGES LOG
#
# v1.3 - 06/12/2004 : new option "l" which compiles the TeX file. 
#
# v1.2 - 19/11/2004 : new option "z" which creates an archive tgz of the repertory
#
# v1.1 - 17/11/2004 : verify that glosstex script is present in repertory indicated
# in variable GlosstexPath.
#
# -----------------------------------------------------------------------------




