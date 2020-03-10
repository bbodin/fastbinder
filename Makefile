C_FILES := $(shell find ./fastbinder/src -name \*.\*pp -not -name RcppExports.cpp)
R_FILES := $(shell find ./fastbinder/R ./fastbinder/tests -name \*.R -not -name RcppExports.R)

R_CMD := R -q

PACKAGE_VERSION=1.0.0

all : fastbinder.Rinstall/fastbinder/libs/fastbinder.so  tests_cpp/test_binders  fastbinder.pdf 


infos :
	@echo "C_FILES=${C_FILES}"
	@echo "R_FILES=${R_FILES}"


check : fastbinder/src/RcppExports.cpp
	${R_CMD} -e  "devtools::check(\"fastbinder\");"

%/NAMESPACE : %/R/fastbinderAPI.R %/DESCRIPTION
	rm -f $*/man/*
	echo "# Generated by roxygen2: do not edit by hand" > $*/NAMESPACE
	${R_CMD} -e  "library(devtools) ; document(\"$*\");"
	find $*/man/ -type f -exec sed -i 's/[\][\][\]%/\\%/gI' {} \;

%/src/RcppExports.cpp  %/R/RcppExports.R : % %/NAMESPACE ${C_FILES} ${R_FILES}
	rm -f $*/src/RcppExports.cpp  $*/R/RcppExports.R
	${R_CMD} -e  "Rcpp::compileAttributes(pkgdir = \"$*\" , verbose=TRUE);"


%_${PACKAGE_VERSION}.tar.gz : ${C_FILES} ${R_FILES} %/src/RcppExports.cpp  %/R/RcppExports.R  
	rm -rf fastbinder/src/*.o ./fastbinder/src/*.so 
	R CMD build ./$*

%.Rinstall/fastbinder/libs/fastbinder.so : %_${PACKAGE_VERSION}.tar.gz 
	mkdir -p $*.Rinstall
	rm $*.Rinstall/* -rf
	R CMD INSTALL  -l $*.Rinstall $<

%_${PACKAGE_VERSION}.pdf : %/NAMESPACE
	${R_CMD} -e  "library(devtools) ; devtools::build_manual(\"$*\"); " || ${R_CMD} -e  "library(devtools) ; devtools::check(\"$*\",manual=TRUE); " || touch $@

%.pdf : %/NAMESPACE
	R CMD Rd2pdf $* --no-preview --force

tests_cpp/% : tests_cpp/%.cpp
	make -C tests_cpp/ $*

deps :
	echo "To be defined."

clean : 
	rm -rf current *~ *.Rinstall *.pdf  *.tar.gz *.Rcheck ./fastbinder/NAMESPACE ./fastbinder/src/*.o ./fastbinder/src/*.so 	./fastbinder/src/*.rds ./fastbinder/src/RcppExports.cpp  ./fastbinder/R/RcppExports.R  ./fastbinder/man/*.Rd tests_cpp/testfastbinder

.PHONY: clean tests_cpp/testfastbinder
.SECONDARY: