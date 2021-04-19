# README
Prepare shared library (TreeManager_C.so) as follows:
	root -l
	.L TreeManager.C++

Change the location of the file at the top of drawHistosFromEvent.C to:
  R__LOAD_LIBRARY(fromUmut/GluedCubes/Analysis/TreeManager_C.so);

Example of how to run: 
  root -l 
  .L drawHistosFromEvent.cpp
  showHighLYEvent(3500)

The argument for the function is the cut you want to make on the light yield.
To show the total light yield distribution to determine a sensible value for this cut just run the "totalLightYield()" function.
