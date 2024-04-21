% NAME:      oldeps.m (script file)
% FUNCTION:  determine/estimate EPS like in the old times
% SYNTAX:    oldeps
% REFERENCE: C. Moler: Numerical Computing with MATLAB, SIAM, 2004
% VERSION:   1. original version
% AUTHOR:    J. Behrens (behrens@ma.tum.de)

%--- switch to long format
format long

%--- this is the crucial calculation with roundoff
a= 4/ 3;

%--- the rest can be done exactly in floating point representation
b= a- 1;
c= 3* b;
e= 1- c
%--- end of oldeps.m