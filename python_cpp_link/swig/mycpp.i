%module mycpp
%{
#include "mycpp.cpp"
%}

%include "std_string.i"
%include "mycpp.cpp"
