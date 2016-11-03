#ifndef UTIL_FUNCS_H 
#define UTIL_FUNCS_H

#include <vector>
#include <string>
#include <utility>



#include "Suffix_t.h"

class ReadsManipulator;       // forward decl.

enum {TUMOUR, HEALTHY, SWITCHED};
enum {LEFT, RIGHT};         // distinguish, paired end type

static const int MIN_SUFFIX = 30;

int computeLCP(Suffix_t &isuf, Suffix_t &jsuf, ReadsManipulator &reads);
// Returns the longest common prefix between isuf and jsuf suffixes

#endif
