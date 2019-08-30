#ifndef INCLUDE_H
#define INCLUDE_H

#include <stdio.h>
#include <string.h>
#include <malloc.h>
#include <stddef.h>
#include <math.h>

#include <R.h> 

#include "memory.h"
#include "matrix.h"
#include "soilwater.h"


#ifndef PI
#define PI 3.141593
#endif

#ifndef E
#define E 2.7182818
#endif

#ifndef MAX
#define MAX(X,Y) ((X)>(Y)?(X):(Y))
#endif

#ifndef MAX
#define MAX(X,Y) ((X)>(Y)?(X):(Y))
#endif

#ifndef TINY
#define TINY 1.0e-20
#endif

#ifndef SWAP
#define SWAP(a, b) {swap=(a); (a)=(b), (b)=swap;}
#endif
 
#endif
