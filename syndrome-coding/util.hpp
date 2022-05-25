#ifndef UTIL_H
#define UTIL_H

/**
 * RNG, from: http://www.cse.yorku.ca/~oz/marsaglia-rng.html
 */
#define znew (z=36969*(z&65535)+(z>>16))
#define wnew (w=18000*(w&65535)+(w>>16))
#define MWC ((znew<<16)+wnew )
#define UC (unsigned char) /*a cast operation*/
typedef unsigned long UL;
static UL z=362436069, w=521288629;

void init_MWC()
{ 
  z=12345;w=65435;
}

#endif
