#ifndef TOOLS_H_
#define TOOLS_H_

#include <time.h>
#include <stdio.h>
#include <vector>

#define eprintf(...) fprintf(stderr, __VA_ARGS__)
#define eputs(str) fputs(str"\n", stderr)

void KEEPTIME(const char *str, clock_t &stamp);

struct Node {
    int u, wei;
    Node (const int &u_ = 0, const int &wei_ = 0): u(u_), wei(wei_){}
    bool operator < (const Node &oth) const{ return wei < oth.wei; }
    bool operator > (const Node &oth) const{ return wei > oth.wei; }
} ;

#endif