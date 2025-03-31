#ifndef DYNAMIC_H
#define DYNAMIC_H

#include <bits/stdc++.h>

#include "../core/Graph.h"
#include "../core/core.h"
#include "../core/glist.h"
#include "../Index/Index.h"
#include "../gadget/treap2.h"

namespace DYNAMIC {
    class Dynamic{

    public:

        core::GLIST *cm;
        Graph *G;
        INDEX::Index *index;
        std::vector<int> ACL_TIME;
        std::vector<std::vector<int> > tree;
        long long candSize = 0, diffSize = 0;

        Dynamic();
        virtual ~Dynamic();
        virtual void init(const char* path) = 0;
        virtual void query(int u, int k, std::vector<std::pair<double, std::vector<int> > > &ans) = 0;
        virtual void remove(int u, int v) = 0;
        virtual void insert(int u, int v) = 0;
        virtual void saveIndex() = 0;
        void checkCore();
        void GET_TIME(double *TT);
    } ;

}

#endif