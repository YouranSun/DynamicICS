#ifndef INDEX_H
#define INDEX_H

#include <bits/stdc++.h>
#include "../gadget/Tools.h"
#include "../gadget/dsu.h"
#include "../core/Graph.h"

namespace INDEX{
    class Index{
    public:

        int level;
        int k_max;
        Graph *G;
        std::vector<int> x;
        std::vector<int> c_;
        std::vector<int> vis;
        std::vector<std::vector<int> > fa;
        std::vector<std::vector<int> > support;
        std::vector<std::vector<int> > pre;
        std::vector<std::vector<std::vector<int>  > > cv;
        DSU::Dsu dsu;

        Index();
        virtual ~Index();

        virtual void construct() = 0;
        virtual void readIndex(FILE *save) = 0;
        virtual void writeIndex(FILE *save) = 0;

        void init();
        void updateCore(int u, int k, std::vector<std::vector<int> > &S, std::vector<int> &U);
        void updateSupport(std::vector<int> &U);
        void dfs(int u, std::vector<std::vector<int> > &tree, std::vector<int> &S);
        void recompute(int k, std::vector<int> &V);
        void computeTree(int k);
        void formAnswer(std::vector<std::pair<int, int> > V, std::vector<std::vector<int> > &tree, std::vector<std::pair<double, std::vector<int> > > &ans);
        void computeIndex();
    } ;
}

#endif