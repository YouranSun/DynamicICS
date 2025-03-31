#ifndef CORE_GADGET_DSU_H
#define CORE_GADGET_DSU_H

#include <vector>

namespace DSU{

class Dsu{
public:
    std::vector<int> fa;
    Dsu(int n = 0) { 
        fa.resize(n);
        for (int i = 0; i < n; ++i) fa[i] = i;
    }

    void init(std::vector<int> vec) {
        for (auto u: vec) fa[u] = u;
        return;
    }

    int find(int x) { return (x == fa[x] ? x : (fa[x] = find(fa[x]))); }

    void merge(int x, int y) {
        x = find(x), y = find(y);
        fa[y] = x;
    }
} ;

}

#endif