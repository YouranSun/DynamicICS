#include <bits/stdc++.h>
#include "../gadget/Tools.h"
#include "../core/Graph.h"

using namespace std;

int main(int argc, char **argv){
    FILE *file = fopen(argv[1], "r");
    Graph *G = new Graph();
    G -> readIndex(file);
    for (int i = 0; i < (G -> n); ++i) fscanf(file, "%d", &(G -> core[i]));

    std::vector<int> all;
    for (int u = 0; u < (G -> n); ++u) all.push_back(u);
    std::sort(all.begin(), all.end(), [&](const int &u, const int &v) { return (G -> core)[u] > (G -> core)[v]; });

    std::vector<int> deg;
    std::vector<int> par;
    std::vector<int> fa;
    deg.resize(G -> n);
    par.resize(G -> n, -1);

    int k_max; fscanf(file, "%d", &k_max);

    for (int k = 1; k <= k_max; ++k) {
        eprintf("Checking: k = %d\n", k);
        std::set<Node> S;
        std::set<int> T;
        for (int i = 0; i < (G -> n); ++i) {
            int u = all[i];
            if ((G -> core)[u] >= k) {
                deg[u] = 0;
                par[u] = u;
                for (auto v: (G -> nei)[u]) if ((G -> core)[v] >= k) ++deg[u];
                S.insert(Node(u, G -> wei[u]));
            }
            else break;
        }
        int cnt; fscanf(file, "%d", &cnt);
        for (int i = 0; i < cnt; ++i) {
            std::vector<int> cv;
            std::vector<int> support;
            std::vector<int> pre;
            int siz; fscanf(file, "%d", &siz);
            for (int j = 0; j < siz; ++j) {
                int u; fscanf(file, "%d", &u);
                cv.push_back(u);
            }
            for (int j = 0; j < siz; ++j) {
                int d; fscanf(file, "%d", &d);
                support.push_back(d);
            }
            for (int j = 0; j < siz; ++j) {
                int p; fscanf(file, "%d", &p);
                pre.push_back(p);
            }
            for (int j = 0; j < siz; ++j) {
                int u = cv[j], d = support[j], p = pre[j];
                if (j == 0) {
                    if (u != (S.begin()) -> u) {
                        eprintf("u = %d %d\n", u, (S.begin() -> u));
                    }
                    if (T.size()) {
                        for (auto v: T) eprintf("%d ", v); eputs("");
                    }
                    assert(T.size() == 0);
                    assert(u == p);
                    assert(u == (S.begin()) -> u);
                }
                else {
                    if (deg[u] >= k) eprintf("u = %d %d d = %d\n", u, deg[u], d);
                    assert(deg[u] < k);
                    T.erase(u);
                }
                S.erase(Node(u, G -> wei[u]));
                if (d != deg[u]) eprintf("u = %d d = %d deg = %d\n", u, d, deg[u]);
                // if (p != par[u]) eprintf("u = %d p = %d par = %d\n", u, p, par[u]);
                assert(d == deg[u]);
                // assert(p == par[u]);
                for (auto v: (G -> nei)[u]) if (S.find(Node(v, G -> wei[v])) != S.end()){
                    --deg[v];
                    if (deg[v] == k - 1) {
                        par[v] = u;
                        T.insert(v);
                    }
                }
            }
        }
    }

    fclose(file);
    return 0;
}