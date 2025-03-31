#include "OrdIndex.h"

namespace ORDINDEX {
    OrdIndex::OrdIndex(){}
    
    OrdIndex::~OrdIndex(){}

    void OrdIndex::construct(){
        eprintf("Index construct starting stage\n");
        for (int u = 0; u < (G -> n); ++u) {
            for (auto v: (G -> nei[u]))
                x[u] += ((G -> core[v]) >= (G -> core[u]));
            c_[u] = (G -> core[u]);
        }
    
        level = -1;
        for (int i = 0; i < (G -> n); ++i) {
            int u = (G -> ord[i]);
            int k = c_[u]; c_[u] = -1;
            std::vector<int> U;
            std::vector<std::vector<int> > S(k + 1, std::vector<int>(1, u));
            ++level;
            updateCore(u, k, S, U);
            updateSupport(U);

            for (int j = 1; j <= k; ++j) cv[j].push_back(S[j]);
        }
    
        std::vector<int> tmp((G -> n), 0);
        for (int k = 1; k <= k_max; ++k) {
            eprintf("Computing index: k = %d\n", k);
            ++level;
            for (auto S: cv[k]) {
                for (auto u: S) {
                    tmp[u] = 0;
                    for (auto v: (G -> nei)[u]) {
                        if (vis[v] == level) ++support[v][k];
                        if ((G -> core)[v] >= k) ++tmp[u];
                    }
                    vis[u] = level;
                }
            }
            for (auto S: cv[k]) {
                for (auto u: S) {
                    for (auto v: (G -> nei)[u]) if ((G -> core)[v] >= k){
                        --tmp[v];
                        if (tmp[v] == k - 1) pre[v][k] = u;
                    }
                }
            }
            for (auto S: cv[k]) pre[S[0]][k] = S[0];
        }
    
        eputs("Construction End");
    }

    void OrdIndex::readIndex(FILE *save){
        eprintf("Itc ReadIndex Start\n");
        G = new Graph();
        G -> readIndex(save);
        dsu = DSU::Dsu(G -> n);
    
        support.resize(G -> n);
        pre.resize(G -> n);
        fa.resize(G -> n);
        for (int i = 0; i < (G -> n); ++i) {
            fscanf(save, "%d", &(G -> core[i]));
            support[i].resize((G -> core[i]) + 1);
            pre[i].resize((G -> core[i]) + 1);
            fa[i].resize((G -> core[i]) + 1, -1);
        }
    
        fscanf(save, "%d", &k_max);
        cv.resize(k_max + 1);
        
        for (int i = 1; i <= k_max; ++i) {
            int cnt; fscanf(save, "%d", &cnt);
            cv[i].resize(cnt);
            for (int j = 0; j < cv[i].size(); ++j) {
                int siz; fscanf(save, "%d", &siz);
                cv[i][j].resize(siz);
                for (auto &u: cv[i][j]) fscanf(save, "%d", &u);
                for (auto u: cv[i][j]) fscanf(save, "%d", &support[u][i]);
                for (auto u: cv[i][j]) fscanf(save, "%d", &pre[u][i]);
            }
        }
    
        c_.resize(G -> n);
        vis.resize(G -> n, -1);
        eprintf("Itc ReadIndex End\n");
    }

    void OrdIndex::writeIndex(FILE *save){
        G -> write(save);
        for (int i = 0; i < (G -> n); ++i) fprintf(save, "%d ", G -> core[i]); fprintf(save, "\n");
        fprintf(save, "%d\n", k_max);
        eprintf("k_max = %d\n", k_max);
        
        for (int i = 1; i <= k_max; ++i) {
            // fprintf(save, "%d %d\n", cv[i].size(), tree[cv[i][j][0]]);
            fprintf(save, "%d\n", (int)cv[i].size());
            for (int j = 0; j < cv[i].size(); ++j) {
                fprintf(save, "%d ", (int)cv[i][j].size());
                for (auto u: cv[i][j]) fprintf(save, "%d ", u);
                for (auto u: cv[i][j]) fprintf(save, "%d ", support[u][i]);
                for (auto u: cv[i][j]) fprintf(save, "%d ", pre[u][i]);
                fprintf(save, "\n");
            }
            fprintf(save, "\n");
        }
        fprintf(save, "%llu %llu\n", (G -> mn_quarter), (G -> mx_quarter));
    
        eputs("Write Index End");
        return;
    }
}