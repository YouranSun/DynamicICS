#include "ICPDynamic.h"
#include "../Index/ICPIndex.h"

namespace ICPDYNAMIC {
    ICPDynamic::ICPDynamic(){
        index = new ICPINDEX::ICPIndex();
    }

    ICPDynamic::~ICPDynamic(){}

    void ICPDynamic::init(const char* path) {
        eprintf("KeyMaintenane init path = %s\n", path);
        FILE *file = fopen(path, "r");
        index -> readIndex(file);
        G = (index -> G);
        for (int i = 0; i < (G -> n); ++i) (G -> core)[i] = 0;
        cm = new core::GLIST(G -> n);
        cm -> ComputeCore(G -> nei, true, G -> core);
        fclose(file);
    
        ACL_TIME.resize((index -> k_max) + 1, 0);
        tree.resize(G -> n, std::vector<int>());
    
        R.resize(G -> n);
        for (int u = 0; u < (G -> n); ++u) R[u].resize((G -> core)[u] + 1);
        for (int i = 1; i <= (index -> k_max); ++i) {
            for (int j = 0; j < (index -> cv)[i].size(); ++j) {
                for (auto u: (index -> cv)[i][j]) {
                    R[u][i] = j;
                }
            }
        }
        // eputs("ok");
        L.resize(G -> n, 0);
    
        return;
    }

    void ICPDynamic::insertionDFS(int u, int k, int R_min) {
        L[u] = -1;
        for (auto v: (G -> nei)[u]) {
            if (R[v].size() <= k || R[v][k] != R_min || L[v] == -1) continue;
            --L[v];
            if (L[v] < k) insertionDFS(v, k, R_min);
        }
    }
    
    bool ICPDynamic::isRecompute(int u, int v, int k) {
        int R_min = std::min(R[u][k], R[v][k]);
        int w = R[u][k] < R[v][k] ? u : v;
        for (auto x: (index -> cv)[k][R_min]) {
            L[x] = 0;
            for (auto y: (G -> nei)[x]) if (R[y].size() > k){
                L[x] += (R[y][k] >= R[x][k]);
            }
        }
        insertionDFS((index -> cv)[k][R_min][0], k, R_min);
        return (L[w] != -1);
    }
    
    void ICPDynamic::recompute(int k) {
        std::vector<int> V;
        cm -> GetAtLeast(k, V);
        candSize += V.size();
        index -> recompute(k, V);
        for (int i = 0; i < (index -> cv)[k].size(); ++i) {
            for (auto u: (index -> cv)[k][i]) {
                if (R[u].size() <= k) R[u].push_back(i);
                else R[u][k] = i;
            }
        }
    }

    void ICPDynamic::insert(int u, int v) {        
        int cu = (G -> core)[u], cv = (G -> core)[v];
        int c_min = std::min(cu, cv);
        cm -> InsertUpd(u, v, (index -> G) -> nei, (index -> G) -> core, (index -> k_max));
        for (int i = 1; i <= c_min; ++i) {
            if (RECOMPUTE || isRecompute(u, v, i)) recompute(i);
        }
        if ((G -> core)[u] > cu || (G -> core[v]) > cv) recompute(c_min + 1);
        return;      
    }

    void ICPDynamic::remove(int u, int v) {
        int cu = (G -> core)[u], cv = (G -> core)[v];
        int cmin = std::min(cu, cv);
        cm -> RemoveUpd(u, v, G -> nei, G -> core, index -> k_max);
        for (int i = 1; i <= cmin; ++i) {
            if (RECOMPUTE) {
                recompute(i);
                continue;
            }
            L[u] = 0;
            for (auto y: (G -> nei)[u]) if (R[y].size() > i) {
                L[u] += (R[y][i] >= R[u][i]);
            }
            L[v] = 0;
            for (auto y: (G -> nei)[v]) if (R[y].size() > i) {
                L[v] += (R[y][i] >= R[v][i]);
            }
            if (L[u] < i || L[v] < i)
                recompute(i);
        }
    }

    void ICPDynamic::queryDfs(int u, int k, std::vector<int> &S) {
        // eprintf("queryDfs u = %d k = %d\n", u, k);
        for (auto v: (index -> cv)[k][R[u][k]]) S.push_back(v);
        for (auto v: tree[u]) /*eprintf("u = %d v = %d\n", u, v),*/ queryDfs(v, k, S);
        return;
    }

    void ICPDynamic::query(int k, int r, std::vector<std::pair<double, std::vector<int> > > &ans){
        ans.clear();
        r = std::min(r, (int)(index -> cv)[k].size());
        for (int i = (int)(index -> cv)[k].size() - 1; i >= (int)(index -> cv)[k].size() - r; --i){
            int u = (index -> cv)[k][i][0];
            tree[u].clear();
        }
        for (int i = (int)(index -> cv)[k].size() - 1; i >= (int)(index -> cv)[k].size() - r; --i){
            // eprintf("i = %d\n", i);
            int u = (index -> cv)[k][i][0];
            if (~(index -> fa)[u][k]) tree[(index -> fa)[u][k]].push_back(u);
            std::vector<int> S;
            queryDfs(u, k, S);
            double sum = 0;
            for (auto u: S) sum += (G -> wei_ori)[u];
            sum /= (int)S.size();
            ans.push_back(std::make_pair(sum, S));
        }
        return;
    }

    void ICPDynamic::saveIndex(){}
}