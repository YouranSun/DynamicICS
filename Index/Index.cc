#include "Index.h"

namespace INDEX {
    Index::Index(){
        level = -1;
        k_max = 0;
        G = nullptr;
        return;
    }
    
    Index::~Index(){
        delete G;
        return;
    }
    
    void Index::init(){
        eputs("Init Starts");
        level = -1;
        x.resize(G -> n, 0);
        c_.resize(G -> n, 0);
        vis.resize(G -> n, -1);
        support.resize(G -> n, std::vector<int>());
        pre.resize(G -> n, std::vector<int>());
        fa.resize(G -> n);
        dsu = DSU::Dsu(G -> n);
        
        k_max = 0;
        for (int i = 0; i < (G -> n); ++i) {
            k_max = std::max(k_max, (G -> core)[i]);
            support[i].resize((G -> core)[i] + 1, 0);
            pre[i].resize((G -> core)[i] + 1, 0);
            fa[i].resize((G -> core)[i] + 1, -1);
        }
        cv.resize(k_max + 1);
        return;
    }
    
    void Index::updateCore(int u, int k, std::vector<std::vector<int> > &S, std::vector<int> &U){
        if(c_[u] != -1) S[c_[u] + 1].push_back(u);
        U.push_back(u);
        vis[u] = level;
    
        for (auto v: (G -> nei[u])) if(c_[u] <= (G -> core[v])) {
            if (c_[v] == -1 || vis[v] == level) continue;
            if (c_[u] == -1 && c_[v] <= k || c_[u] != -1 && c_[v] == c_[u] + 1) {
                x[v] = x[v] - 1;
                if (x[v] < c_[v]) {
                    c_[v] = c_[v] - 1;
                    updateCore(v, k, S, U);
                }
            }
        }
    
        return;
    }
    
    void Index::updateSupport(std::vector<int> &U){
        for (auto u: U) {
            x[u] = 0;
            if(c_[u] == -1) continue;
            for (auto v: (G -> nei[u])) if(c_[u] <= (G -> core[v]) && c_[u] <= c_[v])
                x[u] = x[u] + 1;
        }
        return;
    }
    
    void Index::computeTree(int k){
        ++level;
        std::vector<int> V;
        for (auto S: cv[k]) for (auto u: S) V.push_back(u);
        dsu.init(V);
        for (int i = (int)cv[k].size() - 1; i >= 0; --i) {
            int key = cv[k][i][0];
            for (auto u: cv[k][i]) {
                dsu.merge(key, u);
            }
            if (fa[key].size() <= k) fa[key].push_back(-1);
        }
        for (int i = (int)V.size() - 1; i >= 0; --i) {
            int u = V[i];
            for (auto v: (G -> nei)[u]) {
                if ((G -> core)[v] >= k && vis[v] == level) {
                    int x = dsu.find(u), y = dsu.find(v);
                    // eprintf("k = %d u = %d v = %d x = %d y = %d\n", k, u, v, x, y);
                    if (x != y) {
                        dsu.merge(x, y);
                        fa[y][k] = x;
                    }
                }
            }
            vis[u] = level;
        }
    }
    
    void Index::dfs(int u, std::vector<std::vector<int> > &tree, std::vector<int> &S) {
        S.push_back(u);
        for (auto v: tree[u]) dfs(v, tree, S);
    }
    
    void Index::recompute(int k, std::vector<int> &V) {
        ++level;
        std::sort(V.begin(), V.end());
        dsu.init(V);
        for (auto u: V) vis[u] = level;
        if (k_max >= cv.size()) cv.push_back(std::vector<std::vector<int> >());
        else cv[k].clear();
        for (auto u: V) {
            c_[u] = 0;
            for (auto v: (G -> nei)[u]) if (vis[v] == level) ++c_[u];
        }
    
        std::set<Node> S;
        std::vector<std::pair<int, int> > edgesToAdd;
        for (auto u: V) S.insert(Node(u, G -> wei[u]));
        while (!S.empty()) {
            int u = (*S.begin()).u;
            std::vector<int> que;
            que.push_back(u); S.erase(S.begin());
            int id = cv[k].size();
            for (int i = 0; i < que.size(); ++i) {
                int v = que[i];
                for (auto w: (G -> nei)[v]) {
                    --c_[w];
                    if (c_[w] < k && S.find(Node(w, G -> wei[w])) != S.end()) {
                        edgesToAdd.push_back(std::make_pair(u, v));
                        que.push_back(w);
                        S.erase(Node(w, G -> wei[w]));
                    }
                }
            }
            cv[k].push_back(que);
        }
        
        computeTree(k);
    }
    
    void Index::computeIndex(){
        clock_t T_BEGIN = clock();
        init();
        construct();
        // constructTree();
        clock_t T_END = clock();
        eprintf("TIME USED IN ITC: %f(s)\n", 1.0 * (T_END - T_BEGIN) / CLOCKS_PER_SEC);
        return;
    }

    void Index::formAnswer(std::vector<std::pair<int, int> > V, std::vector<std::vector<int> > &tree, std::vector<std::pair<double, std::vector<int> > > &ans){
        ++level;
        std::vector<int> ver;
        for (auto p: V) ver.push_back(p.first);
        int lst = -1;
        // for (int i = (int)V.size() - 1; i >= 0; --i) {
        //     tree[V[i].first].clear();
        //     if (V[i].second) lst = V[i].first;
        //     else tree[lst].push_back(V[i].first);
        //     dsu.fa[V[i].first] = lst;
        // }
        for (int i = 0; i < V.size(); ++i) {
            int u = V[i].first;
            // for (auto v: (G -> nei)[u]) if (vis[v] == level) {
            //     int x = dsu.find(u), y = dsu.find(v);
            //     // eprintf("u = %d v = %d x = %d y = %d level = %d %d\n", u, v, x, y, level, vis[v]);
            //     if (x != y) {
            //         dsu.fa[y] = x;
            //         tree[x].push_back(y);
            //     }
            // }
            if (V[i].second) {
                std::vector<int> S;
                S.push_back(u);
                // dfs(u, tree, S);
                double sum = (G -> wei_ori)[u];
                for (int j = i - 1; j >= 0 && !V[j].second; --j) {
                    int v = V[j].first;
                    sum += (G -> wei_ori)[v];
                    S.push_back(v);
                }
                // for (auto u: S) sum += (G -> wei_ori)[u];
                // sum /= (int)S.size();
                sum /= (int)S.size();
                ans.push_back(std::make_pair(sum, S));
            }
            vis[u] = level;
        }
    }
}