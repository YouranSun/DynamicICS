#include "OrdDynamic.h"
#include <assert.h>

namespace ORDDYNAMIC {
    int tP;
    int tP2;
    int tQ;
    int tF;
    int tR;
    int tS;
    std::vector<int> inP;
    std::vector<int> inP2;
    std::vector<int> inQ;
    std::vector<int> inW;
    std::vector<int> inR;
    std::vector<int> inS;
    std::vector<int> cover;
    std::vector<int> rank;

    OrdDynamic::OrdDynamic(){
        index = new ORDINDEX::OrdIndex();
    }

    OrdDynamic::~OrdDynamic(){
        for (int i = 0; i <= (index -> k_max); ++i) delete(trpOrd[i].first);
    }

    int OrdDynamic::findNext(int k, Node x) { // get first key node weight >= wei[u]
        return (trpOrd[k].first -> FindNextKey(x, trpOrd[k].second));
    }
    
    void OrdDynamic::init(const char* path) {
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
        trpOrd.resize((index -> k_max) + 1, std::make_pair(nullptr, 0));
        trpId.resize(G -> n);
        for (int u = 0; u < (G -> n); ++u) {
            trpId[u].resize((G -> core[u]) + 1);
        }

        mxKey.resize(G -> n, 0);
        deg.resize(G -> n, 0);

        tP = 0;
        tP2 = 0;
        tQ = 0;
        tR = 0;
        tS = 0;
        inR.resize((G -> n), 0);
        inP.resize((G -> n), 0);
        inP2.resize((G -> n), 0);
        inQ.resize((G -> n), 0);
        inW.resize((G -> n), 0);
        inS.resize((G -> n), 0);
        rank.resize((G -> n), 0);
        cover.resize((G -> n), 0);
        insPrev.resize((G -> n), -1);

        for (int i = 1; i <= (index -> k_max); ++i) {
            ++tP;
            trpOrd[i].first = new TREAP2::Treap();
            for (auto S: (index -> cv)[i]) {
                mxKey[S[0]] = std::max(mxKey[S[0]], i);
                for (int j = 0; j < S.size(); ++j) {
                    int u = S[j];
                    trpId[u][i] = trpOrd[i].first -> NewTreeNode(Node(u, G -> wei[u]), j == 0);
                    trpOrd[i].first -> Insert(trpId[u][i], false, trpOrd[i].second);
                    if (j == 0) mxKey[u] = i;
                    for (auto v: (G -> nei)[u]) if ((G -> core)[v] >= i){
                        if (inP[v] < tP) (index -> fa)[v][i] = u;
                    }
                    inP[u] = tP;
                }
            }
        }
        return;
    }

    void OrdDynamic::printIndex(int k){
        eprintf("Print Index: k = %d root = %d\n", k, trpOrd[k].second);
        std::vector<int> tmp;
        trpOrd[k].first -> TraverseToVec(trpOrd[k].second, tmp);
        for (auto u: tmp) eprintf("(%d %d %d) ", u, (index -> support)[u][k], (index -> pre)[u][k]);
        eputs("");
        trpOrd[k].first -> Print();
        eputs("");
    }

    void OrdDynamic::insertShrink(int k, int s, std::set<Node> &P, int &prev, std::vector<int> &waitList) {
        std::queue<int> que;
        que.push(s);
        while (!que.empty()) {
            int u = que.front();
            que.pop();
            if (inP[u] == false) continue;
            P.erase(Node(u, G -> wei[u]));
            inP[u] = false;
            inW[u] = true;
            insPrev[u] = prev;
            // eprintf("u = %d prev = %d deg = %d\n", u, prev, deg[u]);
            prev = trpId[u][k];
            waitList.push_back(u);
            for (auto v: (G -> nei)[u]) if (inP[v] == tP) {
                // eprintf("  v = %d deg = %d\n", v, deg[v]);
                --deg[v];
                if (deg[v] < k) {
                    if (deg[v] == k - 1) (index -> pre)[v][k] = u;
                    que.push(v);
                }
            }
        }
    }

    void OrdDynamic::insertUpdate(int k, std::vector<int> waitList, int x, int y) {
        if (((G -> core)[x] < k || (G -> core)[y] < k) && !waitList.size()) return;
    
        std::set<Node> P;
        std::priority_queue<std::pair<int, int>, std::vector<std::pair<int, int> >, std::greater<std::pair<int, int> > > Q;
    
        ++tP;
        ++tQ;
    
        for (auto w: waitList) {
            P.insert(Node(w, G -> wei[w]));
            inP[w] = tP;
            deg[w] = 0;
            trpId[w].push_back(trpOrd[k].first -> NewTreeNode(Node(w, G -> wei[w]), false));
            (index -> pre)[w].push_back(0);
        }
    
        for (auto w: waitList) {
            for (auto z: (G -> nei)[w]) {
                if ((G -> core[z]) >= k) {
                    if(inP[z] != tP && inQ[z] != tQ) {
                        Q.push(std::make_pair(trpOrd[k].first -> Rank(trpId[z][k]), z));
                        inQ[z] = tQ;
                    }
                    ++deg[w];
                }
            }
        }
    
        if (inP[x] != tP && inP[y] != tP) {
            int rx = trpOrd[k].first -> Rank(trpId[x][k]), ry = trpOrd[k].first -> Rank(trpId[y][k]);
            if (rx > ry) std::swap(x, y), std::swap(rx, ry);
            ++(index -> support)[x][k];
            Q.push(std::make_pair(rx, x));
            inQ[x] = tQ;
        }
    
        waitList.clear();
    
        int ptr = 0;
        while (!Q.empty()) {
            ++candSize;
            int u = (P.size() ? (*P.begin()).u : -1), v = Q.top().second;

            int nextKey = (~u ? findNext(k, Node(u, G -> wei[u])) : 0);
    
            if (nextKey && (trpOrd[k].first -> Rank(nextKey)) <= (trpOrd[k].first -> Rank(trpId[v][k]))) {
                P.erase(u);
                nextKey = ~nextKey;
                insertShrink(k, u, P, nextKey, waitList);
                mxKey[u] = k;
                (index -> pre)[u][k] = u;
            }
            else {
                int rkv = Q.top().first;
                Q.pop();
                inQ[v] = 0;
                int ndeg = (index -> support)[v][k];
                for (auto w: (G -> nei)[v]) {
                    if (inP[w] == tP) ++ndeg, --deg[w];
                }
                deg[v] = ndeg;
    
                if (mxKey[v] >= k || ndeg < k) {
                    (index -> support)[v][k] = ndeg;
                    int prev = trpId[v][k];
                    for (auto w: (G -> nei)[v]) {
                        if (inP[w] == tP && deg[w] < k) insertShrink(k, w, P, prev, waitList);
                    }
                    continue;
                }
                for (auto w: (G -> nei)[v]) {
                    if (inP[w] == tP) ++deg[w];
                    else if ((G -> core)[w] >= k && inW[w] == false) {
                        int rkw = (trpOrd[k].first -> Rank(trpId[w][k]));
                        if(rkw > rkv && inQ[w] != tQ) {
                            Q.push(std::make_pair(rkw, w));
                            inQ[w] = tQ;
                        }
                    }
                }
                inP[v] = tP;
                P.insert(Node(v, G -> wei[v]));
    
    
            }
        }
    
        while (!P.empty()) {
            int u = (*P.begin()).u;
            int prev = findNext(k, Node(u, G -> wei[u]));
            prev = ~prev;
            insertShrink(k, u, P, prev, waitList);
            mxKey[u] = k;
        }
    
        for (auto u: waitList) {
            if ((index -> support)[u].size() <= k) {
                (index -> support)[u].push_back(deg[u]);
            }
            else {
                (index -> support)[u][k] = deg[u];
                trpOrd[k].first -> Delete(trpId[u][k], trpOrd[k].second);
            }
            trpId[u][k] = (trpOrd[k].first -> TreeNode(trpId[u][k], Node(u, G -> wei[u]), mxKey[u] >= k));
            if (insPrev[u] >= 0) {
                trpOrd[k].first -> InsertAfter(trpId[u][k], insPrev[u], trpOrd[k].second);
            }
            else {
                trpOrd[k].first -> InsertBefore(trpId[u][k], ~insPrev[u], trpOrd[k].second);
            }
            inW[u] = false;
        }
    }
    
    void OrdDynamic::insert(int u, int v) {
        std::vector<std::vector<int> > core_inc((index -> k_max) + 2, std::vector<int>());
        cm -> InsertRem(u, v, (index -> G) -> nei, (index -> G) -> core, core_inc);
        if(core_inc[(index -> k_max) + 1].size()) {
            ++(index -> k_max);
            ACL_TIME.push_back(0);
            trpOrd.push_back(std::make_pair(new TREAP2::Treap(), 0));
        }
        for (int k = 1; k <= (index -> k_max); ++k) {
            clock_t START = clock();
            insertUpdate(k, core_inc[k], u, v);
            clock_t END = clock();
            ACL_TIME[k] += END - START;
        }
    }

    int OrdDynamic::getRank(int k, int u) { 
        if (inR[u] == tR) return rank[u];
        inR[u] = tR;
        return rank[u] = (trpOrd[k].first -> Rank(trpId[u][k])); 
    }
    
    bool OrdDynamic::removeShrink(int k, int s, std::unordered_set<int> &P, int &prev, std::vector<int> &waitList) {
        clock_t BEG = clock();
    
        bool flg = true;
        
        std::queue<int> que;
        ++tP2;
    
        int ndeg = (index -> support)[s][k];
    
        for (auto u: P) {
            deg[u] = (index -> support)[u][k];
            for (auto v: (G -> nei)[u]) {
                if (inS[v] == tS && getRank(k, v) < getRank(k, u) && getRank(k, v) >= getRank(k, s)) ++deg[u];
                if (v == s) ++ndeg;
            }
            if (deg[u] < k) que.push(u), inP2[u] = tP2;
        }
    
        while (!que.empty()) {
            int u = que.front();
            que.pop();
            for (auto v: (G -> nei)[u]) if ((G -> core)[v] >= k) {
                if ((G -> core)[u] >= k && inP[v] == tP){
                    --deg[v];
                    if (deg[v] < k && inP2[v] != tP2) {
                        inP2[v] = tP2;
                        que.push(v);
                    }
                }
                if (v == s) {
                    --ndeg;
                    if (ndeg < k) break;
                }
            }
            if (ndeg < k) break;
        }
        
        if (!flg || ndeg < k) {
            return false;
        }
    
        (index -> support)[s][k] = ndeg;
        for (auto u: P) if (inP2[u] != tP2) {
            std::vector<int> vnei;
            for (auto v: (G -> nei)[u]) if ((v == s || inS[v] == tS) && getRank(k, v) >= getRank(k, s)) vnei.push_back(v);
            std::sort(vnei.begin(), vnei.end(), [&](const int &i, const int &j) { return getRank(k, i) < getRank(k, j); });
            for (auto v: vnei) {
                --deg[u];
                if (deg[u] == k - 1) {
                    (index -> pre)[u][k] = v;
                    que.push(u);
                }
                if (getRank(k, v) > getRank(k, prev)) prev = v;
            }
        }
        prev = trpId[prev][k];
    
        while (!que.empty()) {
            int u = que.front();
            que.pop();
            if (inP[u] != tP) continue;
            inP[u] = false;
            insPrev[u] = prev;
            prev = trpId[u][k];
            cover[u] = false;
            P.erase(u);
            waitList.push_back(u);
            inW[u] = true;
            for (auto v: (G -> nei)[u]) {
                if (inP[v] == tP && inP2[v] != tP2){
                    --deg[v];
                    if (deg[v] == k - 1) {
                        (index -> pre)[v][k] = u;
                        que.push(v);
                    }
                }
            } 
        }
    
        return true;
    }
    
    int OrdDynamic::findPrev(int k, int u) {
        return (trpOrd[k].first -> Get((trpOrd[k].first -> FindPrevKey(getRank(k, u), trpOrd[k].second)))).u;
    }
    
    bool OrdDynamic::checkToQue(int w, int k) {
        if (mxKey[w] >= k) return true;
        int ndeg = (index -> support)[w][k];
        int prevKey = findPrev(k, w);
        // eprintf("prevKey = %d\n", prevKey);
        assert(mxKey[prevKey] >= k);
        int rkw = getRank(k, w);
        int rkKey = getRank(k, prevKey);
        assert(rkKey <= rkw);
        for (auto u: (G -> nei)[w]) if (inW[u] == false && (index -> support)[u].size() > k) {
            if (getRank(k, u) > rkw && inP[u] == tP) {
                --ndeg;
            }
            if (getRank(k, u) < rkw && inQ[u] != tQ && inP[u] != tP && getRank(k, u) >= getRank(k, prevKey)) {
                ++ndeg;
            }
        }
        if (ndeg >= k) return false;
        return true;
    }
    
    void OrdDynamic::removeUpdate(int k, std::vector<int> waitList, int x, int y) {
        if (((G -> core)[x] < k || (G -> core)[y] < k) && !waitList.size()) return;
    
        std::set<std::pair<int, int> > S;
        std::unordered_set<int> P;
        std::priority_queue<std::pair<int, int>, std::vector<std::pair<int, int> > > Q;
    
        ++tP;
        ++tQ;
        ++tR;
        ++tS;
    
        for (auto w: waitList) {
            inQ[w] = tQ;
            Q.push(std::make_pair(getRank(k, w), w));
        }
    
        if ((G -> core)[x] + (inQ[x] == tQ) >= k && (G -> core)[y] + (inQ[y] == tQ) >= k) {
            int rx = getRank(k, x), ry = getRank(k, y);
            if (rx > ry) std::swap(x, y), std::swap(rx, ry);
            --(index -> support)[x][k];
            deg[x] = (index -> support)[x][k];
            if (inQ[x] != tQ) {
                Q.push(std::make_pair(rx, x));
                inQ[x] = tQ;
            }
            if (rx >= getRank(k, (index -> pre)[y][k])) {
                Q.push(std::make_pair(ry, y));
                inQ[y] = tQ;
            }
        }
    
        waitList.clear();
    
        while (!Q.empty()) {
            ++candSize;
            int v = Q.top().second;
            int rkv = Q.top().first;
            inQ[v] = false;
            cover[v] = true;
            if (inS[v] == tS) inS[v] = false, S.erase(std::make_pair(rkv, v));
            Q.pop();
                
            int prevKey = findPrev(k, v);
            while (S.size() && (G -> wei)[prevKey] < (G -> wei)[findPrev(k, (*S.rbegin()).second)]) {
                int w = (*S.rbegin()).second;
                S.erase(std::make_pair(getRank(k, w), w));
                inS[w] = false;
                for (auto u: (G -> nei)[w]) if (inP[u] == tP && getRank(k, w) < getRank(k, u)) {
                    ++(index -> support)[u][k], --(index -> support)[w][k];
                }
            }
    
            int ndeg = (index -> support)[v][k], odeg = ndeg;
            for (auto w: (G -> nei)[v]) {
                if (inP[w] == tP && getRank(k, w) > rkv) {
                    --ndeg;
                    if ((G -> core)[v] >= k && (G -> core)[w] >= k) {
                        ++(index -> support)[w][k];
                        ++deg[w];
                    }
                }
            }
    
            if ((G -> core)[v] >= k) (index -> support)[v][k] = ndeg;
            else (index -> support)[v][k] = 0;
    
            if ((G -> core[v]) >= k && mxKey[v] >= k) {
                deg[v] = ndeg; 
                int prev = v;
                if (removeShrink(k, v, P, prev, waitList)) {
                    cover[v] = false;
                    continue;
                }
            }
    
            bool flag = false;
            std::queue<int> que;
            for (auto w: (G -> nei)[v]) {
                if (inP[w] == tP) {
                    if ((G -> core)[v] >= k && (G -> core)[w] >= k && getRank(k, w) > rkv) {
                        --(index -> support)[w][k];
                        --deg[w];
                    }
                }
                else if ((G -> core)[w] >= k && inW[w] == false){
                    int rkw = getRank(k, w);
                    if (rkw < rkv && inQ[w] != tQ) {
                        if (checkToQue(w, k) || true) {
                            Q.push(std::make_pair(rkw, w));
                            inQ[w] = tQ;
                            if (inS[w] == tS) {
                                S.erase(std::make_pair(rkw, w));
                                inS[w] = false;
                            }
                        } else {
                            inS[w] = tS;
                            S.insert(std::make_pair(rkw, w));
                            int prevKey = findPrev(k, w);
                            if (inQ[prevKey] != tQ && inP[prevKey] != tP) {
                                Q.push(std::make_pair(getRank(k, prevKey), prevKey));
                                inQ[prevKey] = tQ;
                            }
                        }
                    }
                    else if (rkv >= getRank(k, (index -> pre)[w][k])/*(index -> pre)[w][k] == v*/ && inP[w] != tP) {
                        if (inQ[w] != tQ && inW[w] == false) {
                            que.push(w);
                            flag = true;
                        }
                    }
                }
            }
    
            if (v == x && (G -> core[y]) >= k && inW[y] == false && (index -> pre)[y][k] == x && inP[y] != tP) {
                assert(inQ[y] != tQ);
                que.push(y);
                flag = true;
            }
    
            if (!flag) {
                inP[v] = tP;
                P.insert(v);
                if (mxKey[v] >= k) {
                    assert(mxKey[v] == k);
                    --mxKey[v];
                }
                deg[v] = (index -> support)[v][k];
                for (auto w: (G -> nei)[v]) {
                    if (inP[w] == tP && (G -> core)[v] >= k && (G -> core)[w] >= k) {
                        if (getRank(k, w) > rkv) {
                            ++(index -> support)[w][k];
                            ++deg[w];
                        }
                        ++(index -> support)[v][k], ++deg[v];
                    } else if ((G -> core)[w] >= k){
                        int rkw = getRank(k, w);
                        if (rkw < rkv && inQ[w] != tQ) {
                            if (checkToQue(w, k)) {
                                Q.push(std::make_pair(rkw, w));
                                inQ[w] = tQ;
                            } else {
                                inS[w] = tS;
                                S.insert(std::make_pair(rkw, w));
                                int prevKey = findPrev(k, w);
                                if (inQ[prevKey] != tQ && inP[prevKey] != tP) {
                                    Q.push(std::make_pair(getRank(k, prevKey), prevKey));
                                    inQ[prevKey] = tQ;
                                }
                            }
                        }
                    }
                }
            }
            else {
                while (!que.empty()) {
                    int u = que.front(); que.pop();
                    if (inQ[u] == tQ) continue;
                    int rku = getRank(k, u);
                    Q.push(std::make_pair(rku, u));
                    inQ[u] = tQ;
                    if (inS[u] == tS) {
                        S.erase(std::make_pair(rku, u));
                        inS[u] = false;
                    }
                    for (auto w: (G -> nei)[u]) if ((G -> core[w]) >= k && inW[w] == false) {
                        if ((index -> pre)[w][k] == u && inP[w] != tP) {
                            assert(inQ[w] != tQ);
                            que.push(w);
                        }
                    }
                    if (u == x && (G -> core[y]) >= k && inW[y] == false && (index -> pre)[y][k] == x && inP[y] != tP) {
                        assert(inQ[y] != tQ);
                        que.push(y);
                    }
                }
                (index -> support)[v][k] = odeg;
                inQ[v] = tQ;
                Q.push(std::make_pair(rkv, v));
                if (inS[v] == tS) {
                    S.erase(std::make_pair(rkv, v));
                    inS[v] = false;
                }
            }
        }
    
        while (S.size()) {
            int w = (*S.rbegin()).second;
            S.erase(std::make_pair(getRank(k, w), w));
            inS[w] = false;
            for (auto u: (G -> nei)[w]) if (inP[u] == tP && getRank(k, w) < getRank(k, u))
                ++(index -> support)[u][k], --(index -> support)[w][k];
        }
    
        for (auto u: waitList) {
            trpOrd[k].first -> Delete(trpId[u][k], trpOrd[k].second);
            assert((G -> core)[u] >= k);
            (index -> support)[u][k] = deg[u];
            trpId[u][k] = (trpOrd[k].first -> TreeNode(trpId[u][k], Node(u, G -> wei[u]), mxKey[u] >= k));
            trpOrd[k].first -> InsertAfter(trpId[u][k], insPrev[u], trpOrd[k].second);
            inW[u] = false;
        }
    
        for (auto u: P) {
            cover[u] = false;
            trpOrd[k].first -> Delete(trpId[u][k], trpOrd[k].second);
            assert((G -> core)[u] < k);
            (index -> support)[u].pop_back();
            trpId[u].pop_back();
            (index -> pre)[u].pop_back();
        }
    }
    
    void OrdDynamic::remove(int u, int v) {
        std::vector<std::vector<int> > core_dec((index -> k_max) + 1, std::vector<int>());
        cm -> RemoveRem(u, v, G -> nei, G -> core, core_dec);

        for (int k = (index -> k_max); k >= 1; --k) {
            clock_t START = clock();
            removeUpdate(k, core_dec[k], u, v);
            clock_t END = clock();
            ACL_TIME[k] += END - START;
            if (trpOrd[k].second == 0) {
                trpOrd.pop_back();
            }
        }
    }

    void OrdDynamic::query(int k, int r, std::vector<std::pair<double, std::vector<int> > > &ans){
        ans.clear();
        std::vector<std::pair<int, int> > V;
        trpOrd[k].first -> QueryCV(trpOrd[k].second, r, V);
        index -> formAnswer(V, tree, ans);
        return;
    }
    
    void OrdDynamic::saveIndex(){
        (index -> cv).resize((index -> k_max) + 1);
        for (int k = 1; k <= (index -> k_max); ++k) {
            (index -> cv)[k].clear();
            std::vector<int> tmp;
            trpOrd[k].first -> TraverseToVec(trpOrd[k].second, tmp);
            std::vector<int> S;
            for (auto u: tmp) {
                if (mxKey[u] >= k && S.size()) {
                    (index -> cv)[k].push_back(S);
                    S.clear();
                }
                S.push_back(u);
            }
            if (S.size()) (index -> cv)[k].push_back(S);
        }
    }
}