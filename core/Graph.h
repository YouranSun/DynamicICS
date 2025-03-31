#ifndef GRAPH_H_
#define GRAPH_H_

#include <bits/stdc++.h>
#include <ctime>
#include "../gadget/Tools.h"

class Graph{
public:
    int n, m, m2;
    std::vector<int> wei;
    std::vector<double> wei_ori;
    std::vector<int> ord;
    std::vector<int> rnk;
    std::vector<int> seq;
    std::vector<int> core;
    std::vector<std::vector<int> > nei;
    std::vector<std::pair<int, int> > edges;

    long long mn_quarter;
    long long mx_quarter;

    Graph(const int &m2_ = 0){
        m2 = m2_;
        n = m = 0;
        wei = std::vector<int>();
        ord = std::vector<int>();
        rnk = std::vector<int>();
        core = std::vector<int>();
        nei = std::vector<std::vector<int> >();
        edges = std::vector<std::pair<int, int> >();
        return;
    }
    ~Graph(){}

    void init(){
        wei.resize(n, 0);
        ord.resize(n, 0);
        rnk.resize(n, 0);
        seq.resize(n, 0);
        core.resize(n, 0);
        nei.resize(n);
        edges.resize(m);
        return;
    }

    void pageRank(){
        std::vector<double> tmp(n, 1.0);
        const int ITER_TIME = 200;
        for (int T = 0; T < ITER_TIME; ++T) {
            std::vector<double> ntmp(n, 0.0);
            for (int u = 0; u < n; ++u) {
                if (nei[u].size()){
                    for (auto v: nei[u]) ntmp[u] += tmp[v];
                    ntmp[u] /= nei[u].size();
                }
            }
            tmp.swap(ntmp);
        }
        std::sort(ord.begin(), ord.end(), [&](const int &i, const int &j) { return tmp[i] < tmp[j]; });
        for (int i = 0; i < n; ++i) wei[ord[i]] = i;
    }

    void read(const char *path){
        eprintf("Graph Reading Start path = %s\n", path);
        std::set<std::pair<int, int> > deduplicate;
        FILE *file = fopen(path, "r");

        fscanf(file, "%d %d", &n, &m);
        init();

        for (int i = 0; i < m; ++i) {
            int u, v;
            fscanf(file, "%d %d", &u, &v);
            if (u == v) continue;
            if (u > v) std::swap(u, v);
            deduplicate.insert(std::make_pair(u, v));
        }

        for (int i = 0; i < n; ++i) wei[i] = i, ord[i] = i;
        pageRank();
        fclose(file);

        
        // Split into the original edge set and insertion
        m = deduplicate.size();
        // m -= m2;

        int mptr = 0;
        for(auto it: deduplicate) {
            int u = it.first, v = it.second;
            edges[mptr++] = std::make_pair(u, v);
        }
        // std::mt19937 rnd(2226701);
        // std::shuffle(edges.begin(), edges.end(), rnd);

        for (int i = 0; i < m - m2; ++i) {
            int u = edges[i].first, v = edges[i].second;
            nei[u].push_back(v), nei[v].push_back(u);
        }
        m -= m2;
    }
    
    void readCsv(const char *path_v, const char *path_e, int shift) {
        eprintf("%s\n%s\n", path_v, path_e);
        std::ifstream in_v(path_v);
        std::ifstream in_e(path_e);
        std::string line;

        std::map<long long, double> weights;
        while (std::getline(in_v, line)) {
            long long vid;
            double vweight;
            sscanf(line.c_str(), "%lld,%lf", &vid, &vweight);
            weights[vid] = vweight;
        }

        std::vector<long long> vids_in_e;
        std::vector<std::pair<std::pair<long long, long long>, long long> > vid_of_edges;

        mn_quarter = 1e18;
        mx_quarter = 0;

        while(std::getline(in_e, line)) {
            long long uid, vid;
            time_t timestamp;
            sscanf(line.c_str(), "%lld,%lld,%lu", &uid, &vid, &timestamp);
            struct tm timeinfo;
            gmtime_r(&timestamp, &timeinfo);
            long long quarter = timeinfo.tm_year * 12 + timeinfo.tm_mon;
            mn_quarter = std::min(mn_quarter, quarter);
            mx_quarter = std::max(mx_quarter, quarter);
            vid_of_edges.push_back(std::make_pair(std::make_pair(uid, vid), quarter));
            vids_in_e.push_back(uid);
            vids_in_e.push_back(vid);
        }

        sort(vids_in_e.begin(), vids_in_e.end(), [&] (const long long &u, const long long &v) {
            return weights[u] < weights[v];
        });

        vids_in_e.erase(std::unique(vids_in_e.begin(), vids_in_e.end()), vids_in_e.end());
        
        n = vids_in_e.size();
        m = vid_of_edges.size();
        init();

        wei_ori = std::vector<double>(n);

        std::map<long long, int> vid_to_index;
        for (int i = 0; i < vids_in_e.size(); ++i) {
            vid_to_index[vids_in_e[i]] = i;
            ord[i] = i;
            wei[i] = i;
            wei_ori[i] = weights[vids_in_e[i]];
        }

        long long section = (mx_quarter - mn_quarter + 1) / 4;
        mn_quarter += section * shift;
        mx_quarter += section * (-2 + shift);

        std::set<std::pair<int, int> > deduplicate;
        for (auto [e, w]: vid_of_edges) {
            auto [u, v] = e;
            if (w >= mn_quarter && w <= mx_quarter) {
                deduplicate.insert(std::make_pair(vid_to_index[u], vid_to_index[v]));
            }
        }

        m = deduplicate.size();
        // m -= m2;

        int mptr = 0;
        for(auto it: deduplicate) {
            int u = it.first, v = it.second;
            edges[mptr++] = std::make_pair(u, v);
        }
        // std::mt19937 rnd(2226701);
        // std::shuffle(edges.begin(), edges.end(), rnd);

        for (int i = 0; i < m - m2; ++i) {
            int u = edges[i].first, v = edges[i].second;
            nei[u].push_back(v), nei[v].push_back(u);
        }
        m -= m2;
        in_v.close();
        in_e.close();
    }

    void readIndex(FILE *file){
        eprintf("Graph ReadIndex Start\n");
        fscanf(file, "%d%d%d", &n, &m, &m2);
        eprintf("n = %d m = %d m2 = %d\n", n, m, m2);

        nei.resize(n);
        core.resize(n);
        wei.resize(n);
        wei_ori.resize(n);
        for (int i = 0; i < n; ++i) fscanf(file, "%d", &wei[i]);
        edges.resize(m + m2);
        for (int i = 0; i < m + m2; ++i){
            int u, v;
            fscanf(file, "%d%d", &u, &v);
            edges[i] = std::make_pair(u, v);
            if (i < m) nei[u].push_back(v), nei[v].push_back(u);
        }
        for (int i = 0; i < n; ++i) fscanf(file, "%lf ", &wei_ori[i]);

        eprintf("Graph ReadIndex End\n");
        return;
    }

    void write(FILE *file) {
        fprintf(file, "%d %d %d\n", n, m, m2);
        for (int i = 0; i < n; ++i) fprintf(file, "%d ", wei[i]); fputs("\n", file);
        for (int i = 0; i < m + m2; ++i) fprintf(file, "%d %d\n", edges[i].first, edges[i].second);
        for (int i = 0; i < n; ++i) fprintf(file, "%lf ", wei_ori[i]); fputs("\n", file);
        return;
    }
} ;

#endif