#include <unistd.h>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <map>

#include "../gadget/Tools.h"

#include "../core/Graph.h"
#include "../core/core.h"
#include "../core/glist.h"

#include "../Index/Index.h"
#include "../Index/ICPIndex.h"
#include "../Index/OrdIndex.h"

#include "../Dynamic/Dynamic.h"
#include "../Dynamic/ICPDynamic.h"
#include "../Dynamic/OrdDynamic.h"

clock_t stamp;

inline void write(char *savepath, INDEX::Index *index){
    eprintf("write savepath = %s", savepath);
    FILE *save = fopen(savepath, "w");
    index -> writeIndex(save);
    fclose(save);
}

inline void build(char *inpath, char *savepath, const int m2, const int baseline){
    INDEX::Index *index = nullptr;
    if (baseline) {
        index = new ICPINDEX::ICPIndex();
    }
    else {
        index = new ORDINDEX::OrdIndex();
    }

    index -> G = new Graph(m2);
    (index -> G) -> read(inpath);
    KEEPTIME("Finished reading: ", stamp);

    core::GLIST *cm = new core::GLIST((index -> G) -> n);
    cm -> ComputeCore((index -> G) -> nei, true, (index -> G) -> core);
    KEEPTIME("Finished core number: ", stamp);

    index -> computeIndex();
    write(savepath, index);
    KEEPTIME("Index init: ", stamp);

    delete cm;
    delete index;
    return;
}

inline void testInsert(char *savepath, char *outpath, int m2, const int baseline, char *logpath) {
    FILE *logfile = fopen(logpath, "w");

    DYNAMIC::Dynamic *dynamic = nullptr;
    eprintf("baseline = %d\n", baseline);
    if (baseline) {
        dynamic = new ICPDYNAMIC::ICPDynamic();
    }
    else {
        dynamic = new ORDDYNAMIC::OrdDynamic();
    }
    dynamic -> init(savepath);
    // for (int u = 0; u < (km -> G) -> n; ++u) eprintf("%d ", (km -> G) -> wei[u]); eputs("");

    clock_t START_INSERT = clock(), END_INSERT = clock();
    KEEPTIME("Start Testing Insert: ", stamp);
    Graph *G = (dynamic -> index) -> G;
    int ORI = (G -> m);
    for (int &i = (G -> m), &j = (G -> m2); j && m2; ++i, --j, --m2) {
        if ((i - ORI) % 2000 == 0) {
            END_INSERT = clock();
            fprintf(logfile, "%d %lf %lld\n", i - ORI, 1.0 * (END_INSERT - START_INSERT) / CLOCKS_PER_SEC, dynamic -> candSize);
            eprintf("Edge cnt: %d", i - ORI), KEEPTIME(" Time: ", stamp);
        }
        dynamic -> insert(G -> edges[i].first, G -> edges[i].second);        
    }
    KEEPTIME("End Testing Insert: ", stamp);
    eprintf("Total Edge: %d Total Time: %lf\n", (G -> m) - ORI, 1.0 * (END_INSERT - START_INSERT) / CLOCKS_PER_SEC);

    double TT[4] = {0, 0, 0, 0};
    dynamic -> GET_TIME(TT);
    fprintf(logfile, "%lf %lf %lf %lf\n", TT[0], TT[1], TT[2], TT[3]);

    // km -> saveIndex();
    // write(outpath, km -> index, baseline);
    delete dynamic;
    fclose(logfile);
    return;
}

inline void testRemove(char *savepath, char *outpath, int m2, const int baseline, char *logpath) {
    FILE *logfile = fopen(logpath, "w");

    DYNAMIC::Dynamic *dynamic = nullptr;
    if (baseline) {
        dynamic = new ICPDYNAMIC::ICPDynamic();
    }
    else {
        dynamic = new ORDDYNAMIC::OrdDynamic();
    }
    dynamic -> init(savepath);
    // for (int u = 0; u < (km -> G) -> n; ++u) eprintf("%d ", (km -> G) -> wei[u]); eputs("");

    clock_t START_REMOVE = clock(), END_REMOVE = clock();
    KEEPTIME("Start Testing Remove: ", stamp);
    Graph *G = (dynamic -> index) -> G;
    int ORI = (G -> m);
    for (int &i = (G -> m), &j = (G -> m2), ptr = (G -> m) + (G -> m2); i && m2; ++j, --i, --ptr, --m2) {
        if ((ORI - i) % 2000 == 0) {
            END_REMOVE = clock();
            fprintf(logfile, "%d %lf %lld\n", ORI - i, 1.0 * (END_REMOVE - START_REMOVE) / CLOCKS_PER_SEC, dynamic -> candSize);
            eprintf("Edge cnt: %d ", ORI - i), KEEPTIME(" Time: ", stamp);
        }
        dynamic -> remove(G -> edges[ptr - 1].first, G -> edges[ptr - 1].second);        
        if ((ORI - i) > 12000) break;
    }
    KEEPTIME("End Testing Delete: ", stamp);
    eprintf("Total Edge: %d Total Time: %lf\n", (G -> m2), 1.0 * (END_REMOVE - START_REMOVE) / CLOCKS_PER_SEC);

    double TT[4] = {0, 0, 0, 0};
    dynamic -> GET_TIME(TT);
    fprintf(logfile, "%lf %lf %lf %lf\n", TT[0], TT[1], TT[2], TT[3]);

    // km -> saveIndex();
    // write(outpath, km -> index, baseline);
    delete dynamic;
    fclose(logfile);
    return;
}


inline void writeAns(int k, int r, FILE* file, std::vector<std::pair<double, std::vector<int> > > ans) {
    // eprintf("ok\n");
    fprintf(file, "k = %d r = %d\n", k, r);
    for (auto [val, S]: ans){
        fprintf(file, "%lf ", val);
        for (auto u: S) fprintf(file, "%d ", u);
        fprintf(file, "\n");
    }
}

inline void testQuery(char *savepath, char *outpath, char *parapath, const int baseline){
    DYNAMIC::Dynamic *dynamic = nullptr;
    if (baseline) {
        dynamic = new ICPDYNAMIC::ICPDynamic();
    }
    else {
        dynamic = new ORDDYNAMIC::OrdDynamic();
    }
    dynamic -> init(savepath);

    FILE *parafile = fopen(parapath, "r");
    FILE *outfile = fopen(outpath, "w");

    int q; fscanf(parafile, "%d", &q);
    for (int i = 0; i < q; ++i) {
        int k, r; fscanf(parafile, "%d%d", &k, &r);
        std::vector<std::pair<double, std::vector<int> > > ans;
        k = std::min(k, (dynamic -> index) -> k_max);
        KEEPTIME("start query = ", stamp);
        dynamic -> query(k, r, ans);
        writeAns(k, r, outfile, ans);
        eprintf("k = %d r = %d\n", k, r);
        KEEPTIME("query time = ", stamp);
    }

    delete dynamic;

    fclose(parafile);
    fclose(outfile);
}

int main(int argc, char **argv){

    stamp = clock();

    // argv[1]: inpath
    // argv[2]: outpath
    // argv[3]: m2
    // argv[4]: b/i/d/q
    // argv[5]: baseline/advanced
    // argv[6]: parapath
    // argv[6]: log path

    int m2 = 0; sscanf(argv[3], "%d", &m2);
    char *inpath = argv[1], *outpath = argv[2], *type = argv[4], *logpath = argv[6];
    int baseline = (argv[5][0] == 'b' ? 1 : (argv[5][0] == 'a' ? 0 : -1));

    eprintf("%s\n", inpath);
    eprintf("%s\n", outpath);
    eprintf("%s\n", type);

    if(*type == 'b') build(inpath, outpath, m2, baseline);
    else if (*type == 'i') testInsert(inpath, outpath, m2, baseline, logpath);
    else if (*type == 'd') testRemove(inpath, outpath, m2, baseline, logpath);
    else {
        char *parapath = argv[6];
        testQuery(inpath, outpath, parapath, baseline);
    }

    return 0;
}