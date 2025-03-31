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

inline void writeAns(int k, int r, FILE* file, std::vector<std::pair<double, std::vector<int> > > ans) {
    // eprintf("ok\n");
    fprintf(file, "k = %d r = %d\n", k, r);
    for (auto [val, S]: ans){
        fprintf(file, "%lf ", val);
        for (auto u: S) fprintf(file, "%d ", u);
        fprintf(file, "\n");
    }
}

int main(int argc, char **argv){

    // argv[1]: index path
    // argv[2]: query path
    // argv[3]: out path
    
    char *path_index = argv[1];
    char *path_query = argv[2];
    char *path_out = argv[3];

    DYNAMIC::Dynamic *dynamic = nullptr;
    dynamic = new ORDDYNAMIC::OrdDynamic();

    dynamic -> init(path_index);

    FILE *parafile = fopen(path_query, "r");
    FILE *outfile = fopen(path_out, "w");

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
    return 0;
}