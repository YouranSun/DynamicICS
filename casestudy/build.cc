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

int main(int argc, char **argv){

    // argv[1]: vertex path
    // argv[2]: edge path
    // argv[3]: index path
    // argv[4]: shift
    
    char *path_v = argv[1];
    char *path_e = argv[2];
    char *path_index = argv[3];
    int shift = atoi(argv[4]);

    INDEX::Index *index = new ORDINDEX::OrdIndex();

    index -> G = new Graph();
    (index -> G) -> readCsv(path_v, path_e, shift);
    KEEPTIME("Finished reading: ", stamp);

    core::GLIST *cm = new core::GLIST((index -> G) -> n);
    cm -> ComputeCore((index -> G) -> nei, true, (index -> G) -> core);
    KEEPTIME("Finished core number: ", stamp);

    index -> computeIndex();
    write(path_index, index);
    KEEPTIME("Index init: ", stamp);

    delete cm;
    delete index;
    return 0;
}