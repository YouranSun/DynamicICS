#include <bits/stdc++.h>
#include "Dynamic.h"

namespace DYNAMIC{

Dynamic::Dynamic(): cm(nullptr), index(nullptr), G(nullptr){
    return;
}


Dynamic::~Dynamic(){
    delete index;
    delete cm;
    return;
}


void Dynamic::checkCore(){
    core::CoreMaintenance *tmpcm = new core::GLIST(G -> n);
    std::vector<int> tmpcore(G -> n);
    tmpcm -> ComputeCore(G -> nei, true, tmpcore);
    for (int i = 0; i < (G -> n); ++i) {
        eprintf("%d ", tmpcore[i]);
        assert(tmpcore[i] == (G -> core)[i]);
    }
    eputs("");
    delete tmpcm;
}

void Dynamic::GET_TIME(double *TT){
    for (int k = 1; k < ACL_TIME.size(); ++k) {
        eprintf("k = %d %d\n", k, ACL_TIME[k]);
        if (k <= ceil((ACL_TIME.size() - 1) * 0.25)) {
            TT[0] += ACL_TIME[k];
        }
        else if (k <= ceil((ACL_TIME.size() - 1) * 0.5)) {
            TT[1] += ACL_TIME[k];
        }
        else if (k <= ceil((ACL_TIME.size() - 1) * 0.75)) {
            TT[2] += ACL_TIME[k];
        }
        else {
            TT[3] += ACL_TIME[k];
        }
    }
    eputs("");
    return;
}

}