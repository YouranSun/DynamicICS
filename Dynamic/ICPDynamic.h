#ifndef ICPDYNAMIC_H
#define ICPDYNAMIC_H

#include "Dynamic.h"
#include "../Index/ICPIndex.h"

namespace ICPDYNAMIC {
    class ICPDynamic : public DYNAMIC::Dynamic {
    public:
        const int RECOMPUTE = false;

        std::vector<int> L;
        std::vector<std::vector<int> > R;

        ICPDynamic();
        ~ICPDynamic();

        void init(const char* path) override;
        void query(int u, int k, std::vector<std::pair<double, std::vector<int> > > &ans) override;
        void insert(int u, int v) override;
        void remove(int u, int v) override;
        void saveIndex() override;

        void queryDfs(int i, int k, std::vector<int> &S);

        bool isRecompute(int u, int v, int k);
        void recompute(int k);
        void insertionDFS(int u, int k, int R_min);

    } ;
}

#endif