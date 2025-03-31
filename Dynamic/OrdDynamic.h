#ifndef ORDDYNAMIC_H
#define ORDDYNAMIC_H

#include "Dynamic.h"
#include "../Index/OrdIndex.h"

namespace ORDDYNAMIC {
    class OrdDynamic: public DYNAMIC::Dynamic {
    public:

        std::vector<int> deg;
        std::vector<int> mxKey;
        std::vector<int> insPrev;
        std::vector<std::vector<int> > trpId;
        std::vector<std::pair<TREAP2::Treap*, int> > trpOrd;

        
        OrdDynamic();
        ~OrdDynamic();

        void init(const char* path) override;
        void query(int u, int k, std::vector<std::pair<double, std::vector<int> > > &ans) override;
        void insert(int u, int v) override;
        void remove(int u, int v) override;
        void saveIndex() override;

        int findPrev(int k, int u);
        int getRank(int k, int u);
        bool checkToQue(int w, int k);
        bool checkToQue2(int w, int v, int x, int y, int k);
        void removeUpdate(int k, std::vector<int> waitList, int x, int y);
        bool removeShrink(int k, int s, std::unordered_set<int> &P, int &prev, std::vector<int> &waitList);

        void insertUpdate(int k, std::vector<int> waitList, int x, int y);
        int findNext(int k, Node x);
        void insertShrink(int k, int s, std::set<Node> &P, int &prev, std::vector<int> &waitList);

        void printIndex(int k);
    } ;
}

#endif