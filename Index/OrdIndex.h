#ifndef ORDINDEX_H
#define ORDINDEX_H

#include "Index.h"
#include "../gadget/Tools.h"

namespace ORDINDEX {
    class OrdIndex: public INDEX::Index {
    public:
        OrdIndex();
        ~OrdIndex();

        void construct() override;
        void readIndex(FILE *save) override;
        void writeIndex(FILE *save) override;
    } ;
}

#endif