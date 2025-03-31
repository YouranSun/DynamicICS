#ifndef ICPINDEX_H
#define ICPINDEX_H

#include "Index.h"

namespace ICPINDEX{
    class ICPIndex : public INDEX::Index{
    public:
        ICPIndex();
        ~ICPIndex();

        void construct() override;
        void readIndex(FILE *save) override;
        void writeIndex(FILE *save) override;
        
    } ;
}

#endif