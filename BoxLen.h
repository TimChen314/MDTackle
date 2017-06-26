#include "Traj.h"
#include "MolDataStruct.h"
#include "AnalFunc.h"

#ifndef __DENSITY_H__
#define __DENSITY_H__

class BoxLen : public AnalFunc
{
    public:
        BoxLen(boost::shared_ptr<Traj> traj );
        void compute();
        void finish();
    private:
        ofstream m_file;
        vec sum_boxlen;
        vector<vec> boxlen;


        void DumpTailBoxLen();

};


void export_BoxLen();

#endif
