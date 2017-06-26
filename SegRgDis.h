#include "BasicInfo.h"
#include "Traj.h"
#include "MolDataStruct.h"
#include "AnalFunc.h"

#ifndef __SEGRGDIS_H__
#define __SEGRGDIS_H__


class SegRgDis : public AnalFunc
{
    public:
        SegRgDis(boost::shared_ptr<Traj> traj, boost::shared_ptr<BasicInfo> basicinfo, unsigned seglen);
        void compute();
        void finish();
    private:
        static unsigned n_compute;
        double cutoff;
        unsigned int maxbin,seglen;
        boost::shared_ptr<BasicInfo> basicinfo;
        vector< vector<double> > rg_dis;
        vector< vector<unsigned long int> > rg_dis_n;
        void OnceInit();

        void DumpSegRgDis();

};


void export_SegRgDis();

#endif
