#include "BasicInfo.h"
#include "Traj.h"
#include "MolDataStruct.h"
#include "AnalFunc.h"
#include "particles/XmlReader.h"

#ifndef __RGBYTYPE_H__
#define __RGBYTYPE_H__


class RgByType : public AnalFunc
{
    public:
        RgByType(boost::shared_ptr<XmlReader> read, boost::shared_ptr<Traj> traj, const string & type, unsigned length);
        boost::shared_ptr<XmlReader> build;

        void compute();
        void finish();

    private:
        static unsigned n_compute;
        boost::shared_ptr<BasicInfo> basicinfo;
        vector<double> rg;
        vector<double> ree;
        string rg_type;
        unsigned length = 0;

        void DumpRgByType();

};


void export_RgByType();

#endif
