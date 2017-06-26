#include "BasicInfo.h"
#include "DumpXml.h"
#include "Traj.h"
#include "MolDataStruct.h"
#include "AnalFunc.h"

#ifndef __REALCOORCM_H__
#define __REALCOORCM_H__

class Traj2CMBin : public AnalFunc
{
    public:
        Traj2CMBin(boost::shared_ptr<Traj> traj, boost::shared_ptr<BasicInfo> bi , const string & fn);
        Traj2CMBin(boost::shared_ptr<Traj> traj, boost::shared_ptr<BasicInfo> bi);
        void compute();
        void finish();
    private:
        boost::shared_ptr<BasicInfo> m_basicinfo;
        boost::shared_ptr<Traj> m_traj;
        vector<vec> m_pos_cm;
        string m_output;


        //  void DumpTraj2CMBin();

};


void export_Traj2CMBin();

#endif
