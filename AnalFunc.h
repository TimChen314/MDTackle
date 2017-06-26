#include "MolDataStruct.h"  //Include vec class
#include "Traj.h"
#include <cuda_runtime.h>
#include <boost/python.hpp>

#ifndef __ANALFUNC_H__
#define __ANALFUNC_H__

class AnalFunc : boost::noncopyable
{
    protected:
        AnalFunc() {};
        boost::shared_ptr<Traj> traj;
    public:
        virtual void compute() {};
        virtual void finish() {};

};

void export_AnalFunc();

#endif
