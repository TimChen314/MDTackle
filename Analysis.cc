#include "Analysis.h"
#include <iostream>
#include <sys/time.h>

using namespace  std;

Analysis::Analysis(boost::shared_ptr<Traj> tj)
{
    traj=tj;
}

void Analysis::add(boost::shared_ptr<AnalFunc> af)
{
    m_analfunc.push_back(af);
}

void Analysis::run()
{
    timeval start;
    timeval end;

    unsigned nfram_remain=traj->getNfram_remain();
    for(unsigned i=0;i!=nfram_remain;++i) {
        cout<<(100.0*i/nfram_remain) <<"%"<<endl;
        if(i%100==0) 
            gettimeofday(&start,NULL);
        vector<vec> pos;
        vec box;
        traj->getRealCoor(i);
        for(unsigned j=0;j!=m_analfunc.size();++j) {
            m_analfunc[j]->compute();
        }
        if(i%10==0) { 
            gettimeofday(&end,NULL);
            long timeusr=(end.tv_sec-start.tv_sec)*1000000+(end.tv_usec-start.tv_usec);
            cout<<"Time(ms) "<<timeusr/1000.0<<endl;
        }
    }

    for(unsigned j=0;j!=m_analfunc.size();++j) {
        m_analfunc[j]->finish();
        cout<<"+----------------------------------------------------------------------+"<<endl;
    }

}


using namespace boost::python;
void export_Analysis()
{
    class_<Analysis, boost::shared_ptr<Analysis>, boost::noncopyable>
        ("Analysis",init<boost::shared_ptr<Traj> >() )
        .def("add", &Analysis::add)
        .def("run", &Analysis::run)
        ;
}

