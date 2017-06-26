#include "BoxLen.h"
using namespace boost::python;
using namespace boost;


BoxLen::BoxLen(boost::shared_ptr<Traj> tj ) 
{
    traj=tj;
    m_file.open("boxlength.dat");
}

void BoxLen::compute()
{
    vec tmp(traj->getBox());
    sum_boxlen += tmp;
    //test    cout<<tmp.x<<" "<<tmp.y<<" "<<tmp.z<<"\n";
    m_file<<tmp.x<<" "<<tmp.y<<" "<<tmp.z<<"\n";
}

void BoxLen::finish()
{
    sum_boxlen /= double(traj->getNfram_remain());
    DumpTailBoxLen();
}

void BoxLen::DumpTailBoxLen()
{
    vector<string> fn=traj->getFilename();
    m_file<<"#Box length of "<<fn[0]<<" to "<<fn[traj->getNfram_remain()-1]<<endl;
    m_file<<"#The average box length is "<<sum_boxlen.x<<" "<<sum_boxlen.y<<" "<<sum_boxlen.z<<endl;
    m_file << "#Number Density is : " << traj -> getNtot() / (sum_boxlen.x * sum_boxlen.y * sum_boxlen.z) << endl;
    cout<<"BoxLen : Finish to dump boxlength.dat"<<endl;
}

void export_BoxLen()
{
    class_<BoxLen, bases<AnalFunc> , boost::noncopyable>
        ("BoxLen", init<boost::shared_ptr<Traj> >() )
        .def("compute",&BoxLen::compute)
        .def("finish",&BoxLen::finish)
        ;
}


