#include "Traj2CMBin.h"
using namespace boost::python;
using namespace boost;

Traj2CMBin::Traj2CMBin(boost::shared_ptr<Traj> tj , boost::shared_ptr<BasicInfo> bi): m_output("all_CM.bin")
{
    m_traj=tj;
    m_basicinfo=bi;
}

Traj2CMBin::Traj2CMBin(boost::shared_ptr<Traj> tj , boost::shared_ptr<BasicInfo> bi, const string& dump )
{
    m_traj=tj;
    m_basicinfo=bi;
    m_output = dump;
}


void Traj2CMBin::compute()
{
    QuestRealCoorCMAndSetToOri( m_traj->getPos(), m_pos_cm, m_basicinfo->nm_tot, m_basicinfo->MolNA  );
}

void Traj2CMBin::finish()
{
    unsigned timestep = 0;
    vector < vec > vel;
    DumpAs3202Version(m_output, m_pos_cm, vel, m_traj ->getOnceBox(), timestep, (vel.size()!=0));
}


void export_Traj2CMBin()
{
    class_<Traj2CMBin, bases<AnalFunc> , boost::noncopyable>
        ("Traj2CMBin", init<boost::shared_ptr<Traj>, boost::shared_ptr<BasicInfo> , const string& >() )
        .def(init<boost::shared_ptr<Traj>, boost::shared_ptr<BasicInfo> >() )
        ;
}


