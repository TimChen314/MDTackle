#include "AnalFunc.h"
#include "Traj.h"

class Analysis : boost::noncopyable  
{
    public:
        Analysis(boost::shared_ptr<Traj> tj);
        void add(boost::shared_ptr<AnalFunc> af);
        void run();

    private:
        boost::shared_ptr<Traj> traj;
        vector< boost::shared_ptr<AnalFunc> > m_analfunc;
        timeval start;
        timeval end;

};

void export_Analysis();
