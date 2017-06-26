#include "RgByType.h"

#include <sys/time.h>
using namespace boost::python;
using namespace boost;

unsigned RgByType::n_compute = 0;

RgByType::RgByType(boost::shared_ptr<XmlReader> read, boost::shared_ptr<Traj> tj, const string& type, unsigned len)
{
    build = read;
    traj=tj;
    rg_type = type;
    length = len;
    vector<double> tmp_vecf1(traj->getNfram_remain(),0.0);
}


void RgByType::finish() 
{
    DumpRgByType();
}

void RgByType::DumpRgByType()
{
    cout<<"RgByType : Begin to dump !"<<endl;
    ostringstream ss;
    ss << "RgByType_"   << rg_type <<  ".dat" ;
    ofstream h(ss.str().c_str());
    h<<"#From MDTackle.so by Tim Chen:"<<endl;
    if(!rg.size()) {
        cerr<<"***RgByType:: output array is empty***";
        throw runtime_error("RgByType:: output array is empty");
    }
    for(unsigned i=0;i<rg.size();i++)  {
        h<<rg[i] << " " << ree[i] << endl;
    }
    double aver_rg = accumulate(rg.begin(), rg.end(), 0.0) / rg.size();
    double aver_ree = accumulate(ree.begin(), ree.end(), 0.0) / ree.size();
    h << "# average rg and ree is : "  << aver_rg << " " << aver_ree << endl;
    h.close();
}

void RgByType::compute()
{
    n_compute ++ ;
    vector<vec> pos = traj->getPos();
    unsigned int_type = build->getTypeID(rg_type);
    vector<unsigned> m_type = build -> getIdxType();
    vector<vec> target_pos;
    for(unsigned i=0;i!=pos.size(); ++i)
    {
        if(m_type[i]!=int_type)
            continue;
        target_pos.push_back(pos[i]);
    }
    if(n_compute==1)
        cout << "Info: number of beads of type "  << rg_type  << " is "  << target_pos.size()<< endl;
    if(target_pos.size()%length!=0)
    {
        cerr<<"***RgByType:: number of beads of type " << rg_type << " should be divided by chain length ("<<length<<")***";
        throw runtime_error("RgByType:: number of beads of type should be divided by chain length ");
    }

    unsigned n_mol = target_pos.size() / length;
    unsigned idx = 0;
    vec ree2_fram(0, 0, 0);
    vec rg2_fram(0, 0, 0);
    for(unsigned i=0;i!=n_mol;++i)
    {
        vec ree2_tmp(0, 0, 0);
        ree2_tmp  = (target_pos[idx + length - 1] - target_pos[idx]);
        ree2_fram += ree2_tmp * ree2_tmp;

        vec cm(0, 0, 0);
        vec dd;
        vec rg2_tmp(0, 0, 0);
        for(unsigned j=0;j!=length; ++j)
        {
            cm += target_pos[idx];
            ++ idx;
        }
        cm /= length;
        idx -= length;
        for(unsigned j=0;j!=length; ++j)
        {
            dd = target_pos[idx] - cm;
            rg2_tmp += dd * dd;
            ++idx;
        }
        rg2_tmp /= length;
        rg2_fram += rg2_tmp;
    }
    ree2_fram /= n_mol;
    rg2_fram /= n_mol;

    double ree_fram = sqrt(ree2_fram.x + ree2_fram.y + ree2_fram.z);
    double rg_fram = sqrt(rg2_fram.x + rg2_fram.y + rg2_fram.z);

    ree.push_back(ree_fram);
    rg.push_back(rg_fram);
}

/* 
   void RgByType::compute()
   {
   n_compute++;

   vector<vec> pos_np;
   vector<vec> pos = traj->getPos();
   QuestRealCoorCMAndSetToOri( pos, pos_np, basicinfo->nm_tot, basicinfo->MolNA);
   pos_np.resize(basicinfo -> nm[0] );

   unsigned int idx=basicinfo->nm[0] * basicinfo->MolNA[0];
   unsigned int bin;
   vec box = traj->getBox();
   vec boxINV(vec(1.0,1.0,1.0)/box);
   double dr,deltr=cutoff/maxbin;
   for(unsigned int i=basicinfo->nm[0];i!=basicinfo->nm_tot;i++) {
   unsigned mol_type=basicinfo->MolType[i]-1;//mol_type is type index without NP.
   for(unsigned int k=0;k!=basicinfo->MolNA[i]-seglen+1;k++) {
//testdump: cout<<i<<" "<<k<<" "<<basicinfo->nm[0]<<" ha0 "<<basicinfo->MolNA[i]<<endl;
vec pos_c(0.0,0.0,0.0);
vec d(0.0,0.0,0.0);
double rg=0.0;

for(unsigned int j=0;j!=seglen;j++) {
pos_c += pos[idx];
idx++;
}
idx -= seglen;
pos_c /= double(seglen);

for(unsigned int j=0;j!=seglen;j++) {
d=pos_c-pos[idx];
rg += ScalarProduct(d,d);
idx++;
}
rg /= double(seglen);
for(unsigned int j=0;j!=basicinfo->nm[0];j++) {
d = pos_c - pos_np[j];
pbc(box,boxINV,d);
dr = sqrt(ScalarProduct(d,d));
bin=int(dr/deltr);
if(bin<maxbin) {
//testdump: cout<<"ha "<<mol_type<<" "<<bin<<" "<<rg<<endl;
rg_dis[mol_type][bin] += rg;
rg_dis_n[mol_type][bin] += 1;
}
}
idx = idx - seglen + 1 ;
}
idx = idx + seglen - 1 ;
}

}
*/ 


void export_RgByType()
{
    class_<RgByType, bases<AnalFunc> , boost::noncopyable>
        ("RgByType", init<boost::shared_ptr<XmlReader> , boost::shared_ptr<Traj>, const string &, unsigned >() )
        ;
}

