#include "SegRgDis.h"

#include <sys/time.h>
using namespace boost::python;
using namespace boost;

unsigned SegRgDis::n_compute = 0;

SegRgDis::SegRgDis(boost::shared_ptr<Traj> tj, boost::shared_ptr<BasicInfo> bi, unsigned sl)
{
    traj=tj;
    basicinfo=bi;
    seglen = sl;
    OnceInit();

    vector<double> tmp_vecf1(maxbin,0.0);
    vector<unsigned long int> tmp_veci1(maxbin,0);
    for(unsigned i=0;i!=basicinfo->ntype-1; ++i) {
        rg_dis.push_back(tmp_vecf1);
        rg_dis_n.push_back(tmp_veci1);
    }
}


void SegRgDis::finish() 
{
    DumpSegRgDis();
}

void SegRgDis::DumpSegRgDis()
{
    cout<<"SegRgDis : Begin to dump !"<<endl;
    double bin=cutoff/(double)maxbin;
    for(unsigned i=0;i<basicinfo->ntype-1;i++)  {
        string tmp_str1 = "Seg";
        string tmp_str2 = "Rg_dis_";
        string tmp_str3 = ".dat";
        ostringstream dump_name;
        dump_name<<tmp_str1<<seglen<<tmp_str2<<i<<tmp_str3;
        ofstream h(dump_name.str().c_str());
        h<<"#From MDTackle.so by Tim Chen:"<<endl;
        if(!rg_dis[i].size()) {
            cerr<<"***SegRgDis:: output array is empty***";
            throw runtime_error("SegRgDis:: output array is empty");
        }
        for(unsigned j=0;j<rg_dis[i].size();j++)  {
            rg_dis[i][j] /= rg_dis_n[i][j];
            rg_dis[i][j] = sqrt(rg_dis[i][j]);
            double r=(j+0.5)*bin;
            h<<r<<" "<<rg_dis[i][j]<<" "<<rg_dis_n[i][j]<<endl;
        }
        h.close();
    }
}

void SegRgDis::OnceInit()
{
    vec box=traj->getOnceBox();
    int tmpi=int(box.x*0.5+0.5);
    cutoff=tmpi;
    maxbin=tmpi*10;
    //SegRgDis test: cout << ""  
}

void SegRgDis::compute()
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


void export_SegRgDis()
{
    class_<SegRgDis, bases<AnalFunc> , boost::noncopyable>
        ("SegRgDis", init<boost::shared_ptr<Traj>, boost::shared_ptr<BasicInfo>, unsigned >() )
        ;
}

