#include "BasicInfo.h"
#include "MolInfo.h"
#include "XmlModify.h"
#include "TopoInfo.h"

#include <sys/time.h>
#include <string>
#include <vector>
#include <set>

#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>

#include <math.h>
#include <algorithm>

using namespace std;


void MagnifyBox(vec_int magnify_vec, boost::shared_ptr<Reader> build, vector<vec> &a_pos, vector<vec_int> &a_img, vector<vec> &a_vel, vector<unsigned int> &a_type, vector<double> &a_mass, vector<Bond> &a_bond, vector<Angle> &a_angle, BasicInfo &bi)
{
    cout<<"Info: Begin to magnify box."<<endl;
    vec box=build->getBox();
    vector<vec> r_pos=build->getRealPos();
    vector<vec> vel=build->getVel();
    vector<unsigned int> type = build->getIdxType();
    vector<double> mass = build->getMass();
    vector<Bond> bond=build->getBond();
    vector<Angle> angle=build->getAngle();

    double Lx=box.x;
    double Ly=box.y;
    double Lz=box.z;
    double a_Lx=magnify_vec.x*Lx;
    double a_Ly=magnify_vec.y*Ly;
    double a_Lz=magnify_vec.z*Lz;
    double a_LxINV=1.0/a_Lx;
    double a_LyINV=1.0/a_Ly;
    double a_LzINV=1.0/a_Lz;

    vector<vector<vec> > allbox_r_pos;
    vector<vec> tmp_r_pos;
    for(int i=0;i!=magnify_vec.z;++i) {
        for(int j=0;j!=magnify_vec.y;++j) {
            for(int k=0;k!=magnify_vec.x;++k) { // loop along x axis first.
                vector<vec> tmp_r_pos(r_pos);
                for(unsigned l=0;l!=tmp_r_pos.size();++l) {
                    tmp_r_pos[l].x += k*Lx;
                    tmp_r_pos[l].y += j*Ly;
                    tmp_r_pos[l].z += i*Lz;
                }
                allbox_r_pos.push_back(tmp_r_pos);
                //box_r_pos[i+j*vec_int.x+k*vec_int.x*vec_int.y]
            }
        }
    }
    //cout<<bi.a_bidx.size()<<" "<<bi.a_eidx.size()<<endl;
    //cout<<bi.b_bidx.size()<<" "<<bi.b_eidx.size()<<endl;
    //cout<<bi.agl_bidx.size()<<" "<<bi.agl_eidx.size()<<endl;
    for(unsigned i=0;i!=bi.ntype;++i) {
        unsigned n_bond = bi.b_eidx[i] - bi.b_bidx[i] + 1;
        unsigned n_angle = bi.agl_eidx[i] - bi.agl_bidx[i] + 1;
        for(unsigned j=0;j!=allbox_r_pos.size();++j) {
            unsigned pos_size = a_pos.size();
            if(i>0)
                pos_size -= (bi.a_eidx[i-1] + 1);
            a_pos.insert(a_pos.end(), allbox_r_pos[j].begin()+bi.a_bidx[i], allbox_r_pos[j].begin()+bi.a_eidx[i]+1);
            if(vel.size()!=0) a_vel.insert(a_vel.end(), vel.begin()+bi.a_bidx[i], vel.begin()+bi.a_eidx[i]+1);
            a_type.insert(a_type.end(), type.begin()+bi.a_bidx[i], type.begin()+bi.a_eidx[i]+1);
            a_mass.insert(a_mass.end(), mass.begin()+bi.a_bidx[i], mass.begin()+bi.a_eidx[i]+1);
            a_bond.insert(a_bond.end(), bond.begin()+bi.b_bidx[i], bond.begin()+bi.b_eidx[i]+1);
            for(unsigned k=a_bond.size()-n_bond;k!=a_bond.size();++k) {
                a_bond[k].a += pos_size;
                a_bond[k].b += pos_size;
            }
            a_angle.insert(a_angle.end(), angle.begin()+bi.agl_bidx[i], angle.begin()+bi.agl_eidx[i]+1);
            for(unsigned k=a_angle.size()-n_angle;k!=a_angle.size();++k) {
                a_angle[k].a += pos_size;
                a_angle[k].b += pos_size;
                a_angle[k].c += pos_size;
            }
        }
    }
    for(unsigned i=0;i!=a_pos.size();++i) { // move center of mass back.
        a_pos[i].x -= (0.5*magnify_vec.x -0.5)*Lx;
        a_pos[i].y -= (0.5*magnify_vec.y -0.5)*Ly;
        a_pos[i].z -= (0.5*magnify_vec.z -0.5)*Lz;
    }
    for(unsigned i=0;i!=a_pos.size();++i) { // Make image and do PBC.
        vec_int tmp_veci(rint(a_pos[i].x*a_LxINV),rint(a_pos[i].y*a_LyINV),rint(a_pos[i].z*a_LzINV));
        a_img.push_back(tmp_veci);
        a_pos[i].x -= a_Lx*tmp_veci.x;
        a_pos[i].y -= a_Ly*tmp_veci.y;
        a_pos[i].z -= a_Lz*tmp_veci.z;
    }
}




void XmlModify::XmlMagnifyBox( BasicInfo& basicinfo, unsigned xtime, unsigned ytime,unsigned ztime)
{
    //unsigned del_idx,del_num;
    vec_int magnify_vec;
    cout<<"**************************************************************************************************************************"<<endl;
    cout<<"**   Note: New boxs place in positive direction along the axis, and all \"image\" info has reset, and same mol ranged   **"<<endl;
    cout<<"**   together in both old mol and new mol. (e.g. Assumption there are A mol and Bmol, then AABB magnify 2 times will    **"<<endl;
    cout<<"**   become AAAABBBB).                                                                                                  **"<<endl;
    cout<<"**                                                                                                                      **"<<endl;
    cout<<"**************************************************************************************************************************"<<endl;
    magnify_vec.x = xtime;
    magnify_vec.y = ytime;
    magnify_vec.z = ztime;
    if(magnify_vec.x<=0 || magnify_vec.y<=0 || magnify_vec.z<=0) {
        cerr<<"***** Error! Wrong magnification times! *****"<<endl;
        throw runtime_error(" Error! Wrong magnification times!");
    }
    if(vel.size()==0) cout<<"Note: No velocity information in this coor.xml."<<endl;
    if(magnify_vec.x==1 && magnify_vec.y==1 && magnify_vec.z==1) cerr<<"***** Warning! Don't need magnify any more! *****"<<endl;
    //cout<<"Info: Read is over!"<<endl;

    TopoInfo topo(build,basicinfo); // NOTE!! : Actually, we only get basicinfo.b_bidx and basicinfo.agl_bidx in here.

    vector<vec> a_pos;// this pos storage all particles.
    vector<vec_int> a_img;
    vector<vec> a_vel;
    vector<unsigned> a_type;
    vector<double> a_mass;
    vector<Bond> a_bond;
    vector<Angle> a_angle;
    vec a_box(box.x*magnify_vec.x,box.y*magnify_vec.y,box.z*magnify_vec.z);

    MagnifyBox(magnify_vec, build, a_pos, a_img, a_vel, a_type, a_mass, a_bond, a_angle, basicinfo); // magnification based on basicinfo

    DumpXml dump("MagnifiedBox.xml", ExistTerm, a_box, a_pos, a_img, a_vel, a_type, a_mass, a_bond, a_angle, typemapping);


}

