#include "MolInfo.h"

ostream& operator<<(ostream& os, const MolInfo &oneMol)
{
    os<<"+-----------------------------------------------------------+"<<endl;
    os<<"The number of bead in this mol : "<<oneMol.OneMolIdx.size()<<endl;
    os<<"The number of bond : "<<oneMol.m_bond.size()<<endl;
    os<<"The number of angle : "<<oneMol.m_angle.size()<<endl;
    os<<"Elements in OneMolIdx: "<<endl;
    for(set<unsigned>::iterator iter=oneMol.OneMolIdx.begin();iter!=oneMol.OneMolIdx.end();++iter) {
        os<<*iter<<endl;
    }
    os<<"Bond info: "<<endl;
    for(vector<Bond>::const_iterator iter=oneMol.m_bond.begin(); iter!=oneMol.m_bond.end() ;++iter)
        os<<iter->type<<" "<<iter->a<<" "<<iter->b<<endl;
    os<<"+-----------------------------------------------------------+"<<endl;
    return os;
}

MolInfo::MolInfo(set<unsigned> &OMI,vector<vec> &pos, vector<vec_int> &img, vector<vec> &vel, vector<double> &mass, vector<unsigned> &type, \
        set<Bond> &OMB, set<Angle> &OMA)
{
    cout<<"test MolInfo 1"<<" "<<OMI.size()<<" "<<pos.size()<<" "<<img.size()<<" "<<mass.size()<<" "<<type.size()<<" "<<m_bond.size()<<endl;
    OneMolIdx.insert(OMI.begin(),OMI.end());
    for(set<unsigned>::iterator iter=OneMolIdx.begin();iter!=OneMolIdx.end();++iter) {
        unsigned idx=*iter;
        m_pos.push_back(pos[idx]);
        m_img.push_back(img[idx]);
        if(vel.size()!=0) m_vel.push_back(vel[idx]);
        m_mass.push_back(mass[idx]);
        m_type.push_back(type[idx]);
    }
    m_bond.insert(m_bond.begin(),OMB.begin(),OMB.end());
    m_angle.insert(m_angle.begin(),OMA.begin(),OMA.end());

    ReMapIdx();// Remap index in mol[i].m_bond and mol[i].m_angle

}



MolInfo::MolInfo(set<unsigned> &OMI,vector<vec> &pos, vector<vec_int> &img, vector<vec> &vel, vector<double> &mass, vector<unsigned> &type, \
        set<Bond> &OMB)
{
    cout<<"test MolInfo 0"<<" "<<OMI.size()<<" "<<pos.size()<<" "<<img.size()<<" "<<mass.size()<<" "<<type.size()<<" "<<m_bond.size()<<endl;
    OneMolIdx.insert(OMI.begin(),OMI.end());
    for(set<unsigned>::iterator iter=OneMolIdx.begin();iter!=OneMolIdx.end();++iter) {
        unsigned idx=*iter;
        m_pos.push_back(pos[idx]);
        m_img.push_back(img[idx]);
        if(vel.size()!=0) m_vel.push_back(vel[idx]);
        m_mass.push_back(mass[idx]);
        m_type.push_back(type[idx]);
    }
    m_bond.insert(m_bond.begin(),OMB.begin(),OMB.end());

    ReMapIdx();// Remap index in mol[i].m_bond and mol[i].m_angle

}


void MolInfo::ReMapIdx()
{
    map<int,int> MapIdx;// By this method, no matter how the index is arranged we can tackle it.
    set<unsigned>::iterator iter=OneMolIdx.begin();
    for(unsigned j=0;j!=m_pos.size();++j) {
        MapIdx.insert(make_pair(*iter,j));
        ++iter;
    }
    for(unsigned j=0;j!=m_bond.size();++j) { // map the index of bond to the index in the new molecule.
        map<int,int>::iterator iter2=MapIdx.find(m_bond[j].a);
        if(iter2!=MapIdx.end()) {
            m_bond[j].a = iter2->second;
        }
        else
            cerr<<"***** Error! An atom in a bond is not found! *****"<<endl;
        iter2=MapIdx.find(m_bond[j].b);
        if(iter2!=MapIdx.end()) {
            m_bond[j].b = iter2->second;
        }
        else
            cerr<<"***** Error! An atom in a bond is not found! *****"<<endl;
    }
    for(unsigned j=0;j!=m_angle.size();++j) { // map the index of angle to the index in the new molecule.
        map<int,int>::iterator iter2=MapIdx.find(m_angle[j].a);
        if(iter2!=MapIdx.end()) {
            m_angle[j].a = iter2->second;
        }
        else
            cerr<<"***** Error! An atom in a angle is not found! *****"<<endl;
        iter2=MapIdx.find(m_angle[j].b);
        if(iter2!=MapIdx.end()) {
            m_angle[j].b = iter2->second;
        }
        else
            cerr<<"***** Error! An atom in a angle is not found! *****"<<endl;
        iter2=MapIdx.find(m_angle[j].c);
        if(iter2!=MapIdx.end()) {
            m_angle[j].c = iter2->second;
        }
        else
            cerr<<"***** Error! An atom in a angle is not found! *****"<<endl;
    }
}



