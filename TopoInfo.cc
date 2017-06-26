#include "TopoInfo.h"

TopoInfo::TopoInfo(boost::shared_ptr<Reader> build,BasicInfo& basicinfo)
{
    vector<Bond> bond=build->getBond();
    if(bond.size()==0) {
        cerr<<"TopoInfo: the number of bonds can't be 0 ! "<<endl;
        throw runtime_error("the number of bonds can't be 0 ");
    }
    set<string> ExistTerm=build->getExistTerm();
    AvailableTerm.push_back("position");
    AvailableTerm.push_back("angle");
    AvailableTerm.push_back("bond");
    AvailableTerm.push_back("image");
    AvailableTerm.push_back("type");
    AvailableTerm.push_back("mass");
    // There is no need to check velocity.
    for(unsigned i=0;i!=AvailableTerm.size(); ++i) {
        set<string>::iterator iter=ExistTerm.find(AvailableTerm[i]);
        if(iter==ExistTerm.end()) {
            cerr<<"Error! TopoInfo: Can't get mol vector !"<<endl;
            throw runtime_error("Error! TopoInfo: Can't get mol vector !");
        }
    }
    vector<vec> pos=build->getBoxPos();
    vector<vec> vel=build->getVel();
    vector<vec_int> image=build->getImage();
    vector<unsigned> type=build->getIdxType();
    vector<double> mass=build->getMass();
    vector<Angle> angle=build->getAngle();

    getTopoFromBond(bond,pos.size());

    if(!IsTopoConsistency(basicinfo)) {
        cerr<<"***** Error! There are unconsistency between topology from basicinfo and from bond !!! *****"<<endl;
        throw runtime_error("Error! There are unconsistency between topology from basicinfo and from bond");
    }
    if(!IsAllParticleInclude(basicinfo.na_tot,pos.size())) {
        cerr<<"***** Error! Some particle are not included in topology!!! *****"<<endl;
        throw runtime_error(" Error! Some particle are not included in topology!");
    }
    if(!IsMolIdxContinuous()) {
        cerr<<"***** Error! MolIdx is not continuous!!! *****"<<endl;
        throw runtime_error(" Error! MolIdx is not continuous!");
    }

    getMolAngleFromMolIdx(angle);

    list<set<unsigned> >::iterator MolIdx_iter=MolIdx.begin();
    list<set<Bond> >::iterator MolBond_iter=MolBond.begin();
    list<set<Angle> >::iterator MolAngle_iter=MolAngle.begin();
    for(;MolIdx_iter!=MolIdx.end();++MolIdx_iter,++MolBond_iter,++MolAngle_iter) {
        MolInfo oneMol(*MolIdx_iter, pos, image, vel, mass, type, *MolBond_iter, *MolAngle_iter);
        Mol.push_back(oneMol);
    }
    cout<<"TopoInfo: getMol is over!"<<endl;

    AddBondIdxToBI(basicinfo);
}




TopoInfo::TopoInfo(Reader * build)
{
    vector<Bond> bond=build->getBond();
    if(bond.size()==0) {
        cerr<<"TopoInfo: the number of bonds can't be 0 ! "<<endl;
        throw runtime_error("the number of bonds can't be 0 ");
    }
    set<string> ExistTerm=build->getExistTerm();
    AvailableTerm.push_back("position");
    AvailableTerm.push_back("bond");
    AvailableTerm.push_back("image");
    AvailableTerm.push_back("type");
    AvailableTerm.push_back("mass");
    // There is no need to check velocity.
    for(unsigned i=0;i!=AvailableTerm.size(); ++i) {
        set<string>::iterator iter=ExistTerm.find(AvailableTerm[i]);
        if(iter==ExistTerm.end()) {
            cerr<<"Error! TopoInfo: Can't get mol vector !"<<endl;
            throw runtime_error("Error! TopoInfo: Can't get mol vector !");
        }
    }
    vector<vec> pos=build->getBoxPos();
    vector<vec> vel=build->getVel();
    vector<vec_int> image=build->getImage();
    vector<unsigned> type=build->getIdxType();
    vector<double> mass=build->getMass();
    vector<Angle> angle=build->getAngle();

    getTopoFromBond(bond,pos.size());
    if(angle.size()!=0) getMolAngleFromMolIdx(angle);
    list<set<unsigned> >::iterator MolIdx_iter=MolIdx.begin();
    list<set<Bond> >::iterator MolBond_iter=MolBond.begin();
    list<set<Angle> >::iterator MolAngle_iter=MolAngle.begin();
    for(;MolIdx_iter!=MolIdx.end();++MolIdx_iter,++MolBond_iter,++MolAngle_iter) {
        MolInfo * oneMol;
        if(angle.size()!=0) 
            oneMol = new MolInfo(*MolIdx_iter, pos, image, vel, mass, type, *MolBond_iter, *MolAngle_iter);
        else {
            oneMol = new MolInfo(*MolIdx_iter, pos, image, vel, mass, type, *MolBond_iter);
        }
        Mol.push_back(*oneMol);
    }
    cout<<"TopoInfo: getMol is over!"<<endl;
}





// get b_bidx, agl_bidx, etc.  ! Note: In here, there is an implicit assumption that the info of bond in the class "buid" is sequenced by molecule .
void TopoInfo::AddBondIdxToBI(BasicInfo& basicinfo)
{
    unsigned tmp_sum1=0;
    unsigned tmp_sum2=0;
    for(unsigned j=0;j!=basicinfo.ntype;++j) {
        unsigned tmp_i1=0;
        unsigned tmp_i2=0;
        for(unsigned k=basicinfo.m_bidx[j];k<=basicinfo.m_eidx[j];++k) {
            tmp_i1 += Mol[k].m_bond.size();
            tmp_i2 += Mol[k].m_angle.size();
        }
        basicinfo.b_bidx.push_back(tmp_sum1);
        tmp_sum1 += tmp_i1;
        basicinfo.b_eidx.push_back(tmp_sum1-1);
        basicinfo.agl_bidx.push_back(tmp_sum2);
        tmp_sum2 += tmp_i2;
        basicinfo.agl_eidx.push_back(tmp_sum2-1);
    }
}


bool TopoInfo::IsTopoConsistency(BasicInfo &bi )
{
    if(bi.MolNA.size()!=MolIdx.size())
        return false;
    list<set<unsigned> >::iterator MolIdx_iter=MolIdx.begin();
    vector<unsigned>::iterator MolNA_iter=bi.MolNA.begin();
    for(;MolIdx_iter!=MolIdx.end();++MolIdx_iter,++MolNA_iter) {
        if(*MolNA_iter!=MolIdx_iter->size())
            return false;
    }
    return true;
}


bool TopoInfo::IsAllParticleInclude(unsigned bi_na_tot,unsigned pos_size)
{
    if(bi_na_tot!=pos_size)
        return false;
    return true;
}


bool TopoInfo::IsMolIdxContinuous()
{
    vector<unsigned> fullIdx;
    list<set<unsigned> >::iterator MolIdx_iter=MolIdx.begin();
    for(;MolIdx_iter!=MolIdx.end();++MolIdx_iter) {
        fullIdx.insert(fullIdx.end(),MolIdx_iter->begin(),MolIdx_iter->end());
    }
    for(unsigned i=0;i!=fullIdx.size();++i) {
        if(fullIdx[i]!=i) return false;
    }
    return true;
}





void TopoInfo::getMolAngleFromMolIdx(vector<Angle> &m_angle )
{
    // test: cout<<"TopoInfo::getMolAngleFromMolIdx : Subroutine start."<<endl;
    list<Angle> tmp_angle(m_angle.begin(),m_angle.end());

    list<set<unsigned> >::iterator MolIdx_iter=MolIdx.begin();
    for(;MolIdx_iter!=MolIdx.end();++MolIdx_iter) {
        set<Angle> OneMolAngle;

        list<Angle>::iterator iter=tmp_angle.begin();

        while( iter!=tmp_angle.end() && (MolIdx_iter->find(iter->a)!=MolIdx_iter->end() && MolIdx_iter->find(iter->b)!=MolIdx_iter->end() && MolIdx_iter->find(iter->c)!=MolIdx_iter->end()) ) { 
            // if index a, b and c of present angle is found in present MolIdx, then this angle is inserted into present element of MolAngle
            OneMolAngle.insert(*iter);
            iter=tmp_angle.erase(iter);
        } // find all angles in this mol.
        MolAngle.push_back(OneMolAngle);
    }

    if(tmp_angle.size()!=0 ) {
        cerr<<"***** Error! An Angle doesn't belong to any Mol!!! *****"<<endl;
        throw runtime_error(" Error! An Angle doesn't belong to any Mol!");
    }
    if(MolAngle.size()!=MolIdx.size()) {
        cerr<<"***** Error! Number of Mol in MolAngle and MolIdx doesn't match!!! *****"<<endl;
        throw runtime_error(" Error! Number of Mol in MolAngle and MolIdx doesn't match!");
    }
    cout<<"TopoInfo: getMolAngleFromMolIdx is over!"<<endl;
}


void TopoInfo::getTopoFromBond(vector<Bond> &m_bond, unsigned nbeads)
{
    list<Bond> tmp_bond(m_bond.begin(),m_bond.end());
    list<Bond>::iterator iter=tmp_bond.begin(); 

    Bond bond_idx2("11",1,153);
    while(iter!=tmp_bond.end()) {  
        // Find the mol that *iter belong to. One loop for one bond.
        // For the rest, they are one-bead molecule, and will be added into MolIdx later.
        list<set<unsigned> >::reverse_iterator MolIdx_r_iter=MolIdx.rbegin();
        list<set<Bond> >::reverse_iterator MolBond_r_iter=MolBond.rbegin();
        unsigned flag=1;
        // To identify whether this bond is included in MolBond; 1 means included; 0 
        // unincluded and it will be included by add a new mol.
        // test: if(tmp_bond.size()%10==0) cout<<"MolIdx.size "<<tmp_bond.size()<<endl;
        for(;MolIdx_r_iter!=MolIdx.rend();++MolIdx_r_iter,++MolBond_r_iter) {
            if(MolIdx_r_iter->find(iter->a)==MolIdx_r_iter->end() && MolIdx_r_iter->find(iter->b)==MolIdx_r_iter->end() ) { 
                // if ture, present bond doesn't belong to this mol. So we "continue" and try next mol.
                continue;
            }
            else
                flag *= 0;
            if(MolIdx_r_iter->find(iter->a)!=MolIdx_r_iter->end() && MolIdx_r_iter->find(iter->b)!=MolIdx_r_iter->end() ) { 
                if(MolBond_r_iter->find(*iter)==MolBond_r_iter->end()) {
                    MolBond_r_iter->insert(*iter);
                    iter=tmp_bond.erase(iter);
                    break; // find the owner mol of present bond, so "break" and try next bond.
                }
                else {
                    cerr<<"***** Error! The same bond has been found twice!!!  *****"<<endl;
                    throw runtime_error(" Error! The same bond has been found twice!");
                }
            }
            if( MolIdx_r_iter->find(iter->a)!=MolIdx_r_iter->end() && MolIdx_r_iter->find(iter->b)==MolIdx_r_iter->end() ) { 
                // if true, we find a bond of present mol.
                MolIdx_r_iter->insert(iter->b);
                MolBond_r_iter->insert(*iter);
                mergeTwoMol(MolIdx, MolBond, MolIdx_r_iter, MolBond_r_iter, iter->b); 
                // if iter->b (a bead) is found in another mol, merge these two mol together, because they are the same mol.
                iter=tmp_bond.erase(iter);
                break;
            }
            if( MolIdx_r_iter->find(iter->a)==MolIdx_r_iter->end() && MolIdx_r_iter->find(iter->b)!=MolIdx_r_iter->end() ) {  
                // if true, we find a bond of present mol.
                MolIdx_r_iter->insert(iter->a);
                MolBond_r_iter->insert(*iter);
                mergeTwoMol(MolIdx, MolBond, MolIdx_r_iter, MolBond_r_iter, iter->a);
                iter=tmp_bond.erase(iter);
                break;
            }
        }
        if(flag==1) // present bond doesn't belong to any known Mol, nedd to create a new mol.
            MakeMultiBeadsMolAndInsert(iter,tmp_bond);
    }

    set<unsigned> fullIdx;
    unsigned sum_size=0;
    list<set<unsigned> >::iterator MolIdx_iter=MolIdx.begin();
    for(;MolIdx_iter!=MolIdx.end();++MolIdx_iter) {
        sum_size += MolIdx_iter->size();
        fullIdx.insert(MolIdx_iter->begin(),MolIdx_iter->end());
    }
    if(fullIdx.size()!=sum_size) {
        cerr<<"***** Warning! The bond part of input is not ranged in sequence with molecule. *****"<<endl;
    }
    set<unsigned> rawIdx;
    for(unsigned i=0;i!=nbeads;++i)
        rawIdx.insert(i);
    vector<unsigned> OneBeadMolIdx;
    set_difference(rawIdx.begin(), rawIdx.end(), fullIdx.begin(), fullIdx.end(), back_inserter(OneBeadMolIdx) );
    for(vector<unsigned>::iterator OBMI_iter=OneBeadMolIdx.begin();OBMI_iter!=OneBeadMolIdx.end();++OBMI_iter)
        MakeOneBeadMolAndInsert(*OBMI_iter);
    cout<<"TopoInfo: getTopoFromBond : The number of one-bead mols is "<<OneBeadMolIdx.size()<<endl;
    /* test:
       for(vector<unsigned>::iterator iter=OneBeadMolIdx.begin(); iter!=OneBeadMolIdx.end(); ++iter)
       cout<<"TopoInfo: getTopoFromBond : The index of one-bead mols: "<<*iter<<endl;
       */


    if(!CheckMolBondNum(m_bond.size())) {
        cerr<<"***** Error! Bond Number is wrong !!! *****"<<endl;
        throw runtime_error(" Error!  Bond Number is wrong!");
    }
    cout<<"TopoInfo: getTopoFromBond is over! Number of Mol is "<<MolIdx.size()<<". "<<endl;
}




void TopoInfo::MakeMultiBeadsMolAndInsert(list<Bond>::iterator &iter,list<Bond> &tmp_bond )
{
    set<unsigned> OneMolIdx;
    set<Bond> OneMolBond;
    OneMolIdx.insert(iter->a);
    OneMolIdx.insert(iter->b);
    OneMolBond.insert(*iter);
    iter=tmp_bond.erase(iter);// Note: There's no need to break because when this commond has down the next loop begins.
    MolIdx.push_back(OneMolIdx);
    MolBond.push_back(OneMolBond);
}




void TopoInfo::mergeTwoMol(list<set<unsigned> > & MolIdx, list<set<Bond> > & MolBond, list<set<unsigned> >::reverse_iterator & \
        MolIdx_r_iter, list<set<Bond> >::reverse_iterator & MolBond_r_iter, unsigned beadidx)
{
    list<set<unsigned> >::reverse_iterator MolIdx_search=MolIdx.rbegin();
    list<set<Bond> >::reverse_iterator MolBond_search=MolBond.rbegin();
    for(;MolIdx_search!=MolIdx.rend();++MolIdx_search,++MolBond_search) { 
        // beadidx may be belong to another mol, if that's true, these two mols is actually one mol.
        // So first we search if beadidx is belong to other mol. Note that the present mol shouldn't be searched. 
        // And note that reverse_iterator type varible and iterator type varible can't be compared,  
        // hence MolIdx_search must be defined as reverse_iterator type.
        if(MolIdx_search==MolIdx_r_iter) continue; 
        if(MolIdx_search->find(beadidx)!=MolIdx_search->end()) { 
            // merge these two mol together.
            // note : Since we'll break out of this closest loop, whether MolIdx_search is a dangling pointer doesn't matter.
            MolIdx_r_iter->insert(MolIdx_search->begin(),MolIdx_search->end());
            MolBond_r_iter->insert(MolBond_search->begin(),MolBond_search->end());
            MolIdx_r_iter=list<set<unsigned> >::reverse_iterator( MolIdx.erase(--MolIdx_search.base()) );
            // C++ : This part is a tricky, it's because erase only accept iterator type.
            // For more info, see http://blog.csdn.net/kesalin/article/details/24265303
            MolBond_r_iter=list<set<Bond> >::reverse_iterator(MolBond.erase(--MolBond_search.base()) );
            break;
        }
    }
}



bool TopoInfo::CheckMolBondNum(unsigned m_bond_size )
{
    unsigned n_bond=0;
    list<set<Bond> >::iterator MolBond_iter=MolBond.begin();
    for(;MolBond_iter!=MolBond.end();++MolBond_iter) {
        n_bond += MolBond_iter->size();
    }
    if(n_bond!=m_bond_size)
        return false;
    else
        return true;
}

void TopoInfo::DumpMol() const
{
    ofstream f("MolFromTopo.tackle_info");
    for(unsigned i=0;i!=Mol.size();++i) 
        f<<"### Mol "<<i+1<<"'th "<<endl<<Mol[i]<<endl;
}

