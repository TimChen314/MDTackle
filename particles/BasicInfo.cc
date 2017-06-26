#include "BasicInfo.h"
using namespace boost::python;

void BasicInfo::MakeInfo() // Make MolNA, MolType, AtomSubj, a_bidx, a_eidx, m_bidx, m_eidx
{
    unsigned tmp_sum=0;
    unsigned tmp_idx=0;
    for(unsigned i=0;i<nm.size();++i) {
        for(unsigned int j=tmp_sum;j< tmp_sum + nm[i] ;++j ) {
            MolNA.push_back(na[i]);
            MolType.push_back(i);
            mol_a_bidx.push_back(tmp_idx);
            tmp_idx += na[i];
            mol_a_eidx.push_back(tmp_idx-1);
            for(unsigned int k=0;k<na[i];k++) {
                AtomSubj.push_back(j);
            }
        }
        tmp_sum += nm[i];
    }

    tmp_idx=0;
    for(unsigned i=0;i<nm.size();++i) {
        a_bidx.push_back(tmp_idx);
        tmp_idx += na[i]*nm[i];
        a_eidx.push_back(tmp_idx-1);
    }
    tmp_idx=0;
    for(unsigned i=0;i<nm.size();++i) {
        m_bidx.push_back(tmp_idx);
        tmp_idx += nm[i];
        m_eidx.push_back(tmp_idx-1);
    }
}


BasicInfo::BasicInfo()
{
    unsigned nm_nano,length_nano,nkind;
    cout<<"######################################################################"<<endl;
    cout<<"Enter the number of Nano particles:"<<endl;
    cin>>nm_nano;
    cout<<"Enter the atom's number in a nanoparticle: "<<endl;
    cin>>length_nano;
    na_tot = length_nano*nm_nano;
    nm_tot = nm_nano;
    na_nano = nm_nano*length_nano;
    nm.push_back(nm_nano);
    na.push_back(length_nano);
    cout<<"Enter the number of kind of polymers: "<<endl;
    cin>>nkind;
    for(unsigned i=0;i<nkind;i++)
    {
        unsigned int nm_chain,length;
        cout<<"Enter the number of "<<i+1<<"'s polymer chains:"<<endl;
        cin>>nm_chain;
        cout<<"Enter the length of "<<i+1<<"'s polymer chain:"<<endl;
        cin>>length;
        nm.push_back(nm_chain);
        na.push_back(length);
        na_tot += (length*nm_chain);
        nm_tot += nm_chain;
    }
    cout<<"######################################################################"<<endl;

    setInitPara();
    MakeInfo();// Make MolNA, MolType, AtomSubj, a_bidx, a_eidx, m_bidx, m_eidx

}




void BasicInfo::setInitPara()
{
    init_para.push_back(nm[0]);
    init_para.push_back(na[0]);
    init_para.push_back(nm.size() - 1);
    for(unsigned i=0; i!=nm.size() - 1;  ++ i )
    {
        init_para.push_back(nm[i]);
        init_para.push_back(na[i]);
    }

    ntype=nm.size();
    if(ntype==1)
        cerr<<"***** Warning ! Only NP in this system. *****"<<endl<<endl;
    else if(ntype==0)
    {
        cerr << "***** Error! No NP in this system *****" << endl;  
        throw runtime_error("Error! No NP in this system.");
    }
    else if(ntype > 100)
    {
        cerr<<"***** Warning ! Only NP in this system. *****"<<endl<<endl;
    }
}


void export_BasicInfo()
{
    class_<BasicInfo, boost::shared_ptr<BasicInfo> >
        ("BasicInfo",init<>() )
        ;
}


