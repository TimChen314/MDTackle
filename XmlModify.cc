#include "XmlModify.h"
#include <numeric>  // for generic algorithm accumulate
using namespace boost::python;
using namespace boost;


XmlModify::XmlModify(boost::shared_ptr<Reader> read)
{
    build = read;
    box=build->getBox();
    pos=build->getBoxPos();
    vel=build->getVel();
    image=build->getImage();
    type=build->getIdxType();
    typemapping = build -> getTypeMapping();
    mass=build->getMass();
    bond=build->getBond();
    angle=build->getAngle();
    ExistTerm=build->getExistTerm();
}


void XmlModify::XmlCutTail()
{
    string fname("XmlCutTail.xml");

    unsigned length,nm_be_cut;
    cout<<"Enter the length of polymer: "<<endl;	
    cin>>length;
    cout<<"Enter the number of polymer that want to be cut: "<<endl;
    cin>>nm_be_cut;

    unsigned na_be_cut=nm_be_cut*length;
    unsigned int N = build->getNumParticles();
    unsigned na = N - na_be_cut;

    DumpXml dump(fname, ExistTerm, box, VecResize(pos,na), VecResize(image,na), VecResize(vel,na), VecResize(type,na), VecResize(mass,na), \
            VecResize(bond, bond.size()-nm_be_cut*(length-1)), VecResize(angle,angle.size()-nm_be_cut*(length-2)), typemapping);
} 




void XmlModify::XmlCutOutOneSphere(unsigned na_sphere)
{
    const vec SphereCenter(0.0,0.0,0.0); // NOTE!!! : if SphereCenter is not set to (0.0,0.0,0.0), then pbc should be considered.
    double radius=-1; // if there is something wrong, the result will be obvious.

    unsigned nshell=int(box.x*0.5); // note: shell thickness is 1.0
    vector<unsigned> na_shell(nshell,0);
    for(unsigned i=0;i!=pos.size();++i) {
        vec d = pos[i] - SphereCenter; // operator??
        double dr = sqrt( ScalarProduct(d,d) );
        unsigned bin = int(dr);  // note: shell thickness is 1.0
        if(bin<nshell) 
            ++na_shell[bin];
    }

    unsigned sum=0;
    double upshell=-1,lowshell=-1; // if there is something wrong, the result will be obvious.
    for(unsigned i=0;i!=nshell;++i) {
        sum += na_shell[i]; // note: sum is the num of atom wtihin shpere
        if(na_sphere<=sum) { // note: radius is between lowshell and upshell. The value of lowshell also is the index of corresponding shell.
            lowshell = i;
            upshell = i + 1;
            break;
        }
        if(i==nshell-1) {
            cerr<<"*****Too many beads in sphere. *****"<<endl;
            throw runtime_error("Too many beads in sphere.");
        }
    }
    // test: for(unsigned i=0;i!=nshell;++i) cout<<na_shell[i]<<" na_shell "<<endl;

    unsigned nsubshell = na_shell[lowshell]*10;
    double dsubshell = 1.0/nsubshell;  
    // note: dsubshell is the shell thickness of subshell,  
    // and nsubshell is defined as na_shell[lowshell]*10 to ensure dsubshell is small enough(generally, a 
    // relatively large nsubshell will be good.).
    vector<unsigned> na_subshell(nsubshell,0);
    for(unsigned i=0;i!=pos.size();++i) {
        vec d = pos[i] - SphereCenter; 
        double dr = sqrt( ScalarProduct(d,d) ) -lowshell;
        if(dr>0.0 && dr<1.0 ) {
            unsigned bin = int(dr/dsubshell);  
            ++na_subshell[bin];
        }
    }
    unsigned tmp_na_sphere=accumulate(na_shell.begin(), na_shell.begin()+lowshell, 0);
    // note: for now, tmp_na_sphere must be < na_sphere.
    // test: for(unsigned i=0;i!=nsubshell;++i) cout<<na_subshell[i]<<" "<<endl;
    // test: cout<<tmp_na_sphere<<" "<<tmp_na_sphere+na_shell[lowshell]<<" "<<accumulate(na_subshell.begin(), na_subshell.end(), 0)<<endl;
    for(unsigned i=0;i!=nsubshell;++i) {
        tmp_na_sphere += na_subshell[i];
        if(tmp_na_sphere == na_sphere) {
            radius = lowshell+dsubshell*(i+1);
            // test: cout<<lowshell<<" "<<dsubshell*(i+1)<<endl;
            break;
        }
        if(i==nsubshell-1) {
            cerr<<"*****Error ! Can't get the value of radius correctly. *****"<<endl;
            throw runtime_error("Error ! Can't get the value of radius correctly.");
        }
    }
    cout<<"Results: Radius of sphere is "<<radius<<endl;


    // check if necessary terms exists.
    set<string>::iterator it1 = ExistTerm.find("position");
    set<string>::iterator it2 = ExistTerm.find("type");
    set<string>::iterator it3 = ExistTerm.find("bond");
    if( !(it1!=ExistTerm.end()&&it2!=ExistTerm.end()&&it3!=ExistTerm.end()) ) {
        cerr<<"*****Error! necessary term has not been found. *****"<<endl;
        throw runtime_error("Error! necessary term has not been found.");
    }
    // tackle the info
    vector<vec> spos; // pos of sphere
    vector<unsigned> sidx; // index of bead in sphere 
    for(unsigned i=0;i!=pos.size();++i) {
        vec d = pos[i] - SphereCenter; // note: Use SphereCenter in here
        double dr = sqrt( ScalarProduct(d,d) );
        if(dr<radius) {
            spos.push_back(pos[i]);
            sidx.push_back(i); // note: sidx is not ranged in sequence.
        }
    }
    if(sidx.size()!=na_sphere) {
        cerr<<"*****Error ! Number of beads in sphere is wrong!. *****"<<endl;
        throw runtime_error("Error ! Number of beads in sphere is wrong!.");
    }

    /* get the unessential term:
       set<string> UselessTerm(ExistTerm);
       UselessTerm.erase("position");
       UselessTerm.erase("image");
       UselessTerm.erase("bond");
       UselessTerm.erase("angle"); // Generally, we will use this sphere to crosslinking a ball, so \
       we only need the bond
       */


    vector<vec_int> simage;
    vector<double> smass;
    vector<unsigned int> stype;
    for(unsigned i=0;i!=sidx.size();++i) // Whether the image exist doesn't matter
        simage.push_back(vec_int(0,0,0));
    if(ExistTerm.find("type")!=ExistTerm.end()) {
        for(unsigned i=0;i!=sidx.size();++i)
            stype.push_back(type[ sidx[i] ]);
    }
    if(ExistTerm.find("mass")!=ExistTerm.end()) {
        for(unsigned i=0;i!=sidx.size();++i)
            smass.push_back(mass[ sidx[i] ]);
    }


    map<unsigned,unsigned> idx2sidx;
    for(unsigned i=0;i!=sidx.size();++i) {
        idx2sidx.insert( make_pair(sidx[i],i) ); // idx2sidx[idx]=sidx; while for sidx, sidx[sidx]=idx;
    }
    set<Bond> sbond; // note: we want to preserve the old bonds, and in next step we will create new bonds.
    set<unsigned> set_sidx(sidx.begin(),sidx.end()); 
    // note: Below we'll use "find" algorithm, so use of 
    // "set" will be more efficiency.
    for(vector<Bond>::iterator bondit = bond.begin(); bondit !=bond.end(); ++bondit) {
        if( set_sidx.find(bondit->a)!=set_sidx.end()&&set_sidx.find(bondit->b)!=set_sidx.end() ) {
            sbond.insert(Bond(bondit->type,idx2sidx[bondit->a],idx2sidx[bondit->b] ));
        }
    }
    cout<<"CutOutOneSphere: Sphere is dumped to CutOutOneSphere_uncrosslinked.xml "<<sbond.size()<<endl;
    vector<vec> tmp_vel; // don't need vel and angle info, just use them to hold the agrument place.
    vector<Angle> tmp_angle;
    vector<Bond> v_sbond(sbond.begin(),sbond.end());
    cout<<"XmlModify::XmlCutOutOneSphere stype.size()"<<stype.size()<<endl;
    ReaderData rawdata(box,spos,simage,tmp_vel,stype,smass,v_sbond,tmp_angle, typemapping);

    string fname_raw("CutOutOneSphere.xml");
    DumpXml rawdump(fname_raw, rawdata);


    double threshold=0.55; // when the distance between two monomers < threshold, these two monomers will bonding.
    double binsize=0.01; // threshold will increase a "binsize" in one time.
    for(unsigned i=0;i!=50;++i) {
        for(unsigned j=0;j!=spos.size()-1;++j) {
            for(unsigned k=j+1;k!=spos.size();++k) {
                vec tmpvec(spos[j]-spos[k]);
                if(sqrt(ScalarProduct(tmpvec,tmpvec)) > threshold) continue;
                ostringstream os;
                os<<stype[j]<<stype[k];
                string bondtype(os.str());
                if(bondtype=="22") 
                    bondtype="11";
                else if(bondtype=="21")
                    bondtype="12";
                //sbond.insert(Bond(bondtype,j,k));
                sbond.insert(Bond("33",j,k));
            }
        }

        v_sbond.assign(sbond.begin(),sbond.end());
        ReaderData data(box,spos,simage,tmp_vel,stype,smass,v_sbond,tmp_angle, typemapping);
        TopoInfo topo(&data);
        topo.DumpMol();
        if(topo.getMolNum()==0) {
            cerr<<"*****Can't build the topology of OneSphere molecule. *****"<<endl;
            throw runtime_error("Can't build the topology of OneSphere molecule.");
        }
        else if(topo.getMolNum()==1) 
            break;
        threshold += binsize;
    }
    cout<<"XmlCutOutOneSphere: threshold is "<<threshold<<endl;

    v_sbond.assign(sbond.begin(),sbond.end());
    ReaderData data(box,spos,simage,tmp_vel,stype,smass,v_sbond,tmp_angle, typemapping);

    string fname("CutOutOneSphere_crosslinked.xml");
    DumpXml dump(fname, data); 
}

void export_XmlModify()
{
    class_<XmlModify, boost::shared_ptr<XmlModify> >
        ("XmlModify", init<boost::shared_ptr<Reader> >())
        .def("XmlCutTail", &XmlModify::XmlCutTail)
        .def("XmlMagnifyBox", &XmlModify::XmlMagnifyBox)
        .def("XmlCutOutOneSphere", &XmlModify::XmlCutOutOneSphere)
        ; 
}

