#include "MolDataStruct.h"
#include "DumpXml.h"


DumpXml::DumpXml(const string& fname, set<string>& ExistTerm, vec &box, vector<vec> &pos, vector<vec_int> &img, vector<vec> &vel, vector<unsigned> &type, vector<double> &mass, vector<Bond> &bond, vector<Angle> &angle, vector < string > & typemapping)
{
    for(set<string>::iterator iter=ExistTerm.begin();iter!=ExistTerm.end(); ++iter) {
        vector<string>::iterator iter_v=find(AvailableTerm.begin(), AvailableTerm.end(), (*iter) );
        if(iter_v==AvailableTerm.end() ) {
            cerr<<"Error! DumpXml can't dump "<<(*iter)<<" info!"<<endl;
            throw runtime_error("Error! DumpXml can't dump info.");
        }
    }
    vector<string> str_type;
    for(auto v : type)
        str_type.push_back(typemapping[v]);
    write(fname, box, pos, img, vel, str_type, mass, bond, angle);
}

//  void DumpXml::DumpXml_new(const string& fname, Reader& build)



DumpXml::DumpXml(const string& fname, Reader& build)
{
    set<string> ExistTerm(build.getExistTerm());
    cout << "Info DumpXml: " ;
    for(set<string>::iterator iter=ExistTerm.begin();iter!=ExistTerm.end(); ++iter) {
        cout <<  *iter << " ";    
        vector<string>::iterator iter_v=find(AvailableTerm.begin(), AvailableTerm.end(), (*iter) );
        if(iter_v==AvailableTerm.end()) {
            cerr<<"Error! DumpXml can't dump "<<(*iter)<<" info!"<<endl;
            throw runtime_error("Error! DumpXml can't dump info.");
        }
    }
    cout << endl; 

    cout<<"Info: Begin to dump "<<fname<<" !"<<endl;
    ofstream f(fname.c_str());
    if (!f.good())
    {
        cerr << endl << "***Error! Unable to open dump file for writing " << endl << endl;
        throw runtime_error("Error writting Xml dump file");
    }

    // Begin to Dump 
    DumpHead(f, build.getBoxPos());
    DumpBox(f, build.getBox());
    auto it = ExistTerm.begin();
    // if certain term exist, then dump it. 
    if( (it=find(ExistTerm.begin(), ExistTerm.end(), "position")) != ExistTerm.end()) DumpNode < vec > (f, *it, build.getBoxPos());
    if( (it=find(ExistTerm.begin(), ExistTerm.end(), "image")) != ExistTerm.end()) DumpNode < vec_int > (f, *it, build.getImage());
    if( (it=find(ExistTerm.begin(), ExistTerm.end(), "velocity")) != ExistTerm.end()) DumpNode < vec > (f, *it, build.getVel());
    if( (it=find(ExistTerm.begin(), ExistTerm.end(), "type")) != ExistTerm.end()) DumpNode < string > (f, *it, build.getStringType());
    if( (it=find(ExistTerm.begin(), ExistTerm.end(), "mass")) != ExistTerm.end()) DumpNode < double > (f, *it, build.getMass());
    if( (it=find(ExistTerm.begin(), ExistTerm.end(), "bond")) != ExistTerm.end()) DumpNode < Bond > (f, *it, build.getBond());
    if( (it=find(ExistTerm.begin(), ExistTerm.end(), "angle")) != ExistTerm.end()) DumpNode < Angle > (f, *it, build.getAngle());
    DumpTail(f);

}



void DumpXml::DumpHead(ofstream& f, vector < vec > & pos)
{
    f << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << "\n";
    f << "<polymer_xml version=\"1.3\">" << "\n";
    f << "<configuration time_step=\"" << 0 << "\" "
        << "dimensions=\"" << "3" << "\" "
        << "natoms=\"" << pos.size()  << "\" "
        << ">" << "\n";
}
void DumpXml::DumpBox(ofstream & f, const vec& box)
{
    f << "<box lx=\""   << box.x << "\" ly=\""  << box.y << "\" lz=\""  <<  box.z << "\"/>\n"; 
}


void DumpXml::DumpTail(ofstream & f)
{
    f << "</configuration>" << "\n";
    f << "</polymer_xml>" << "\n";
} 




void DumpXml::write(const string& fname, vec &box, vector<vec> &pos, vector<vec_int> &img, vector<vec> &vel, vector<string> &type, vector<double> &mass, vector<Bond> &bond, vector<Angle> &angle)
{
    cout<<"Info: Begin to dump "<<fname<<" !"<<endl;
    ofstream f(fname.c_str());

    if (!f.good())
    {
        cerr << endl << "***Error! Unable to open dump file for writing " << endl << endl;
        throw runtime_error("Error writting Xml dump file");
    }


    f << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << "\n";
    f << "<polymer_xml version=\"1.3\">" << "\n";
    f << "<configuration time_step=\"" << 0 << "\" "
        << "dimensions=\"" << "3" << "\" "
        << "natoms=\"" << pos.size()  << "\" "
        << ">" << "\n";
    f << "<box " << "lx=\""<< box.x << "\" ly=\""<< box.y << "\" lz=\""<< box.z << "\"/>" << "\n";

    f << "<position num=\"" << pos.size() << "\">" << "\n";
    for(unsigned int i=0;i!=pos.size();i++)
        f<<pos[i].x<<" "<<pos[i].y<<" "<<pos[i].z<<"\n";
    f <<"</position>" << "\n";
    f << "<image num=\"" << img.size() << "\">" << "\n";
    for(unsigned int i=0;i!=img.size();i++)
        f<<img[i].x<<" "<<img[i].y<<" "<<img[i].z<<"\n";
    f <<"</image>" << "\n";
    if(vel.size()!=0) {
        f << "<velocity num=\"" << vel.size() << "\">" << "\n";
        for(unsigned int i=0;i!=vel.size();i++)
            f<<vel[i].x<<" "<<vel[i].y<<" "<<vel[i].z<<"\n";
        f <<"</velocity>" << "\n";
    }
    f <<"<type num=\"" << type.size() << "\">" << "\n";
    for(unsigned int i=0;i!=type.size();i++)
        f<<type[i]<<"\n";
    f <<"</type>" << "\n";
    f <<"<mass num=\"" << mass.size() << "\">" << "\n";
    for(unsigned int i=0;i!=mass.size();i++)
        f<<mass[i]<<"\n";
    f <<"</mass>" << "\n";
    f << "<bond num=\"" << bond.size() << "\">" << "\n";
    for(unsigned int i=0;i!=bond.size();i++)
        f<<bond[i].type<<" "<<bond[i].a<<" "<<bond[i].b<<"\n";
    f << "</bond>" << "\n";
    if(angle.size()!=0) {
        f << "<angle num=\"" << angle.size() << "\">" << "\n";
        for(unsigned int i=0;i!=angle.size();i++)
            f<<angle[i].type<<" "<<angle[i].a<<" "<<angle[i].b<<" "<<angle[i].c<<"\n";
        f << "</angle>" << "\n";
    }
    f << "</configuration>" << "\n";
    f << "</polymer_xml>" << "\n";

    if (!f.good())
    {
        cerr << endl << "***Error! Unexpected error writing Xml dump file" << endl << endl;
        throw runtime_error("Error writting crosslink_bondinfo dump file");
    }

}



