#include "MolDataStruct.h"
#include "Reader.h"

#include <sys/time.h>
#include <string>
#include <vector>
#include <map>
#include <set>

#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>

#include <math.h>
#include <algorithm>

using namespace std;

#ifndef __DUMPXML_H__
#define __DUMPXML_H__


class DumpXml
{
    public:
        DumpXml(const string& fname, set<string>& ExistTerm, vec &box, vector<vec> &pos, vector<vec_int> &img, vector<vec> &vel, vector<unsigned> &type, vector<double> &mass, vector<Bond> &bond, vector<Angle> &angle, vector < string > & typemapping);
        DumpXml(const string& fname, Reader& build);  

    private:
        vector<string> AvailableTerm = {"position", "box",  "velocity",  "image",  "type",  "mass",  "bond",  "angle"}; // the info node supported by this program.
        void write(const string& fname, vec &box, vector<vec> &pos, vector<vec_int> &img, vector<vec> &vel, vector<string> &type, vector<double> &mass, vector<Bond> &bond, vector<Angle> &angle);
        void DumpHead(ofstream& f, vector < vec > & pos);
        void DumpBox(ofstream& f, const vec & box);
        void DumpTail(ofstream& f);
        template  < class T > void DumpNode(ofstream& f, const string & name, const vector < T > & vec_data);

};

    template  < class T > 
void DumpXml::DumpNode(ofstream & f, const string & name, const vector < T > & vec_data)
{
    f << "<"  << name << " num=\"" << vec_data.size() << "\">" << "\n";
    for(unsigned i=0;i!=vec_data.size();  ++ i )
    {
        f << vec_data[i] << "\n" ; 
    }
    f << "</"  << name << ">" << "\n";
}

#endif
