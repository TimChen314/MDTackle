#include "MolDataStruct.h"
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


#ifndef __MOLINFO_H__
#define __MOLINFO_H__

class MolInfo
{
    public:
        MolInfo(set<unsigned> &OMI, vector<vec> &pos, vector<vec_int> &img, vector<vec> &vel, vector<double> &mass, vector<unsigned> &type, \
                set<Bond> &OMB, set<Angle> &OMA); // OMI is the acronym of OneMolIdx
        MolInfo(set<unsigned> &OMI, vector<vec> &pos, vector<vec_int> &img, vector<vec> &vel, vector<double> &mass, vector<unsigned> &type, \
                set<Bond> &OMB);


        set<unsigned> OneMolIdx;  
        vector<Bond> m_bond;
        vector<Angle> m_angle;
        // storage the bead index of a mol, NOTE that these beads must be     
        // arranged continuously(i.e. the bead index range of first mol should be like 0~250).
        vector<vec> m_pos;
        vector<vec_int> m_img;
        vector<vec> m_vel;
        vector<double> m_mass;
        vector<unsigned> m_type;

    private:
        void ReMapIdx(); 
        // Remap index in mol[i].m_bond and mol[i].m_angle, even if m_bond.size() = 0 or 
        // m_angle.size() = 0, this subroutine will be ok.
        friend ostream& operator<<(ostream& os, const MolInfo &oneMol);

};





#endif



