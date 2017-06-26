#include <string>
#include <vector>
#include <tuple>

#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>

#include <math.h>
#include <algorithm>

#include <boost/python.hpp>
using namespace std;


// There are three level in BasicInfo: 1.molecular type 2. moleculer info 3.atom info
/* From "molecular type" we can get "moleculer info"(e.g. na, nm, m_bidx, m_eidx), 
 * and backward from "moleculer info" we can get "molecular type" (e.g. MolType) ;
 * from "moleculer info" we can get "atom info" (e.g. MolNA), and backward from "atom info" we can get "moleculer info" (e.g. AtomSubj) ;
 * from "molecular type" we can get "atom info" (e.g. a_bidx, a_eidx).
 */

#ifndef __BASICINFO_H__
#define __BASICINFO_H__

class BasicInfo
{
    public:
        BasicInfo();
        ~BasicInfo() {};

        unsigned na_tot;
        unsigned nm_tot;
        unsigned na_nano;
        unsigned ntype;
        vector<unsigned> na;
        vector<unsigned> init_para;
        vector<unsigned> nm;
        vector<unsigned> a_bidx; // a_bidx[i] storage the beginning atom index of i'th kind of molecule.
        vector<unsigned> a_eidx; // a_eidx[i] storage the ending atom index of i'th kind of molecule.
        vector<unsigned> m_bidx; // m_bidx[i] storage the beginning molecular index of i'th kind of molecule.
        vector<unsigned> m_eidx; // m_eidx[i] storage the ending molecular index of i'th kind of molecule.
        vector<unsigned> b_bidx; // b_bidx[i] storage the beginning bond index of i'th kind of molecule. 
        vector<unsigned> b_eidx; 
        vector<unsigned> agl_bidx; // agl_bidx[i] storage the beginning agl index of i'th kind of molecule.
        vector<unsigned> agl_eidx;
        vector<unsigned int> AtomSubj; //get AtomSubjection, so we can know a atom belong to which molecule.
        vector<unsigned int> MolNA; //get MolNA, so we can know how many atoms in this molecule
        vector<unsigned int> MolType; //get MolType, so we can know which type this molecule is.
        vector<unsigned int> mol_a_bidx;
        vector<unsigned int> mol_a_eidx;

    private:
        void setInitPara();
        void MakeInfo();// Make MolNA, MolType, AtomSubj, a_bidx, a_eidx, m_bidx, m_eidx

};

void export_BasicInfo();

#endif


