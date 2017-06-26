#include "BasicInfo.h"
#include "MolDataStruct.h"
#include "MolInfo.h"
#include "Reader.h"

#ifndef __TopoInfo_H__
#define __TopoInfo_H__

class TopoInfo 
{
    public:
        TopoInfo(boost::shared_ptr<Reader> builder, BasicInfo& basicinfo); // NOTE!! : Actually, we only get basicinfo.b_bidx and basicinfo.agl_bidx in here.
        TopoInfo(Reader * builder);
        unsigned getMolNum() const {return Mol.size();}
        void DumpMol() const ;

    private:
        void MakeMultiBeadsMolAndInsert(list<Bond>::iterator &iter,list<Bond> &tmp_bond );
        void MakeNewAngleAndInsert(list<Angle>::iterator &iter,list<Angle> &tmp_angle );
        inline void MakeOneBeadMolAndInsert(unsigned idx);
        void getTopoFromBond(vector<Bond> &m_bond, unsigned nbeads); // MolBond & MolIdx is established.
        void getMolAngleFromMolIdx(vector<Angle> &m_angle ); // MolAngle is established.
        void mergeTwoMol(list<set<unsigned> > & MolIdx, list<set<Bond> > & MolBond, list<set<unsigned> >::reverse_iterator & \
                MolIdx_r_iter, list<set<Bond> >::reverse_iterator & MolBond_r_iter, unsigned beadidx);
        // if iter->b (a bead) is found in another mol, merge these two mol together, because they are the same mol.
        // Exactly, the info of another mol (its Idx and Bond info) will be added into present mol, and that 
        // another mol will be deleted.
        bool IsTopoConsistency(BasicInfo &bi );
        bool IsAllParticleInclude(unsigned bi_na_tot,unsigned pos_size);
        bool IsMolIdxContinuous();
        bool CheckMolBondNum(unsigned m_bond_size );
        void AddBondIdxToBI(BasicInfo& basicinfo); // b_bidx, b_eidx, agl_bidx and agl_eidx, their index need to be modified.

        list<set<unsigned> > MolIdx; 
        list<set<Bond> > MolBond;
        list<set<Angle> > MolAngle;
        vector<MolInfo> Mol;
        vector<string> AvailableTerm; // the info node supported by this program.

};

inline void TopoInfo::MakeOneBeadMolAndInsert(unsigned idx)
{
    set<unsigned> OneMolIdx;
    set<Bond> OneMolBond;
    OneMolIdx.insert(idx);
    MolIdx.push_back(OneMolIdx);
    MolBond.push_back(OneMolBond);
}


#endif

