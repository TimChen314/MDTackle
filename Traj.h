#include "Reader.h"
#ifndef __TRAJ_H__
#define __TRAJ_H__

struct filetype
{
    enum Enum
    {
        xml,
        bin,
        bin_gz
    };
};

class Traj
{
    public:
        Traj(string h, unsigned nf,unsigned nfe,unsigned nda);

        void getBoxCoor(unsigned i);
        void getRealCoor(unsigned i);
        vector<string>& getFilename() {return filename; }
        unsigned getNfram_remain();
        void OnceRead();
        vec& getOnceBox() {return oncebox; }
        vec& getBox() {return box; }
        vector<vec>& getPos() {
            return pos; 
        }
        unsigned getNtot() { return na_tot; }


    private:
        void MakeFileName();
        string head;
        unsigned ndata;
        unsigned nfram;
        unsigned nfram_exc;
        unsigned nfram_remain;
        unsigned na_tot;
        unsigned framidx ;
        double cutoff;
        vec oncebox;
        vec box;
        vector<vec> pos;
        vector<string> filename;
        filetype::Enum ft;

};

void export_Traj();

#endif
