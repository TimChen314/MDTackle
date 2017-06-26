#include "Traj.h"
#include "BinReader.h"
#include "MolDataStruct.h"
#include "XmlReader.h"
using namespace boost::python;
using namespace boost;

Traj::Traj(string h, unsigned nf, unsigned nfe, unsigned nda)
    : head(h + ".")
    , ndata(nda)
    , nfram(nf)
    , nfram_exc(nfe)
    , framidx(0)
{
    nfram_remain = nfram - nfram_exc;
    MakeFileName();
    OnceRead();
}

void Traj::OnceRead() //get the info that only need once.
{
    getRealCoor(0);
    oncebox = box;
    na_tot = pos.size();
}

unsigned Traj::getNfram_remain()
{
    return nfram_remain;
}

void Traj::MakeFileName()
{
    ostringstream flnam;
    string tail = ".bin";

    flnam << head << setfill('0') << setw(10) << (nfram_exc + 1) * ndata << tail;
    ifstream input_test;
    input_test.open(flnam.str().c_str());
    if (!input_test) {
        tail += ".gz";
        flnam.str("");

        ostringstream flnam2;
        ifstream input_test2;
        flnam << head << setfill('0') << setw(10) << (nfram_exc + 1) * ndata << tail;
        input_test2.open(flnam.str().c_str());

        if (!input_test2) {
            flnam.str("");
            ostringstream flnam3;
            ifstream input_test3;
            tail = ".xml";
            flnam << head << setfill('0') << setw(10) << (nfram_exc + 1) * ndata << tail;
            input_test3.open(flnam.str().c_str());
            if (!input_test3)
                cerr << endl
                     << "***** ERROR!!! First trajectory file does not exist ******" << endl
                     << endl;
            else {
                ft = filetype::xml;
            }
            input_test3.close();
        } else {
            ft = filetype::bin_gz;
        }
        input_test2.close();

    } else {
        ft = filetype::bin;
    }
    flnam.str("");

    for (unsigned int i = nfram_exc + 1; i <= nfram; i++) {
        flnam << head << setfill('0') << setw(10) << i * ndata << tail;
        filename.push_back(flnam.str());
        flnam.str("");
    }
}

void Traj::getRealCoor(unsigned i)
{
    framidx = i;
    if (ft == filetype::bin || ft == filetype::bin_gz) {
        BinReader build(filename[i].c_str());
        pos = build.getRealPos();
        box = build.getBox();
        return;
    } else if (ft == filetype::xml) {
        XmlReader build(filename[i].c_str());
        box = build.getBox();
        pos = build.getRealPos();
    } else {
        cerr << "***** Error! Unknown file type *****" << endl;
        throw runtime_error("Error! Unknown file type ");
    }
}

void Traj::getBoxCoor(unsigned i)
{
    if (ft == filetype::bin || ft == filetype::bin_gz) {
        BinReader build(filename[i].c_str());
        pos = build.getBoxPos();
        box = build.getBox();
        return;
    } else if (ft == filetype::xml) {
        XmlReader build(filename[i].c_str());
        box = build.getBox();
        pos = build.getBoxPos();
    } else {
        cerr << "***** Error! Unknown file type *****" << endl;
        throw runtime_error("Error! Unknown file type ");
    }
}

void export_Traj()
{
    class_<Traj, boost::shared_ptr<Traj> >("Traj", init<const string&, unsigned, unsigned, unsigned>());
}
