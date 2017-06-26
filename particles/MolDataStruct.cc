#include "MolDataStruct.h"

void FoldBackToBox(vector<vec>& pos, vec& box)
{
    vec boxINV(vec(1.0, 1.0, 1.0) / box);
    for (unsigned i = 0; i != pos.size(); ++i) {
        pbc(box, boxINV, pos[i]);
    }
}

void QuestRealCoorCMAndSetToOri(vector<vec>& pos, vector<vec>& pos_cm, const unsigned nm_tot, const vector<unsigned int>& MolNA)
{
    unsigned int idx = 0;
    for (unsigned i = 0; i < nm_tot; i++) {
        vec sum(0.0, 0.0, 0.0);
        for (unsigned j = 0; j < MolNA[i]; j++) {
            sum += pos[idx];
            idx++;
        }
        pos_cm.push_back(vec(sum.x / MolNA[i], sum.y / MolNA[i], sum.z / MolNA[i]));
    }

    // set to origin
    vec boxcenter(0, 0, 0);
    for (auto v : pos)
        boxcenter += v;
    boxcenter /= pos.size();
    for (auto& vcm : pos_cm)
        vcm -= boxcenter;
    for (auto& v : pos)
        v -= boxcenter;
}

void Vector_vec_To_Vector_float(const vector<vec>& vv, vector<float>& vx, vector<float>& vy, vector<float>& vz)
{
    if (vx.size() != 0 || vy.size() != 0 || vz.size() != 0) {
        cerr << "Error! Vector_vec_To_Vecter in MolDataStruct" << endl;
        throw runtime_error("Info : Vector_vec_To_Vecter in MolDataStruct, none empty vector!");
    }
    vec tmp;
    for (auto v : vv) {
        tmp = v;
        vx.push_back(tmp.x);
        vy.push_back(tmp.y);
        vz.push_back(tmp.z);
    }
}

void DumpAs3202Version(const string fname, const vector<vec>& pos, const vector<vec>& vel, const vec& box, const unsigned timestep, const bool m_output_velocity)
{
    bool m_output_position = true;
    if (vel.size() == 0 && m_output_velocity)
        cerr << endl
             << "***Error! Velocity info does not exist! " << endl
             << endl;
    unsigned int dimensions = 3;
    unsigned int np = pos.size();

    string ext = fname.substr(fname.size() - 4, fname.size());
    if (ext != string(".bin"))
        cerr << "Warning! File " << fname << " is not a bin type file." << endl;

    // setup the file output for compression
    ofstream f(fname);
    if (!f.good()) {
        cerr << endl
             << "***Error! Unable to open dump file for writing: " << fname << endl
             << endl;
        throw runtime_error("Error writing  binary dump file");
    }

    // Begin to dump
    int version = 3202;
    unsigned int magic = 0x444d4f48;
    f.write((char*)&magic, sizeof(unsigned int));
    f.write((char*)&version, sizeof(int));
    f.write((char*)&m_output_position, sizeof(bool));
    f.write((char*)&m_output_velocity, sizeof(bool));
    f.write((char*)&timestep, sizeof(unsigned int));
    f.write((char*)&dimensions, sizeof(unsigned int));
    float Lx = (float)box.x;
    float Ly = (float)box.y;
    float Lz = (float)box.z;
    f.write((char*)&Lx, sizeof(float));
    f.write((char*)&Ly, sizeof(float));
    f.write((char*)&Lz, sizeof(float));
    f.write((char*)&np, sizeof(unsigned int));

    if (m_output_position) {
        vector<float> posx;
        vector<float> posy;
        vector<float> posz;
        Vector_vec_To_Vector_float(pos, posx, posy, posz);
        f.write((char*)&(posx[0]), np * sizeof(float));
        f.write((char*)&(posy[0]), np * sizeof(float));
        f.write((char*)&(posz[0]), np * sizeof(float));
    }
    if (m_output_velocity) {
        vector<float> velx;
        vector<float> vely;
        vector<float> velz;
        Vector_vec_To_Vector_float(vel, velx, vely, velz);
        f.write((char*)&(velx[0]), np * sizeof(float));
        f.write((char*)&(vely[0]), np * sizeof(float));
        f.write((char*)&(velz[0]), np * sizeof(float));
    }

    if (!f.good()) {
        cerr << endl
             << "***Error! Unexpected error writing  dump file" << endl
             << endl;
        throw runtime_error("Error writing  dump file");
    }
}

void ProcessBar(unsigned percent)
{
    if (percent > 100) {
        cerr << "Error in ProcessBar: " << endl;
        throw runtime_error("Error! percentage is out of range.");
    }
    char bar[102];
    for (unsigned i = 0; i != percent; ++i)
        bar[i] = '#';
    bar[percent] = 0;
    const char* lable = "|/-\\";
    printf("[%-100s][%d%%][%c]\r", bar, percent, lable[percent % 4]);
    fflush(stdout);
    if (percent == 100)
        cout << endl;
}
long unsigned getAutoCorrelationProcessPercent(long unsigned ifram, long unsigned nfram_remain)
{
    return (long unsigned)((2 * nfram_remain - 1 - ifram) * ifram * (long unsigned)100 / (nfram_remain * (nfram_remain - 1)));
}

void printOnePercentTime(timeval& start)
{
    timeval end;
    gettimeofday(&end, NULL);
    long timeusr = (end.tv_sec - start.tv_sec) * 1000000 + (end.tv_usec - start.tv_usec);
    ostringstream ss;
    double second = timeusr / 1000000.0;
    ss.precision(4);
    ss << "Info : 1% process time is " << second << " second. Total time is about " << second * 100 << " second.";
    cout << setiosflags(ios::left) << setw(112) << ss.str() << endl;
}
