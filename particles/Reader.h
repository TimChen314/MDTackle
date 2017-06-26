#include "MolDataStruct.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <sstream>
#include <algorithm>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <time.h>

#include <iomanip>
using namespace std;

#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/python.hpp>

#ifndef __READER_H__
#define __READER_H__

class Reader : boost::noncopyable
{
    public:
        Reader(){};
        virtual ~Reader() {};

        virtual unsigned int getNumParticles() const = 0;
        /*
         *  This is a way to implement virtual function.
         *  virtual unsigned int getNumParticles() const = 0;
         *  In this way, getNumParticles() becomes a pure virtual function. But you need corresponding change in Reader.cc
         */
        virtual unsigned int getNumParticleTypes() const = 0;
        virtual vec getBox() const = 0;
        virtual vector< vec >& getRealPos() = 0 ;
        virtual vector<vec >& getBoxPos() = 0;


        virtual const std::vector< vec_int >& getImage() const { return m_image_array; }
        virtual const std::vector< vec >& getVel() const { return m_vel_array; }
        virtual const std::vector< Bond >& getBond() const { return m_bonds; }
        virtual const std::vector< Angle >& getAngle() const { return m_angles; }
        virtual const std::vector< unsigned int >& getIdxType() const { return m_idx_type_array; }
        virtual const std::vector< string >& getStringType() const { return m_string_type_array; }
        virtual const std::vector<double>& getMass() const { return m_mass_array; }

        virtual const std::vector< std::string >& getTypeMapping() const { return m_type_mapping; }
        virtual const std::vector< std::string >& getBondTypeMapping() const { return m_bond_type_mapping; }
        virtual const std::vector< std::string >& getAngleTypeMapping() const { return m_angle_type_mapping; }
        virtual const std::vector< std::string >& getDihedralTypeMapping() const { return m_dihedral_type_mapping; }
        virtual const std::set< std::string >& getExistTerm() const { return m_ExistTerm; }
        virtual void printExistTerm() const 
        {
            cout << "Info Reader : "; 
            for(auto &term : m_ExistTerm)
                cout << term  << " " ;
            cout << endl; 

        }

    protected:
        void CheckTerm();

        set<string> m_ExistTerm; //ct: added at 2014.02.23
        vec m_box;   //!< Simulation box read from the file
        unsigned int m_num_dimensions;              //!< number of spatial dimensions
        unsigned int m_timestep;                    //!< The time stamp

        std::vector< vec > m_pos_array;             //!< positions of all particles loaded
        std::vector< vec > m_realpos_array;             //!< positions of all particles loaded
        std::vector< vec_int > m_image_array;       //!< images of all particles loaded
        std::vector< vec > m_vel_array;             //!< velocities of all particles loaded
        std::vector< double > m_mass_array;         //!< masses of all particles loaded
        std::vector< unsigned int > m_idx_type_array;   //!< type values for all particles loaded
        std::vector< string > m_string_type_array;   //!< type values for all particles loaded
        std::vector< Bond > m_bonds;                //!< Bonds read in from the file
        std::vector< Angle > m_angles;              //!< Angle read in from the file
        std::vector< Dihedral > m_dihedrals;        //!< Dihedral read in from the file
        std::vector< vec > m_orientation;             //!< orientation of all particles loaded  
        std::vector< vec > m_accel;

        std::vector< unsigned int > m_body_array;   //!< body values for all particles loaded  
        std::vector< double > m_diameter_array;     //!< diameters of all particles loaded
        std::vector< double > m_charge_array;       //!< charge of the particles loaded
        std::vector< vec4D> m_force_array;           //!< force of all particles loaded
        std::vector< vec4D> m_net_torque_array;  //!< net_torque of all particles loaded
        std::vector< vec4D> m_rotation_array;  //!< rotation of all particles loaded
        std::vector< vec4D> m_orientation_array;  //!< orientation of all particles loaded
        std::vector< double > m_virial_array;           //!< virial of all particles loaded

        std::vector<std::string> m_Molecular_mapping;
        std::vector<unsigned int> m_Molecular_id;
        std::vector<std::string> m_type_mapping;          //!< The created mapping between particle types and ids
        std::vector<std::string> m_bond_type_mapping;     //!< The created mapping between bond types and ids
        std::vector<std::string> m_angle_type_mapping;    //!< The created mapping between angle types and ids
        std::vector<std::string> m_dihedral_type_mapping; //!< The created mapping between dihedral types and ids


};

void export_Reader();

#endif
