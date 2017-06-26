#include "xmlParser.h"
#include "MolDataStruct.h"
#include "Reader.h"
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
#include <time.h>

#include <iomanip> 
using namespace std;

#include <boost/bind.hpp>
#include <boost/function.hpp>

#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/lexical_cast.hpp>


#ifndef __XML_BULDER_H__
#define __XML_BULDER_H__



class XmlReader : public Reader
{
    public:
        //! Loads in the file and parses the data
        XmlReader(const std::string &fname);

        virtual vec getBox() const;
        virtual unsigned int getNumParticles() const;
        virtual unsigned int getNumParticleTypes() const;
        virtual vector< vec >& getRealPos() ;
        virtual vector<vec >& getBoxPos() ;

        virtual unsigned int getNumDimensions() const;
        virtual unsigned int getTimeStep() const;
        virtual void setTimeStep(unsigned int ts);
        virtual unsigned int getNumBondTypes() const;
        virtual unsigned int getNumAngleTypes() const;
        virtual unsigned int getNumDihedralTypes() const;
        unsigned int getTypeID(const std::string& name);



    private:
        //! Helper function to read the input file
        void readFile(const std::string &fname);
        //! Helper function to parse the box node
        void parseBoxNode(const XMLNode& node);
        //! Helper function to parse the position node
        void parsePositionNode(const XMLNode& node);
        //! Helper function to parse the image node
        void parseImageNode(const XMLNode& node);
        //! Helper function to parse the velocity node
        void parseVelocityNode(const XMLNode& node);
        //! Helper function to parse the mass node
        void parseMassNode(const XMLNode& node);
        //! Helper function to parse diameter node
        void parseDiameterNode(const XMLNode& node);
        //! Helper function to parse the type node
        void parseTypeNode(const XMLNode& node);
        //! Helper function to parse the body node
        void parseBodyNode(const XMLNode& node);		
        //! Helper function to parse the bonds node
        void parseBondNode(const XMLNode& node);
        //! Helper function to parse the angle node
        void parseAngleNode(const XMLNode& node);
        //! Helper function to parse the dihedral node
        void parseDihedralNode(const XMLNode& node);
        //! Helper function to parse the improper node
        void parseMoleIdNode(const XMLNode& node);
        //! Parse charge node
        void parseChargeNode(const XMLNode& node);

        void parseOrientationNode(const XMLNode& node);

        //! Helper function for identifying the bond type id
        unsigned int getBondTypeID(const std::string& name);
        //! Helper function for identifying the angle type id
        unsigned int getAngleTypeID(const std::string& name);
        //! Helper function for identifying the dihedral type id
        unsigned int getDihedralTypeID(const std::string& name);
        //! Helper function for identifying the improper type id
        unsigned int getMoleculerId(const std::string& name);

        std::map< std::string, boost::function< void (const XMLNode& ) > > m_parser_map; //!< Map for dispatching parsers based on node type

        bool m_box_read;    //!< Stores the box we read in



};

void export_XmlReader();

#endif



