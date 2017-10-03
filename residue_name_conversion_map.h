#ifndef RESIDUE_NAME_CONVERSION_MAP_H
#define RESIDUE_NAME_CONVERSION_MAP_H
#include <map>
#include "io.h"
#include <fstream>
#include <iostream>
#include <cstdlib> // for exit()

class Residue_Name_Conversion_Map
{
public:

    //////////////////////////////////////////////////////////
    //                    TYPE DEFINITION                   //
    //////////////////////////////////////////////////////////

    typedef std::vector<std::string> StringVector;
    typedef std::vector<std::pair<std::string,std::string>> PairVector;

    //////////////////////////////////////////////////////////
    //                       CONSTRUCTOR                    //
    //////////////////////////////////////////////////////////
    /*! \fn
      * Default constructor
      */
    Residue_Name_Conversion_Map();
    Residue_Name_Conversion_Map(std::string input_file);

    //////////////////////////////////////////////////////////
    //                       ACCESSOR                       //
    //////////////////////////////////////////////////////////

    std::string ConvertViaMap(std::string residue_number);

    //////////////////////////////////////////////////////////
    //                          MUTATOR                     //
    //////////////////////////////////////////////////////////

    void SetInputFile(std::string input_file);
    void ReadInputFile();
    //////////////////////////////////////////////////////////
    //                      DISPLAY FUNCTION                //
    //////////////////////////////////////////////////////////

    void PrintMap();

private:
    //////////////////////////////////////////////////////////
    //                       ATTRIBUTES                     //
    //////////////////////////////////////////////////////////

    std::string input_file_;
    PairVector map_;

};

#endif // RESIDUE_NAME_CONVERSION_MAP_H
