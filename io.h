#ifndef IO_H
#define IO_H
#include <iostream>
#include <string>
#include <fstream>      // std::ifstream
#include <stdio.h>  /* defines FILENAME_MAX */
#include <sys/param.h> // for MIN function
#include <stdio.h>  /* defines FILENAME_MAX */
#include <sstream>
#include <vector>
#include <iterator>

template<typename Out>
inline void split(const std::string &s, char delim, Out result) {
    std::stringstream ss;
    ss.str(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        *(result++) = item;
    }
}

inline std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, std::back_inserter(elems));
    return elems;
}

#ifdef WINDOWS
    #include <direct.h>
    #define GetCurrentDir _getcwd
#else
    #include <unistd.h>
    #define GetCurrentDir getcwd
#endif

bool fileExists(std::__cxx11::string fileName);
std::string SplitFilename (const std::string& str);
std::string Find_Program_Installation_Directory();
std::string Find_Program_Working_Directory();

#endif // IO_H
