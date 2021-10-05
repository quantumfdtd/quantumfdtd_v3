/*

   paramreader.cpp

   Copyright (c) Michael Strickland

   GNU General Public License (GPLv3)
   See detailed text in license directory

*/

#include <fstream>
#include <iostream>
#include <cstdlib>
#include <complex>
#include <cstring>

using namespace std;

#include "mpisolve.h"

// this workhorse examines a key to see if it corresponds to a var we are setting
// and then attempts to set the var corresponding to key by converting value to the
// appropriate type.  lots of hardcoding here
void setParameter(const char *key, const char *value)
{
    // integer params
    if (strcmp(key, "NUM") == 0)
        NUM = atoi(value);
    if (strcmp(key, "STEPS") == 0)
        STEPS = atoi(value);
    if (strcmp(key, "UPDATE") == 0)
        UPDATE = atoi(value);
    if (strcmp(key, "SNAPUPDATE") == 0)
        SNAPUPDATE = atoi(value);
    if (strcmp(key, "SNAPDUMP") == 0)
        SNAPDUMP = atoi(value);
    if (strcmp(key, "POTENTIAL") == 0)
        POTENTIAL = atoi(value);
    if (strcmp(key, "INITCONDTYPE") == 0)
        INITCONDTYPE = atoi(value);
    if (strcmp(key, "INITCONDAXIS") == 0)
        INITCONDAXIS = atoi(value);
    if (strcmp(key, "INITSYMMETRY") == 0)
        INITSYMMETRY = atoi(value);
    if (strcmp(key, "SAVEWAVEFNCS") == 0)
        SAVEWAVEFNCS = atoi(value);
    if (strcmp(key, "DUMPSNAPS") == 0)
        DUMPSNAPS = atoi(value);
    if (strcmp(key, "KINTERM") == 0)
        KINTERM = atoi(value);
    if (strcmp(key, "SAVEPOT") == 0)
        SAVEPOT = atoi(value);
    if (strcmp(key, "SAVEDECAY") == 0)
        SAVEDECAY = atoi(value);
    if (strcmp(key, "EXTPOT") == 0)
        EXTPOT = strdup(value);
    if (strcmp(key, "DATAFOLD") == 0)
        DATAFOLD = strdup(value);
    if (strcmp(key, "POTCRITR") == 0)
        POTCRITR = atoi(value);
    // double/float params
    if (strcmp(key, "A") == 0)
        A = atof(value);
    if (strcmp(key, "EPS") == 0)
        EPS = atof(value);
    if (strcmp(key, "SIG") == 0)
        SIG = atof(value);
    if (strcmp(key, "SIGMA") == 0)
        SIGMA = atof(value);
    if (strcmp(key, "MASS") == 0)
        MASS = atof(value);
    if (strcmp(key, "TOLERANCE") == 0)
        TOLERANCE = atof(value);
    if (strcmp(key, "POTFLATR") == 0)
        POTFLATR = atof(value);
    return;
}

//
// This routine assumes that paramters are in text file with
// each parameter on a new line in the format
//
// PARAMKEY	PARAMVALUE
//
// The PARAMKEY must begin the line and only tabs and spaces
// can appear between the PARAMKEY and PARAMVALUE.
//
// Lines which begin with 'commentmarker' defined below are ignored
//
void readParametersFromFile(char *filename, int echo)
{

    string commentmarker = "//";
    char space = ' ';
    char tab = '\t';

    const int maxline = 1024; // maximum line length used in the buffer for reading
    char buffer[maxline];
    ifstream paramFile(filename);
    if (!paramFile)
    {
        cerr << "ERROR: THE PARAM FILE CANNOT BE READ!!!" << endl;
        exit(-1);
    }

    while (!paramFile.eof())
    {
        paramFile.getline(buffer, maxline, '\n');
        string line = buffer;
        int length = strlen(buffer);
        if (line.substr(0, commentmarker.length()) != commentmarker && line.length() > 0)
        {
            char key[maxline] = "", value[maxline] = "";
            int founddelim = 0;
            for (int i = 0; i < length; i++)
            {
                if (buffer[i] == space || buffer[i] == tab)
                    founddelim = 1;
                else
                {
                    if (founddelim == 0)
                        key[strlen(key)] = buffer[i];
                    else
                        value[strlen(value)] = buffer[i];
                }
            }
            if (strlen(key) > 0 && strlen(value) > 0)
            {
                setParameter(key, value);
                if (echo)
                    cout << key << " = " << value << endl;
            }
        }
    }

    return;
}

//
// Read parameters from commandline
//
void readParametersFromCommandLine(int argc, char **argv, int echo)
{
    int optind = 1;
    while (optind < argc)
    {
        if (argv[optind][0] == '-')
        {
            string key = argv[optind];
            key = key.substr(1, key.length()); // remove '-'
            string value = argv[optind + 1];   // load value
            if (echo)
                cout << key << " = " << value << endl;
            setParameter(key.c_str(), value.c_str());
            optind++;
        }
        optind++;
    }
    return;
}
