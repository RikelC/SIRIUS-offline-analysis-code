#ifndef MyGlobal_h
#define MyGlobal_h 1


#define Si  2   // channel number
#define PMT 1
#define SED 0  //channel number
#define timestamp_unit 10 //ns
#define si_foil_distance 15 //cm
#define pre_trig_buffer 512
#define NBOARDS 12
#define NCHANNELS 16
#define NSTRIPS NBOARDS * NCHANNELS
#define TRACE_SIZE 991
#define sampling_period 5 //ns
#define timestamp_unit 10 //ns
#define Resistance 700. //All the ASIC resistors are 700 k Ohms

#include <sys/types.h>
#include <dirent.h>
#include<unistd.h>
#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
#include <iterator>
#include <cmath>
#include <mutex>
#include "TTree.h"
#include "TChain.h"
#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60
using namespace std;
static unsigned short int NUMEXO2[NBOARDS] = {169, 170, 171, 172, 173, 174, 175, 177, 179, 180, 181, 182};
unsigned short int  getboardIndex(unsigned short int board);
double fit_function(const double *x, const double *p);
double crystalball_function(double x, double alpha, double n, double sigma, double mean);
double crystalball_function(const double *x, const double *p);
void printProgress (unsigned long long int i, TTree* t );
void printProgress (unsigned long long int i, TChain* t );
vector<bool> check_if_files_exist(vector<string> list);
/*
  //Singleton class
class MyGlobal
{
protected:
int value_;
private:
static MyGlobal * instance;
static std::mutex mutex_;
MyGlobal(){}

public:

MyGlobal(MyGlobal &other) = delete;

MyGlobal& operator=(const MyGlobal &) = delete;

static MyGlobal *GetInstance();

double fit_function(const double *x, const double *p);
double crystalball_function(double x, double alpha, double n, double sigma, double mean);
double crystalball_function(const double *x, const double *p);
void printProgress (unsigned long long int i, TTree* t );
void printProgress (unsigned long long int i, TChain* t );
vector<bool> check_if_files_exist(vector<string> list);
};

*/

#endif
