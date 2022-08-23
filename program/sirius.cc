//Offline analysis code for SIRIUS test
//author: Rikel
//date: 01 July 2021
//
//

#include <iostream>
#include "signalTreatment.hh"
#include "TRint.h"
#include "TROOT.h"
#include <string>
#include <vector>
#include "TError.h"
#include "ReadDataFiles.hh"
using namespace std;
vector<string> fileList();
int main(int argc, char *argv[])
{

  TRint *app = new TRint("Analysis", &argc, argv);
  ReadDataFiles *a = new ReadDataFiles();
  gErrorIgnoreLevel = kFatal; 
// gErrorIgnoreLevel = kPrint, kInfo, kWarning, kError, kBreak, kSysError, kFatal;

  //a->set_DSSD_gain(1.5);//in pF

  vector<string> b = fileList();

  unsigned short trapPar[6] = {20,250,10,0,250,10};
  //a->set_Trapezoidal_optimization_limits(trapPar);
    a->perform();
  // a->test(fileList());
  // a->perform(1,fileList());
  //a->get_baseline_noise_info_of("../data/Tree_0051_r0-Numexo2.root");
  // a->drawEntry(10,"../data/Tree_0051_r0-Numexo2.root");
  //a->drawRun2(1,"../data/Tree_0051_r0-Numexo2.root");
  //a->compare_trap_algorithms(0,"../data/Tree_0038_r0-Numexo2.root");
  //a->reOptimize(3,"../data/Tree_0042_r0-Numexo2.root");
  //   a->sort_data_with_min_trapParametes2(3,"../data/Tree_0041_r0-Numexo2.root");
  // a->sort_data_with_min_trapParametesNCalibParameters(3,"../data/Tree_0044_r0-Numexo2.root");
  // test t;
  // cout<<"toto:  "<< t.invoke(20,10,&test::add)<<endl;
   app->Run();
   delete a;
  return 0;
}
vector<string> fileList(){
  vector<string> file;
  //file.push_back("../data/Tree_0049_r0-Numexo2.root");
  //file.push_back("../data/Tree_0041_r0-Numexo2.root");
  //file.push_back("../data/Tree_0041_r1-Numexo2.root");
  //file.push_back("../data/Tree_0044_r0-Numexo2.root");
  //file.push_back("../data/Tree_0043_r0-Numexo2.root");
  return file;
}
