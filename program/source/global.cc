#include "global.hh"

unsigned short int  getboardIndex(unsigned short int board)
{
int n = sizeof(NUMEXO2)/sizeof(NUMEXO2[0]);
 
auto itr = std::find(NUMEXO2, NUMEXO2 + n, board);

return distance(NUMEXO2, itr);
};

vector<bool> check_if_files_exist(vector<string> list){
  //check if the files exist in the directory
  string dir = list[0];
  size_t pos = dir.find_last_of("/");
  size_t len = dir.length();
  dir.erase(pos, len);
  vector<string> files = vector<string>();
  size_t nFiles = list.size();
  vector<bool> flags;
  DIR *dp;
  struct dirent *dirp;

  if((dp = opendir(dir.c_str())) == NULL)
  {
    cerr << "Error opening " << dir << endl;
  }
    
  while ((dirp = readdir(dp)) != NULL)
    files.push_back(string(dirp->d_name));
  closedir(dp);
  for(size_t j = 0; j < nFiles; j++){
    bool status = false;
    pos = list[j].find_last_of("/");
    len = list[j].length();
    for (size_t i = 0;i < files.size();i++){
      if(list[j].compare(pos+1, len, files[i]) ==0){
        //cout << list[j]<<"  exists in"<<dir.c_str()<< endl;
        status = true;
      }  
    }
    flags.push_back(status);
  }
  files.clear();
      
  //char dirpath[100];
  
    // printing current working directory
    //printf("Current directory %s\n", getcwd(dirpath, 100));
  
    // using the command
    //chdir(dir.c_str());
  
    // printing current working directory
    //    printf("New current directory %s\n", getcwd(dirpath, 100));

  return flags;
}

void printProgress (unsigned long long int i, TTree* t )
{
  static int percentage = 0;

  static int noEntries_dividedby100 = t->GetEntries()/100;
  if( percentage > 100) {
    percentage = 0;
    //noEntries_dividedby100 = t->GetEntries()/100;<
  }
  if(i == percentage * noEntries_dividedby100){
    int val = static_cast<int>(percentage);
    int lpad = static_cast<int>(percentage * PBWIDTH / 100);
    int rpad = PBWIDTH - lpad;
    printf ("\rProgress:%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
    fflush (stdout);
    percentage++;
  }

}
void printProgress (unsigned long long int i, TChain* t )
{
  static int percentage = 0;

  static int noEntries_dividedby100 = t->GetEntries()/100;
  if( percentage > 100) {
    percentage = 0;
    //noEntries_dividedby100 = t->GetEntries()/100;<
  }
  if(i == percentage * noEntries_dividedby100){
    int val = static_cast<int>(percentage);
    int lpad = static_cast<int>(percentage * PBWIDTH / 100);
    int rpad = PBWIDTH - lpad;
    printf ("\rProgress:%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
    fflush (stdout);
    percentage++;
  }

}

double crystalball_function(double x, double alpha, double n, double sigma, double mean) {
  // evaluate the crystal ball function
  if (sigma < 0.)     return 0.;
  double z = (x - mean)/sigma; 
  if (alpha < 0) z = -z; 
  double abs_alpha = std::abs(alpha);
  // double C = n/abs_alpha * 1./(n-1.) * std::exp(-alpha*alpha/2.);
  // double D = std::sqrt(M_PI/2.)*(1.+ROOT::Math::erf(abs_alpha/std::sqrt(2.)));
  // double N = 1./(sigma*(C+D));
  if (z  > - abs_alpha)
    return std::exp(- 0.5 * z * z);
  else {
    //double A = std::pow(n/abs_alpha,n) * std::exp(-0.5*abs_alpha*abs_alpha);
    double nDivAlpha = n/abs_alpha;
    double AA =  std::exp(-0.5*abs_alpha*abs_alpha);
    double B = nDivAlpha -abs_alpha;
    double arg = nDivAlpha/(B-z);
    return AA * std::pow(arg,n);
  }
}

double crystalball_function(const double *x, const double *p) {
  if ((!x) || (!p)) return 0.; 
  //par list: [Constant], [Alpha], [N], [Sigma], [Mean]
  return (p[0] * crystalball_function(x[0], p[3], p[4], p[2], p[1]));
}

double fit_function(const double *x, const double *p) {

  return crystalball_function(&x[0], &p[0]) + crystalball_function(&x[0], &p[5]) ;
}

/*
MyGlobal* MyGlobal::instance = 0;
std::mutex MyGlobal::mutex_;
MyGlobal *MyGlobal::GetInstance(){
  //for mutithreading
  std::lock_guard<std::mutex> lock(mutex_);
    if (instance == nullptr)
    {
    instance = new MyGlobal();
    }
    return instance;
}

vector<bool> MyGlobal::check_if_files_exist(vector<string> list){
  //check if the files exist in the directory
  string dir = list[0];
  size_t pos = dir.find_last_of("/");
  size_t len = dir.length();
  dir.erase(pos, len);
  vector<string> files = vector<string>();
  size_t nFiles = list.size();
  vector<bool> flags;
  DIR *dp;
  struct dirent *dirp;

  if((dp = opendir(dir.c_str())) == NULL)
  {
    cerr << "Error opening " << dir << endl;
  }
    
  while ((dirp = readdir(dp)) != NULL)
    files.push_back(string(dirp->d_name));
  closedir(dp);
  for(size_t j = 0; j < nFiles; j++){
    bool status = false;
    pos = list[j].find_last_of("/");
    len = list[j].length();
    for (size_t i = 0;i < files.size();i++){
      if(list[j].compare(pos+1, len, files[i]) ==0){
        //cout << list[j]<<"  exists in"<<dir.c_str()<< endl;
        status = true;
      }  
    }
    flags.push_back(status);
  }
  files.clear();
  
  return flags;
}

void MyGlobal::printProgress (unsigned long long int i, TTree* t )
{
  static int percentage = 0;

  static int noEntries_dividedby100 = t->GetEntries()/100;
  if( percentage > 100) {
    percentage = 0;
    //noEntries_dividedby100 = t->GetEntries()/100;<
  }
  if(i == percentage * noEntries_dividedby100){
    int val = static_cast<int>(percentage);
    int lpad = static_cast<int>(percentage * PBWIDTH / 100);
    int rpad = PBWIDTH - lpad;
    printf ("\rProgress:%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
    fflush (stdout);
    percentage++;
  }

}
void MyGlobal::printProgress (unsigned long long int i, TChain* t )
{
  static int percentage = 0;

  static int noEntries_dividedby100 = t->GetEntries()/100;
  if( percentage > 100) {
    percentage = 0;
    //noEntries_dividedby100 = t->GetEntries()/100;<
  }
  if(i == percentage * noEntries_dividedby100){
    int val = static_cast<int>(percentage);
    int lpad = static_cast<int>(percentage * PBWIDTH / 100);
    int rpad = PBWIDTH - lpad;
    printf ("\rProgress:%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
    fflush (stdout);
    percentage++;
  }

}

double MyGlobal::crystalball_function(double x, double alpha, double n, double sigma, double mean) {
  // evaluate the crystal ball function
  if (sigma < 0.)     return 0.;
  double z = (x - mean)/sigma; 
  if (alpha < 0) z = -z; 
  double abs_alpha = std::abs(alpha);
  // double C = n/abs_alpha * 1./(n-1.) * std::exp(-alpha*alpha/2.);
  // double D = std::sqrt(M_PI/2.)*(1.+ROOT::Math::erf(abs_alpha/std::sqrt(2.)));
  // double N = 1./(sigma*(C+D));
  if (z  > - abs_alpha)
    return std::exp(- 0.5 * z * z);
  else {
    //double A = std::pow(n/abs_alpha,n) * std::exp(-0.5*abs_alpha*abs_alpha);
    double nDivAlpha = n/abs_alpha;
    double AA =  std::exp(-0.5*abs_alpha*abs_alpha);
    double B = nDivAlpha -abs_alpha;
    double arg = nDivAlpha/(B-z);
    return AA * std::pow(arg,n);
  }
}

double MyGlobal::crystalball_function(const double *x, const double *p) {
  if ((!x) || (!p)) return 0.; 
  //par list: [Constant], [Alpha], [N], [Sigma], [Mean]
  return (p[0] * crystalball_function(x[0], p[3], p[4], p[2], p[1]));
}

double MyGlobal::fit_function(const double *x, const double *p) {

  return crystalball_function(&x[0], &p[0]) + crystalball_function(&x[0], &p[5]) ;
}
*/
