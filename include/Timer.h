#ifndef TIMER
#define TIMER

#include <iomanip>
#include <chrono>
#include <vector>
#include <string>
#include <iostream>
#include <algorithm>

class Timer{
public:

  Timer();
  void Start();
  void Stop();

  void StartSplit(std::string name );
  void StopSplit();

  void PrintTotalTime();
  void Report();

private:
  bool isRunning = false;
  bool isSplitActive = false;

  int currentSplitIndex = -1;

  std::chrono::high_resolution_clock::time_point start;
  std::chrono::high_resolution_clock::time_point stop;
  std::chrono::high_resolution_clock::time_point splitStart;
  std::chrono::high_resolution_clock::time_point splitStop;

  std::vector< double > totalTimes;
  std::vector< std::string > names;

};

Timer::Timer(){
  totalTimes = std::vector< double >(); 
  names = std::vector< std::string >();
}

void Timer::StartSplit(std::string name){
  auto position = std::find(names.begin(), names.end(), name);
  if (position == names.end())
  {
    names.push_back(name);
    totalTimes.push_back(0);
    position = std::find(names.begin(), names.end(), name);
  }
  currentSplitIndex = position - names.begin();
 
  if(isSplitActive) StopSplit();
  splitStart = std::chrono::high_resolution_clock::now();  
  isSplitActive = true;
}

void Timer::StopSplit(){
  splitStop = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(splitStop - splitStart);
  totalTimes.at(currentSplitIndex) += duration.count()/1000000000.0; 
  isSplitActive = false;
}

void Timer::Start(){
  start = std::chrono::high_resolution_clock::now();
  isRunning = true;
}

void Timer::Stop(){
  if(isSplitActive) StopSplit();
  stop = std::chrono::high_resolution_clock::now();
  isRunning = false;
}

void Timer::PrintTotalTime(){
  if(isRunning){
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>( std::chrono::high_resolution_clock::now() - start);
    std::cout << "Time: " << duration.count()/1000000000.0 << " seconds (but still running!)." << std::endl;
  } else {
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
    std::cout << "Time: " << duration.count()/1000000000.0 <<  " seconds." << std::endl;
  }
}

void Timer::Report(){
  PrintTotalTime();
  std::cout << "Split Times (in seconds)..." << std::endl;
  double total = 0;
  unsigned int maxNameLength = 0;
  for(unsigned int i = 0; i<names.size(); i++){
    total += totalTimes.at(i);
    unsigned int l = strlen(names.at(i).c_str());
    if(l > maxNameLength) maxNameLength = l;
  }
  for(unsigned int i = 0; i<names.size(); i++){
    std::cout << std::left << std::setw(maxNameLength) <<  names.at(i) << " " 
    << std::setprecision(3) << std::setw(10) << totalTimes.at(i) << " (" << std::setprecision(3) << 100*totalTimes.at(i)/total <<" %)" << std::endl;
  }
  std::cout << "Split Total: " << total << "s" << std::endl;
}

#endif
