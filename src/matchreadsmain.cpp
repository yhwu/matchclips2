/*************************************************************************
BEGIN OF LICENSE
Copyright (c) Yinghua Wu and Hongzhe Li (matchclips project).

This program is free software; you can redistribute and/or modify
the codes written by the authors under the terms of the GNU General 
Public License as published by the Free Software Foundation 
(www.fsf.org); either version 2 of the License, or (at your option) 
any later version. 

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; See the GNU General Public License for 
more details.

A copy of the GNU General Public License is available at
http://www.fsf.org/licensing/licenses
END OF LICENSE

CONTACT: wu_yinghua@hotmail.com; hongzhe@upenn.edu
*************************************************************************/
/**** system headers ****/
#include <iostream>
#include <time.h>
#include <cstdlib>
#include <unistd.h>
using namespace std;

#include "functions.h"

#ifndef _BUILD_TIME
#define _BUILD_TIME string(__DATE__)+" "+string(__TIME__);
#endif

void match_MS_SM_reads(int argc, char* argv[]);

int main(int argc, char* argv[])
{
  time_t begin_T,end_T;
  int pid = getpid();
  
  string build_time= _BUILD_TIME;
  string updateFile="https://raw.github.com/yhwu/matchclips2/master/src/UPDATED";
  //  check_github_update(build_time, updateFile);
  check_github_update(_BUILD_TIME, 
		      "https://raw.github.com/yhwu/matchclips2/master/src/UPDATED");
  
  exit(0);
  
  time(&begin_T);
  match_MS_SM_reads(argc, argv);  
  
  cerr << procpidstatus(pid,"VmPeak") ;
  time(&end_T);
  cerr << "#Time elapsed: " << difftime(end_T,begin_T) << " seconds\n";
  exit(0);
} 

