#ifndef Advanced_Header_Ndivide
#define Advanced_Header_Ndivide

#include <iostream>
#include <string>
#include <map>
#include <array>
#include <limits>
#include "TMath.h"


    bool divide_pt = false;
	const int trackbin                 = 6;
    const int ptbin                    = 1;
    float     pt_edges[ ptbin +1]      = {0, 100};

    const int MS                       = 1;
    const int dEta_N                   = 3;
    int       dEta_bin_color[dEta_N]   = {2,3,4};
    int       dPt_bin_style[ptbin]     = {20};
    double    dEta_bin_cut[dEta_N]     = {20,30,40};

    double xntrk[trackbin];
    double xntrke[trackbin];

#endif
