
    bool divide_pt = true;
	int trackbin                       = 8;
    const int ptbin                    = 3;
    float     pt_edges[ ptbin +1]      = {0.3, 1, 2, 100};

    const int MS                       = 1;
    const int dEta_N                   = 3;
    int       dEta_bin_color[dEta_N]   = {2,3,4};
    int       dPt_bin_style[ptbin]     = {20,21,22};
    double    dEta_bin_cut[dEta_N]     = {20,30,40};

    double xntrk[trackbin];
    double xntrke[trackbin];

    for(int i = 0; i < trackbin; i++){
    xntrk[i]  =(80/trackbin)*(i+1)-(80/trackbin)/2;
    //xntrke[i] =(80/trackbin)/2;
    xntrke[i] = 0;
    }
