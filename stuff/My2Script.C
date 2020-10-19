void My2Script() {
gROOT->ProcessLine(".L My2Class.C");
gROOT->ProcessLine("My2Class k");
gROOT->ProcessLine("k.Loop()");
}
