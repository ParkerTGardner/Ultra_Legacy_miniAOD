void MyScript() {
gROOT->ProcessLine(".L MyClass.C+");
gROOT->ProcessLine("MyClass m");
gROOT->ProcessLine("m.Loop()");
}
