{
  gROOT->ProcessLine("ArtServiceHelper::load(\"mytools.fcl\")");
  gROOT->ProcessLine("DuneToolManager* ptm = DuneToolManager::instance(\"mytools.fcl\")");
  gROOT->ProcessLine(".L getResponseLabel.C");
  gROOT->ProcessLine(".L $WIRECELL_KERNEL_DIR/root/cesmearTH1.C");
  gROOT->ProcessLine(".L CeNoise.h");
  cout << "Tool manager is available at ptm" << endl;
}
