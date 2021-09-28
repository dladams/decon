// run.C
// David Adams
// November 2020
//
// This is the script used to generate most of the plots I showed at the
// protoDUNE sim/reco meetings in November 2020.
//
// Charge is deposited uniformly across one wire cell at MIP rate of 60 ke/cm
// with specified angle (thtxzdeg) w.r.t. the wire plane.
// The APA response (charge/tick) for that cell is obtained using the wirecell
// response histograms intepreted with the code in package wirecell-kernel and
// smeared with the CE response from dunetpc.
// There is also the option to add diffusion smearing of 1 or 2 ticks
// roughly corresponding respectively to 1 or 3.6 m drift.
// Noise is added using CeNoise which is random in each tick followed by
// smearing with the CE response function and then scaled to have RMS 1.0.
// This is scaled by the noise level (noise) in ke.
//
// thtxzdeg: horizontal angle w.r.t. the APA plane
// sresp specifies the response simulation
//       = nosig - No charge deposit (noise only)
//       = noconv - Charge deposit only
//       = conv1d:CVAL:CVALS
//       = conv1d:CVALS for 1D convolution
//       = CVALS for 2D convolution
//       CVALS = NSAM:NNBR:NBIN:NSIG:NOISE
//         NSAM = # time samples
//         NBIN = # bins/wire [1]
//         NSIG = longitudinal diffusion Gaussian sigma in ticks [0]
//         NOISE = noise level in e [0]
//         Blank or omitted values take the defaults indicated in brackets.
// sdeco specifies the deconvolution
//   sdeco: nodeco = no deconvolution
//          deco1d:OPT = 1D deconvolution with option OPT
//            OPT = dft - DFT deconvolution
//                  dm is direct matrix deconvolution
//                  fm is filtered matrix deconvolution
//                  cm is chi-square matrix deconvolution
//          deco2d: 2D deconvolution
// sroif specifies the ROI finder: E.g. thr1p for thrSignalFinder1p or noroif for none
// run is the run number and is used to seed the noise
// nevt is the number of events to process (1 track with 3 views for each event)
//
// Products of the script include:
// Waveform plots - One for each event overlaying the three views.
// Power plots - One for each view, before and after reonstruction.

//************************************************************************
// Helpers.
//************************************************************************

int fetchTool(string tnam, TpcDataTool*& ptoo, string prefix) {
  ptoo = ptm->getShared<TpcDataTool>(tnam);
  if ( ptoo == nullptr ) {
    cout << prefix << "ERROR: Unable to find tool " << tnam << endl;
    return 1;
  }
  cout << prefix << "Found tool " << tnam << endl;
  return 0;
}

//************************************************************************
// Main program.
//************************************************************************

int run(string strk, string sresp, string sdeco, string sroif, unsigned int irun =0, unsigned int nevt =1) {
  string myname = "run: ";
  using Index = unsigned int;
  using IndexVector = std::vector<Index>;

  // Define some parameters.
  Index dbg = 1;
  float dzcell = 4.79;         // Cell size in mm
  float dxdt = 0.8;            // Drift speed in mm/tick
  float dqdsMip = 6.0;         // Rate of MIP charge deposit in ke/mm
  Index nsam =   400;          // Number of time samples
  Index rebaselineFlag = 1;    // Rebaselining flag: 0=none, 1=non-signal ticks, 2=tool
  bool mostDistantResponse = false;   // Only for debugging
  float sigThresh = 0.01;   // BG region is evaluated from true signal using this
  Index sigFence = 10;       // threshold and fence. Used for baseline and BG suppresssion
  float samScale = -1.0;   // BG suppression. <0 disables. 0 zeroes. >0 is exp time constant.
  IndexVector sigOff = {0, 10, 5};
  Index nWfPlot = 5;  // # waveform plot for each event: 0, 1, 3 or 5
  float thrFac = 1.0;    // Threshold for area calculations in wf plots is thrFac*noise.

  // Decode the track description:  thtxzdeg:mipfac:wDepMin:wDepMax:wfyfac:nWfPlot
  //   THTXZ - Angle in deg in XZ plane (paralel to plane is 0 deg)
  //   MIPFAC - Rate of energy deposit in MIP [1.0]
  //   WMIN - Wire where track starts (0/1 is left/right edge of central wire) [0.0]
  //   WMAX - Wire where track ends [1.0]
  //   WFYFAC - Range for wf plots is scaled by this factor [1.0]
  float thtxzdeg = 0.0;
  float mipfac =  1.00;         // dQ/ds = dqds * mipfac
  float wDepMin =  0.00;        // Track deposits energy over the wire range
  float wDepMax =  1.00;        // (wDepMin, wDepMax) if max > min.
  float wfyfac = 1.0;
  if ( strk.size() ) {
    vector<string> rsubs = StringManipulator(strk).split(":", true);
    rsubs.resize(6, "");    // Pad with blanks
    if ( rsubs[0].size() ) thtxzdeg = std::stod(rsubs[0]);
    if ( rsubs[1].size() ) mipfac = std::stod(rsubs[1]);
    if ( rsubs[2].size() ) wDepMin = std::stod(rsubs[2]);
    if ( rsubs[3].size() ) wDepMax = std::stod(rsubs[3]);
    if ( rsubs[4].size() ) wfyfac = std::stod(rsubs[4]);
    if ( rsubs[5].size() ) nWfPlot = std::stoi(rsubs[5]);
  }
  string strkLab = "tht" + StringManipulator::floatToString( thtxzdeg, 6, true, "p", "m");
  strkLab +=      "-mip" + StringManipulator::floatToString(   mipfac, 6, true, "p", "m");
  strkLab +=      "-wmn" + StringManipulator::floatToString(  wDepMin, 6, true, "p", "m");
  strkLab +=      "-wmx" + StringManipulator::floatToString(  wDepMax, 6, true, "p", "m");
  ostringstream ssthtxz;
  ssthtxz << std::fixed << std::setprecision(1) << thtxzdeg;
  string sthtxz = ssthtxz.str();
  ostringstream ssmipfac;
  ssmipfac << std::fixed << std::setprecision(1) << mipfac;
  string smipfac = ssmipfac.str();
  ostringstream ssdepRange;
  ssdepRange << "[" << wDepMin << "," << wDepMax << "]";
  string sdepRange = ssdepRange.str();
  float thtxz = thtxzdeg*acos(-1.0)/180.0;
  float dqds = dqdsMip * mipfac;
  cout << myname << "THETA_XZ: " << thtxzdeg << " deg = " << thtxz << " rad" << endl;
  cout << myname << "dQ/ds: " << mipfac << " MIP = " << dqds << " ke/mm" << endl;
  cout << myname << "Wire range: (" << wDepMin << ", " << wDepMax << ")" << endl;
  
  // Build channel labels for the waveform plots.
  std::vector<string> chanLabel(nWfPlot);
  if ( nWfPlot == 1 ) {
    chanLabel[0] = " central";
  } else if ( nWfPlot == 3 ) {
    chanLabel[0] = " left";
    chanLabel[1] = " central";
    chanLabel[2] = " right";
  } else if ( nWfPlot == 5 ) {
    chanLabel[0] = " next left";
    chanLabel[1] = " left";
    chanLabel[2] = " central";
    chanLabel[3] = " right";
    chanLabel[4] = " next right";
  } else {
    cout << myname << "ERROR: Invalid WF plot number: " << nWfPlot << endl;
    return 1;
  }
  cout << myname << "WF plot number: " << nWfPlot << endl;

  // Decode the response sconv:nsam:nnbr:nbinPerWire:ssigl:enoise
  Index nnbr = 0;
  Index nbinPerWire = 1;
  string ssigl = "0";
  Index enoise = 0;
  bool doSig = mipfac != 0.0;
  bool doConvolution = nbinPerWire > 0;
  bool use2dConvolution = true;
  Index badUInt = 9999;
  if ( sresp.size() ) {
    vector<string> rsubs = StringManipulator(sresp).split(":", true);
    Index nsub = rsubs.size();
    string sconv = rsubs[0];
    if ( sconv == "nosig" ) {
      doSig = false;
      mipfac = 0.0;
      doConvolution = false;
    } else if ( sconv == "noconv" ) {
      doConvolution = false;
    } else if ( sconv == "conv1d" ) {
      use2dConvolution = false;
    } else if ( sconv == "conv2d" ) {
      use2dConvolution = true;
    } else {
      cout << myname << "ERROR: Invalid convolution specifier: " << sconv << endl;
      return 1;
    }
    if ( nsub > 1 ) {
      StringManipulator smval(rsubs[1]);
      if ( smval.size() ) nsam = smval.toUnsignedInt(badUInt);
    }
    if ( nsub > 2 ) {
      StringManipulator smval(rsubs[2]);
      if ( smval.size() ) nnbr = smval.toUnsignedInt(badUInt);
    }
    if ( nsub > 3 ) {
      StringManipulator smval(rsubs[3]);
      if ( smval.size() ) nbinPerWire = smval.toUnsignedInt(badUInt);
    }
    if ( nsub > 4 ) ssigl = rsubs[4];
    if ( nsub > 5 ) {
      StringManipulator smval(rsubs[5]);
      if ( smval.size() ) enoise = smval.toUnsignedInt(badUInt);
    }
    if ( sconv == "UnknownConversion" || nnbr == badUInt || nbinPerWire == badUInt || enoise == badUInt ) {
      cout << myname << "ERROR: Invalid response configuration: " << sresp << " --> "
           << sconv << ":" << nnbr << ":" << nbinPerWire << ":" << ssigl << ":" << enoise << endl;
      return 11;
    }
  } else {
    cout << myname << "ERROR: Response configuration must be provided." << endl;
    return 12;
  }
  Index ncha = 2*nnbr + 1;
  bool showNsam = nsam > 0;    // Flag to include nsam in file name.
  cout << myname << "Nsam: " << nsam << endl;
  cout << myname << "Ncha: " << ncha << endl;
  cout << myname << "Convolution: ";
  if ( use2dConvolution ) cout << "2D";
  else if ( doConvolution ) cout << "1D";
  else cout << "none";
  cout << endl;
  if ( doConvolution ) cout << myname << "Bins/wire: " << nbinPerWire << endl;
  cout << myname << "Smearing: " << ssigl << " tick" << endl;
  cout << myname << "RMS noise: " << enoise << " e/tick" << endl;;

  // Get the track slope and the nominal signal widh (ticks).
  // Use the latter to set the signal position isig to center the v-plane signal.
  float dxdz = tan(thtxz);
  float wsignom = dzcell*dxdz/dxdt;
  float wresponse = 20.0;
  float wfac = 5.6;  // Tuned to give full coverage for 3 neighbors at 85 deg.
  if ( nsam == 0 ) nsam = 10*(int(wfac*(wsignom + wresponse))/10);
  Index vplaneOff = 5;
  Index isig = nsam/2 + vplaneOff - 0.5*wsignom;   // Center the signal.

  // Set range for waveform plots
  float wfxmin = isig - 45;    // X-range for waveform plots
  float wfxmax = isig + 45;
  float wfymin = -3.0*wfyfac;  // Y-range for waveform plots
  float wfymax =  8.0*wfyfac;
  if ( false ) {
    wfymax = 50;
    wfymin = -wfymax;
  }
  bool doWideXRange = true;
  if ( doWideXRange ) {
    if ( nsam > 800 ) {
      wfxmin = isig - 400;
      wfxmax = isig + 400;
    } else {
      wfxmin = 0;
      wfxmax = nsam;
    }
  }

  // Decode the deconvolution specifier.
  vector<string> dsubs = StringManipulator(sdeco).split(":", true);
  string dnam = dsubs[0];
  bool do1dDeconvolution = false;
  bool do2dDeconvolution = false;
  vector<string> deconToolNames;
  vector<vector<string>> decon2dToolNames;
  string deconTitle;
  int deconError = 0;
  string sdecon = "No decon";
  if ( dnam == "nodeco" ) {
    cout << myname << "No deconvolution." << endl;
  } else if ( dnam == "deco1d" ) {
    do1dDeconvolution = true;
    string dopt = dsubs.size() == 2 ? dsubs[1] : "BAAAD";
    if ( dopt == "dft" ) {
      cout << myname << "1D DFT deconvolution." << endl;
      deconToolNames.push_back("decon");
      deconTitle = "1D DFT decon of";
      sdecon = "DFT decon";
    } else if ( dopt == "dm" ) {
      cout << myname << "1D direct deconvolution." << endl;
      deconToolNames.push_back("dmdecon");
      deconTitle = "1D Direct matrix decon of";
      sdecon = "DM decon";
    } else if ( dopt == "fm" ) {
      cout << myname << "1D filtered matrix deconvolution." << endl;
      deconToolNames.push_back("fmdecon");
      deconTitle = "1D Filtered matrix decon of";
      sdecon = "FM decon";
    } else if ( dopt == "cm" ) {
      cout << myname << "1D chi-square matrix deconvolution." << endl;
      deconToolNames.push_back("cmdecon");
      deconTitle = "1D Chi-square matrix decon of";
      sdecon = "#chi^{2} decon";
    } else deconError = 2;
  } else if ( dnam == "deco2d" ) {
    sdecon = "2D decon";
    do2dDeconvolution = true;
    deconTitle = "2D decon of";
    if ( dsubs.size() != 1 ) deconError = 3;
    //cout << myname << "ERROR: 2D deconvolution is not yet supported." << endl;
    //return 13;
    decon2dToolNames.resize(3);
    cout << myname << "2D DFT deconvolution." << endl;
    decon2dToolNames[0].push_back("roi2d");
    decon2dToolNames[0].push_back("decon2dx");
    decon2dToolNames[0].push_back("roiToAdc");
    decon2dToolNames[1].push_back("roi2d");
    decon2dToolNames[1].push_back("decon2du");
    decon2dToolNames[1].push_back("roiToAdc");
    if ( rebaselineFlag == 2 ) decon2dToolNames[1].push_back("rebaseline");
    decon2dToolNames[2].push_back("roi2d");
    decon2dToolNames[2].push_back("decon2dv");
    decon2dToolNames[2].push_back("roiToAdc");
    if ( rebaselineFlag == 2 ) decon2dToolNames[2].push_back("rebaseline");
  } else deconError = 1;
  if ( deconError ) {
    cout << "Invalid deconvolution option: " << sdeco << endl;
    return 14;
  }
  bool doDeconvolution = do1dDeconvolution || do2dDeconvolution;
  string deconName;
  for ( string dsub : dsubs ) deconName += dsub;

  // Decode the signal finder specifier.
  dsubs = StringManipulator(sroif).split(":", true);
  dnam = dsubs[0];
  string sroifTool = "roifError";
  string sroifLabel = "???";
  if ( dnam == "decon" ) {
    sroifTool = "deconSignalFinder";
    sroifLabel = "deconSF";
  } else if ( dnam == "nodecon" ) {
    sroifTool = "nodeconSignalFinder";
    sroifLabel = "nodeconSF";
  } else if ( dnam == "noROI" ) {
    sroifLabel = "no SF";
    sroifTool = "thrSignalFinderNone";
  } else if ( dnam.substr(0, 3) == "thr" ) {
    string sthr = dnam.substr(3);
    string sfthr = sthr;
    // Replace '.' with 'p'.
    for ( string::size_type ipos=0; ipos<sfthr.size(); ++ipos ) {
      if ( sfthr[ipos] == '.' ) sfthr[ipos] = 'p';
    }
    sroifLabel = "thr=" + sthr;
    sroifTool = "thrSignalFinder" + sfthr;
  } else {
    cout << myname << "ERROR: Invalid ROI finder specifier: " << dnam << endl;
    return 1;
  }

  // Derived parameters.
  float noise = 0.001*enoise;
  int ncel = nnbr;
  vector<int> apaConCells; // Cells contributing to central signal.
  if ( mostDistantResponse ) {
    if ( ncel > 0 ) apaConCells.push_back(-ncel);
    apaConCells.push_back( ncel);
  } else {
    for ( int icel=-ncel; icel<=ncel; ++icel ) apaConCells.push_back(icel);
  }

  string responseTitle;
  if ( ! doConvolution ) {
    cout << myname << "No convolution." << endl;
    responseTitle = "Response from";
    nbinPerWire = 1;
  } else if ( doDeconvolution ) {
    responseTitle = deconTitle;
  } else if ( ! use2dConvolution ) {
    cout << myname << "1D convolution." << endl;
    responseTitle = "1D convolution of";
  } else {
    cout << myname << "Binned 2D convolution." << endl;
    responseTitle = "Binned 2D convolution of";
  }

  cout << myname << "Neighbor count: " << ncel << endl;
  cout << myname << "# contributing cells: " << apaConCells.size() << endl;
  cout << myname << "# bins/wire: " << nbinPerWire << endl;
  cout << myname << "# samples: " << nsam << endl;
  cout << myname << "# channels " << ncha << endl;
  cout << myname << "Signal tick: " << isig << endl;

  // Fetch noise scaled to 1.0 for this run.
  CeNoise cen(irun);

  // Create global label and name string.
  ostringstream sslab;
  string space = "  ";
  sslab << thtxzdeg << "#circ";
  sslab << space;
  if ( doSig ) {
    sslab << mipfac << " MIP";
    sslab << space;
    sslab << "[" << wDepMin << ", " << wDepMax << "]";
  } else sslab << "Noise only";
  sslab << space;
  sslab << nsam << "#times" << ncha << "/" << nbinPerWire;
  sslab << space;
  sslab << "#sigma_{#lower[-0.5]{l}}=" << ssigl;
  sslab << space;
  sslab << enoise << " e";
  sslab << space;
  sslab << sdecon;
  sslab << space;
  sslab << sroifLabel;
  string sglobalLabel = sslab.str();

  // Fetch the tools.
  // Tools that use Root objects in their dtors are created as private and deleted
  // here to ensure those objects are still present.
  string line = "************************************************************";
  cout << myname << line << endl;
  cout << myname << "Fetching tools." << endl;
  cout << myname << line << endl;
  cout << myname << "Fetching waveform plotter." << endl;
  TpcDataTool* pwf = ptm->getShared<TpcDataTool>("wf");
  cout << myname << "Fetching 1D reponse convolution tools." << endl;
  TpcDataTool* pconvo = ptm->getShared<TpcDataTool>("convo");
  TpcDataTool* pconvo1 = ptm->getShared<TpcDataTool>("convo1");
  TpcDataTool* pconvo2 = ptm->getShared<TpcDataTool>("convo2");
  TpcDataTool* pconvo3 = ptm->getShared<TpcDataTool>("convo3");
  TpcDataTool* pconvo4 = ptm->getShared<TpcDataTool>("convo4");
  cout << myname << "Fetching 2D response convolution tools." << endl;
  map<Index, map<string, TpcDataTool*>> pcon2ds;
  for ( Index nbin : {1, 5, 10, 20} ) {
    for ( string sview : {"x", "u", "v"} ) {
      string tnam = "conv2dNbin" + to_string(nbin) + sview;
      if ( fetchTool(tnam, pcon2ds[nbin][sview], myname) ) return 16;
    }
  }
  vector<TpcDataTool*> deconTools(deconToolNames.size(), nullptr);
  cout << myname << "Fetching deconvolution tool(s)." << endl;
  auto idt = deconTools.begin();
  for ( string deconToolName : deconToolNames ) {
    if ( fetchTool(deconToolName, *(idt++), myname) ) return 17;
  }
  cout << myname << "Fetching 2D deconvolution tool(s)." << endl;
  vector<vector<TpcDataTool*>> decon2dTools(decon2dToolNames.size());
  for ( Index ipla=0; ipla<decon2dToolNames.size(); ++ipla ) {
    decon2dTools[ipla].resize(decon2dToolNames[ipla].size());
    auto idt = decon2dTools[ipla].begin();
    for ( string deconToolName : decon2dToolNames[ipla] ) {
      cout << myname << "...fetching " << deconToolName << endl;
      if ( fetchTool(deconToolName, *(idt++), myname) ) return 17;
    }
  }
  cout << myname << "Fetching longitudinal diffusion smearing tool." << endl;
  TpcDataTool* pconvg = nullptr;
  string slabsmr = "no smearing";
  if ( ssigl == "1" || ssigl == "2" ) {
    string sconvg = "convg" + ssigl;
    cout << myname << "Fetching filter smearing tool " << sconvg << endl;
    pconvg = ptm->getShared<TpcDataTool>(sconvg);
    if ( pconvg == nullptr ) {
      cout << myname << "ERROR: Unable to find tool " << sconvg << endl;
      return 18;
    }
    slabsmr = "#sigma_{smear} = " + ssigl + " ticks";
  } else if ( ssigl != "0" ) {
    cout << myname << "ERROR: Invalid ssigl: " << ssigl << endl;
    return 19;
  }
  cout << myname << "Fetching DFT before plotter." << endl;
  auto pplotDftBefore = ptm->getPrivate<TpcDataTool>("plotDftBefore");
  if ( ! pplotDftBefore ) return 20;
  cout << myname << "Fetching DFT after plotter." << endl;
  auto pplotDftAfter  = ptm->getPrivate<TpcDataTool>("plotDftAfter");
  if ( ! pplotDftAfter ) return 21;
  cout << myname << "Fetching FFT calculator." << endl;
  TpcDataTool* pfft = ptm->getShared<TpcDataTool>("adcFFT");
  if ( pfft == nullptr ) return 22;
  cout << myname << "Fetching signal finder " << sroifTool << endl;
  TpcDataTool* psigfind = sroifTool.size() ? ptm->getShared<TpcDataTool>(sroifTool) : nullptr;
  if ( sroifTool.size() && psigfind == nullptr ) {
    cout << myname << "Unable to find tool " << sroifTool << endl;
    return 23;
  }
  cout << myname << "Fetching sample RMS calculator." << endl;
  TpcDataTool* psamrms = ptm->getShared<TpcDataTool>("adcChannelSamplelRmsPlotter");
  if ( psamrms == nullptr ) return 24;
  cout << myname << "Fetching signal RMS calculator." << endl;
  TpcDataTool* psigrms = ptm->getShared<TpcDataTool>("adcChannelSignalRmsPlotter");
  if ( psigrms == nullptr ) return 25;
  cout << myname << "Fetching not-signal RMS evaluators." << endl;
  vector<TpcDataTool*> rmsEvaluators;
  TpcDataTool* prmsEvaluator = nullptr;
  for ( string sint : {"SigFrac", "", "10", "20", "30", "40", "50", "60"} ) {
    string stvar = sint;
    if ( sint[0] != 'S' ) stvar = "NotSignalRms" + sint;
    string stnam = "adcChannel" + stvar + "Plotter";
    cout << myname << "... " << stnam << endl;
    TpcDataTool* pnsgrms = ptm->getShared<TpcDataTool>(stnam);
    if ( pnsgrms == nullptr ) {
      cout << myname << "Tool not found: " << stnam << endl;
      return 25;
    }
    if ( sint == "" ) prmsEvaluator = pnsgrms;
    rmsEvaluators.push_back(pnsgrms);
  }
  cout << myname << "Fetching rebaseline tool." << endl;
  TpcDataTool* prbl = ptm->getShared<TpcDataTool>("rebaseline");
  //TpcDataTool* prbl = ptm->getShared<TpcDataTool>("adcPedestalFit");
  if ( prbl == nullptr ) return 26;
  cout << myname << "Fetching ROI viewer." << endl;
  string rvname = "roiViewer";
  if ( thtxzdeg > 80 ) rvname += "Wide2";
  else if ( thtxzdeg > 70 ) rvname += "Wide";
  auto proiview = ptm->getPrivate<TpcDataTool>(rvname);
  if ( proiview == nullptr ) return 27;
  // ROI tree maker.
  TpcDataTool* proiTreeMaker = ptm->getShared<TpcDataTool>("adcRoiTreeMaker");
  // End fetching tools.

  // Create charge deposition.
  cout << myname << line << endl;
  // Charge is deposited in time-wire bins.
  // We integrate across each 0.5 us time bin starting at the left edge of bin 0.
  string sttlPrefix;   // Used to build title for waveform plots.
  map<int, vector<float>> binsamples;     // Q[ibin][isam]
  map<int, vector<float>> celsamples;     // Q[icel][isam]
  float qdeptot = 0.0;
  if ( doSig ) {
    cout << myname << "Depositing charge." << endl;
    cout << myname << "Adding signal." << endl;
    //   w is wire coordinate in wire spacing units
    //   t is drift coordinate in ticks (not time or distance)
    //   wbWid = wire bin width [wire spacings]
    //   tbWid = time bin width [ticks]
    //   Wire bin 0 is at the left edge of the central cell.
    //   Time bin 0 is tick 0
    //   (wtoff, ttoff) is a point on the track
    //   dtdw derived from dxdz is the track angle
    int iwir1 = -ncel*nbinPerWire;
    int iwir2 = (ncel + 1)*nbinPerWire;
    for ( int iwir=iwir1; iwir<iwir2; ++iwir ) {
      binsamples[iwir].resize(nsam, 0.0);
    }
    sttlPrefix = "";
    float wbWid = 1.0/nbinPerWire;
    float tbWid = 1.0;
    float wbOff = 0.0;
    float dtdw = dxdz*dzcell/dxdt;    // dimensionless
    float dzds = cos(thtxz);
    float dxds = sin(thtxz);
    float dqdw = dqds*dzcell/dzds;
    float dqdt = dqds*dxdt/dxds;
    float wtOff = 0.0;
    float ttOff = isig;
    // Loop over wire bins.
    for ( auto& kwir : binsamples ) {
      int iwir = kwir.first;
      if ( dbg >= 2 ) cout << myname << "Wire bin " << iwir << endl;
      vector<float>& wsamples = kwir.second;
      // Check what fraction of track is depositing energy.
      float xwir1 = wbWid*iwir;
      float xwir2 = xwir1 + wbWid;
      if ( wDepMax > wDepMin ) {
        if ( xwir1 < wDepMin ) xwir1 = wDepMin;
        if ( xwir2 > wDepMax ) xwir2 = wDepMax;
      }
      if ( xwir2 <= xwir1 ) continue;
      // Handle track parallel to APA
      if ( dxdz == 0.0 ) {
        float fdep = (xwir2 - xwir1)/wbWid;
        float qdep = fdep*dqdw*wbWid;
        wsamples[isig] = qdep;
        qdeptot += qdep;
        if ( dbg >= 2 )cout << myname << "  Q[" << iwir << "][" << isig << "] = " << qdep << endl;
        continue;
      }
      float tmin = ttOff + dtdw*(iwir)*wbWid;
      //if ( tmin < tDepMin ) tmin = tDepMin;
      float tmax = ttOff + dtdw*(iwir+1)*wbWid;
      //if ( tmax > tDepMax ) tmax = tDepMax;
      int it1 = tmin - 1.0;
      int it2 = tmax + 1.0;
      // Loop over ticks.
      for ( int it=it1; it<=it2; ++it ) {
        if ( it < 0 ) {
          cout << myname << "WARNING: Skipping tick " << it << " for cell " << iwir << endl;
          continue;
        }
        if ( it > nsam ) {
          cout << myname << "WARNING: Skipping tick " << it << " for cell " << iwir << endl;
          continue;
        }
        float t1 = it > tmin ? it : tmin;
        float t2 = it+1 < tmax ? it+1 : tmax;
        if ( t2 > t1 ) {
          float qdep = dqdt*(t2 - t1);
          wsamples[it] = qdep;
          qdeptot += qdep;
          if ( dbg >= 2 ) cout << myname << "  Q[" << iwir << "][" << it << "] = " << qdep << endl;
        }
      }
    }  // End loop over wire bins.
    // If we are not doing 2D convolution or not doing convolution at all, then construct
    // cell-level data from bin-level data.
    if ( !doConvolution || !use2dConvolution ) {
      cout << myname << "Filling samples from sample bins." << endl;
      for ( int icel=-ncel; icel<=ncel; ++icel ) {
        celsamples[icel].resize(nsam, 0.0);
      }
      if ( binsamples.size() ) {
        for ( int icel=-ncel; icel<=ncel; ++icel ) {
          if ( find(apaConCells.begin(), apaConCells.end(), icel) == apaConCells.end() ) continue;
          int ibin = icel*nbinPerWire;
          vector<float>& dest = celsamples[icel];
          for ( int iwbi=0; iwbi<nbinPerWire; ++iwbi, ++ibin ) {
            cout << myname << "Adding bin " << ibin << " to cell " << icel << endl;
            std::transform(dest.begin(), dest.end(), binsamples[ibin].begin(), dest.begin(), std::plus<float>());
          }
        }
      }
    } else {
      cout << myname << "Not filling samples from sample bins." << endl;
    }
  }
  cout << myname << "Total deposited charge: " << qdeptot << " e" << endl;
  // Display bin samples.
  if ( dbg >= 2 ) {
    cout << myname << "Bin sample counts, charge: " << endl;
    for ( auto& kwir : binsamples ) {
      int iwir = kwir.first;
      float qsum = 0.0;
      for ( float q : binsamples[iwir] ) qsum += q;
      cout << myname << "Wire bin " << setw(4) << iwir << " has " << binsamples[iwir].size()
           << " samples with total charge " << qsum << " e." << endl;
    }
  }
  if ( celsamples.size() ) {
    cout << myname << "Cell sample counts: " << endl;
    for ( auto& kwir : celsamples ) {
      int iwir = kwir.first;
      float qsum = 0.0;
      for ( float q : celsamples[iwir] ) qsum += q;
      cout << myname << setw(6) << iwir << ": " << celsamples[iwir].size() << ", " << qsum << endl;
    }
  } else {
    cout << myname << "Cell sample array has not been filled." << endl;
  }

  // Assign channel numbers, colors, styles, etc. for the three views.
  // Note channel number will be offset by the track angle.
  LineColors lc;
  vector<int> chans = {2100, 100, 1100};
  vector<int> mycols = {lc.blue(), lc.red(), lc.green()};
  vector<int> mywids = {2, 2, 4};
  vector<int> mystys = {1, 2, 3};
  vector<string> mylabs = {"X", "U", "V"};
  vector<string> views = {"x", "u", "v"};

  // Loop over the three views unless we are only showing the deposition.
  unsigned int npla = doConvolution||doDeconvolution ? 3 : 1;

  // Assign the central channel for each plane.
  // Also create and fill the data map used when closing the roi viewer.
  IndexVector centralChannels;
  DataMap dmRoiview;
  for ( unsigned int ipla=0; ipla<npla; ++ipla ) {
    Index icha = chans[ipla] + thtxzdeg;
    centralChannels.push_back(icha);
    string scha	= to_string(icha);
    while ( scha.size() < 5 ) scha = "0" + scha;
    string resHistName = "hsa_ch" + scha;
    dmRoiview.setFloat(resHistName, qdeptot);
  }

  // Loop over events.
  vector<TH1*> myhsts(npla, nullptr);
  Index maxEvtPlot = 20;    // # events for which waveform plots are made
  unique_ptr<TPadManipulator> pmanWfs;
  for ( Index ievt=1; ievt<nevt+1; ++ievt ) {
    cout << myname << line << endl;
    cout << myname << "######## Event " << ievt << endl;
    cout << myname << line << endl;
    // Create event info.
    DuneEventInfo evt;
    evt.run = irun;
    evt.event = ievt;
    cout << myname << "Noise counts: " << cen.count() << ", " << cen.next()
         << ".  Last = " << cen.last() << endl;
    // Create pad for waveform plots.
    if ( ievt < maxEvtPlot ) {
      if ( nWfPlot > 0 ) pmanWfs.reset(new TPadManipulator(nWfPlot*700,500));
      if ( nWfPlot > 1 ) pmanWfs->split(nWfPlot, 1);
    } else pmanWfs.reset();
    vector<vector<string>> slabsWf(nWfPlot);
    for ( unsigned int ipla=0; ipla<npla; ++ipla ) {
      TpcData tpd(1);
      Index icha = centralChannels[ipla];
      string sview = views[ipla];
      // Create ADC channel data with the energy deposition
      AdcChannelDataMap& acm = *tpd.getAdcData()[0];
      vector<Index> ichas;
      ichas.push_back(icha);
      for ( int icel : apaConCells ) {
        Index jcha = Index(int(icha) + icel);
        if ( jcha == icha ) continue;
        ichas.push_back(jcha);
      }
      cout << myname << "Creating " << ichas.size() << " channels." << endl;
      Index ifmb = 0;
      Index ifch = 0;
      Index chanStat = 0;
      for ( Index jcha : ichas ) {
        AdcChannelData& jacd = acm[jcha];
        jacd.setEventInfo(irun, ievt);
        jacd.sampleUnit = "ke/tick";
        jacd.setChannelInfo(jcha, ifmb, ifch, chanStat);
        jacd.samples.resize(nsam);
      }
      AdcChannelData& acd = acm[icha];
      IndexVector plotChannels;    // Channels to be plotted.
      if ( ! doConvolution ) {
        // Fill samples for the single cell.
        for ( Index isam=0; isam<celsamples[0].size(); ++isam ) {
          if ( isam >= nsam ) break;
          acd.samples[isam] = celsamples[0][isam];
        }
        //if ( doConvolution ) {
        //  cout << myname << "Convoluting with 1D detector response " << icha << endl;
        //  pconvo->update(acd).print();
        //}
        if ( pconvg != nullptr ) {
          cout << myname << "Adding longitudinal diffusion to single channel " << icha << endl;
          pconvg->update(acd).print();
        }
      } else if ( ! doSig ) {
        cout << myname << "ERROR: There is no signal to convolute." << endl;
        return 23;
      } else if ( use2dConvolution ) {
        // Fill bin samples for the signal cells.
        cout << myname << "Convoluting in 2D" << endl;
        // Fill data bins from energy deposit bins.
        for ( int icel : apaConCells ) {
          Index jcha = Index(int(icha) + icel);
          AdcChannelData& jacd = acm[jcha];
          jacd.samples.clear();
          jacd.binSamples.clear();
          jacd.binSamples.resize(nbinPerWire);
          for ( Index jbin=0; jbin<nbinPerWire; ++jbin ) {
            int ibin = icel*nbinPerWire + jbin;
            if ( binsamples[ibin].size() != nsam ) {
              cout << myname << "ERROR: Energy deposit bin has unexpected size." << endl;
              return 28;
            }
            //cout << myname << "  Filling chan/bin " << jcha << "/" << jbin << " from edep bin " << ibin << endl;
            vector<float>& dest = jacd.binSamples[jbin];
            dest.resize(nsam);
            jacd.binSamples[jbin].resize(nsam);
            for ( Index isam=0; isam<nsam; ++isam ) dest[isam] = binsamples[ibin][isam];
          }  // End loop over channel bins
        }
        if ( pcon2ds.count(nbinPerWire) == 0 ) {
          cout << myname << "Tool for 2d convolution with " << nbinPerWire << " bins not found." << endl;
          return 29;
        }
        TpcDataTool* pcon2d = pcon2ds[nbinPerWire][sview];
        pcon2d->updateMap(acm).print();
      } else {
        // Fill samples for the signal cells.
        cout << myname << "Convoluting with neighbors" << endl;
        for ( int icel : apaConCells ) {
          cout << myname << "...cell " << icel << endl;
          TpcDataTool* pcon = abs(icel) == 0 ? pconvo  :
                                  abs(icel) == 1 ? pconvo1 :
                                  abs(icel) == 2 ? pconvo2 :
                                  abs(icel) == 3 ? pconvo3 :
                                  abs(icel) == 4 ? pconvo4 : nullptr;
          if ( pcon == nullptr ) {
            cout << myname << "WARNING: Convolution tool not found for cell " << icel << endl;
            continue;
          }
          AdcChannelData acdtmp;
          acdtmp.setChannelInfo(icha + icel, ifmb, ifch + icel, chanStat);
          acdtmp.samples = celsamples[icel];
          if ( pconvg != nullptr ) {
            cout << myname << "Adding longitudinal diffusion to channel " << acdtmp.channel() << endl;
            pconvg->update(acdtmp).print();
          }
          pcon->update(acdtmp);
          for ( Index isam=0; isam<acdtmp.samples.size(); ++isam ) acd.samples[isam] += acdtmp.samples[isam];
        }  // End loop over neighbors
      }
      // Find the truth BG region for each channel.
      using IndexPair = std::pair<Index, Index>;
      map<Index, vector<IndexPair>> chanBgRanges;
      for ( auto& kacdPair : acm ) {
        AdcChannelData& kacd = kacdPair.second;
        Index kcha = kacd.channel();
        chanBgRanges[kcha];
        Index isam1=0;
        Index isam2 = isam1;
        while ( isam2 < nsam && fabs(kacd.samples[isam2]) < sigThresh ) ++isam2;
        Index dsig = isam2 < nsam ? sigFence : 0;
        if ( isam2 - isam1 > dsig ) chanBgRanges[kcha].emplace_back(isam1, isam2 - dsig);
        if ( isam2 < nsam ) {
          isam2 = nsam;
          isam1 = isam2;
          while ( isam1 > 0 && fabs(kacd.samples[isam1-1]) < sigThresh ) --isam1;
          dsig = sigFence + sigOff[ipla];
          if ( isam2 - isam1 > dsig ) {
            if ( chanBgRanges[kcha].size() == 0 ) chanBgRanges[kcha].emplace_back(0, 0);
            chanBgRanges[kcha].emplace_back(isam1 + dsig, isam2);
          }
        }
        if ( dbg >= 1 ) {
          bool first = true;
          cout << myname << "BG ranges for channel " << setw(4) << kcha << ": ";
          for ( IndexPair ran : chanBgRanges[kcha] ) {
            if ( first ) first = false;
            else cout << ", ";
            cout << "[" << ran.first << ", " << ran.second << ")";
          }
        }
        cout << endl;
      }
      // Add noise.
      if ( noise > 0.0 ) {
        cout << myname << "Adding noise " << icha << endl;
        for ( auto& kacdPair : acm ) {
          AdcChannelData& kacd = kacdPair.second;
          for ( float& sam : kacd.samples ) sam += noise*cen.get();
          kacd.dftmags.clear();
          kacd.dftphases.clear();
        }
      }
      // Zero BG regions.
      if ( samScale >= 0.0 ) {
        if ( noise > 0.0 ) {
          cout << myname << "Suppressing noise " << icha << endl;
          for ( auto& kacdPair : acm ) {
            AdcChannelData& kacd = kacdPair.second;
            Index kcha = kacd.channel();
            bool left = true;
            Index nsup = 0;
            for ( IndexPair ran : chanBgRanges[kcha] ) {
              for ( Index isam=ran.first; isam<ran.second; ++isam ) {
                float fac = 0.0;
                if ( samScale && ran.second - ran.first < nsam ) {
                  Index dsam = left ? ran.second - isam : isam - ran.first;
                  fac = exp(-float(dsam)/samScale);
                }
                kacd.samples[isam] *= fac;
                ++nsup;
              }
              left = false;
            }
            cout << myname << "Channel " << kcha << " suppressed sample count: " << nsup << "/" << nsam << endl;
          }
        }
      }
      if ( true ) {
        pfft->update(acd);
        pplotDftBefore->viewMap(acm);
        // Deconvolution.
        if ( do2dDeconvolution ) {
          cout << myname << "Deconvoluting in 2D." << endl;
          for ( TpcDataTool* pdecon : decon2dTools[ipla] ) pdecon->updateTpcData(tpd).print();
          pfft->update(acd);
          pplotDftAfter->viewMap(acm);
        } else if ( doDeconvolution ) {
          cout << myname << "Deconvoluting channel " << icha << endl;
          for ( TpcDataTool* pdecon : deconTools ) pdecon->update(acd).print();
          if ( rebaselineFlag == 2 ) {
            cout << myname << "Rebaselining channel " << icha << " with tool." << endl;
            prbl->update(acd).print();
          }
          pfft->update(acd);
          pplotDftAfter->viewMap(acm);
        } else {
          cout << myname << "Skipping deconvolution of channel " << icha << endl;
        }
        tpd.print(myname + "TpcData: ");
        if ( rebaselineFlag == 1 ) {
          cout << myname << "Rebaselining with truth BG regions." << endl;
          for ( auto& kacdPair : acm ) {
            AdcChannelData& kacd = kacdPair.second;
            Index kcha = kacd.channel();
            Index nsamtot = 0;
            float samtot = 0.0;
            for ( IndexPair ran : chanBgRanges[kcha] ) {
              for ( Index isam=ran.first; isam<ran.second; ++isam ) {
                samtot += kacd.samples[isam];
                ++nsamtot;
              }
            }
            if ( nsamtot ) {
              float samped = samtot/float(nsamtot);
              for ( float& sam : kacd.samples ) sam -= samped;
            } else {
              cout << myname << "WARNING: Channel << kcha is not rebaselined." << endl;
            }
          }
        }
      }
      // Find ROI and update ROI/noise histograms.
      if ( psigfind != nullptr ) {
        cout << myname << "Finding signal" << endl;
        psigfind->updateMap(acm).print();
      }
      cout << myname << "Updating ROI histograms" << endl;
      proiview->view(acd).print();
      cout << myname << "Updating ROI tree" << endl;
      float sampleRmsAll = 0.0;
      vector<float> wfPlotNoiseLevels(nWfPlot, 0.0);
      for ( Index itoo=0; itoo<rmsEvaluators.size(); ++itoo ) {
        TpcDataTool* ptoo = rmsEvaluators[itoo];
        DataMap dm = ptoo->updateMap(acm);
        if ( dm.getString("metricName") == "nsgRms" ) {
          const DataMap::IntVector& chas = dm.getIntVector("metricChannels");
          const DataMap::FloatVector& vals = dm.getFloatVector("metricValues");
          if ( vals.size() != chas.size() ) return 1;
          Index kchaMet = 0;
          for ( Index iplt=0; iplt<nWfPlot; ++iplt ) {
            Index ichaPlot = icha - nWfPlot/2 + iplt;
            while ( chas[kchaMet]<ichaPlot && kchaMet<chas.size() ) ++kchaMet;
            if ( kchaMet == chas.size() ) {
              cout << myname << "Noise metric nsgRms does not have channel " << ichaPlot << endl;
              return 1;
            }
            wfPlotNoiseLevels[iplt] = vals[kchaMet];
            ++kchaMet;
          }
        }
        dm.print();
      }
      proiTreeMaker->viewMap(acm).print();
      DataMap dsamrmsNew = prmsEvaluator->view(acd);
      float sampleRmsNew = dsamrmsNew.getFloat("metricValue");
      cout << myname << "Evaluating background RMS" << endl;
      // Invert signal finder to get background noise level.
      for ( Index isam=0; isam<acd.signal.size(); ++isam ) acd.signal[isam] = acd.signal[isam] ? false : true;
      DataMap dsamrmsOld = psamrms->view(acd);
      float sampleRmsOld = dsamrmsOld.getFloat("metricValue");
      cout << myname << "Old sample sample RMS: " << sampleRmsOld << endl;
      cout << myname << "New signal sample RMS: " << sampleRmsNew << endl;
      DataMap dsigrms = psigrms->view(acd);
      cout << myname << "Signal RMS: " << dsigrms.getFloat("metricValue") << endl;
      //dsamrms.print();
      // Fetch the waveform and add it to the waveform plot.
      TH1* ph = nullptr;
      if ( pmanWfs ) {
        //IndexVector wfPlotChans = nnbr==0 ? {icha} : {icha-1, icha, icha+1};
        //for ( Index ichaPlot : wfPlotChans ) {
        //for ( Index ichaPlot : (nnbr==0 ? {icha} : {icha-1, icha, icha+1}) ) {
        cout << myname << "Creating waveform plot." << pmanWfs.get() << endl;
        string hopt = "h";
        if ( ipla ) hopt += " same";
        for ( Index iplt=0; iplt<nWfPlot; ++iplt ) {
          Index ichaPlot = icha - nWfPlot/2 + iplt;
          TPadManipulator* pmanPlot = pmanWfs->man(iplt);
          const AdcChannelData& acdPlot = acm[ichaPlot];
          DataMap dm = pwf->view(acdPlot);
          string hnam = "prepared";
          TH1* ph = dm.getHist(hnam);
          if ( ph == nullptr ) {
            cout << myname << "Waveform hist not found: " << hnam << endl;
            return 30;
          }
          ph->Print();
          ph->SetLineColor(mycols[ipla]);
          ph->SetLineWidth(mywids[ipla]);
          ph->SetLineStyle(mystys[ipla]);
          string slab;
          float thr = 3*noise;
          //thr = 0.5;
          float sampleNoise = wfPlotNoiseLevels[iplt];
          if ( nevt ) {
            float noiseThreshold = 3.0*sampleNoise;
            if ( true ) {
              getResponseLabelNoise(ph, slab, sampleNoise, thrFac, 2);
            } else {
              noiseThreshold = 0.0;
              getResponseLabelThreshold(ph, slab, 3*sampleNoise, 2);
            }
            slab = mylabs[ipla] + ": " + slab;
            slabsWf[iplt].push_back(slab);
            cout << myname << slab << endl;
          }
          pmanPlot->add(ph, hopt);
        }  // End loop over plotted channels
      }  // End fill of waveform plots.
    }  // End loop over planes.
    // Decorate and scale the waveform plots.
    if ( pmanWfs ) {
      ostringstream sslab;
      sslab << "Threshold is ";
      if ( thrFac == 0.0 ) sslab << "0";
      else if ( thrFac == 1.0 ) sslab << "N";
      else sslab << thrFac << "#timesN";
      for ( Index iplt=0; iplt<nWfPlot; ++iplt) {
        slabsWf[iplt].push_back(sslab.str());
        TPadManipulator& man = *pmanWfs->man(iplt);
        // Decorate and print the waveform plot for this event.
        man.addAxis();
        man.setRangeX(wfxmin, wfxmax);
        man.setRangeY(wfymin, wfymax);
        float xlab = 0.12;
        float ylab = 0.86;
        float tsiz = 0.035;
        for ( string slab : slabsWf[iplt] ) {
          TLatex* ptxt = new TLatex(xlab, ylab, slab.c_str());
          ptxt->SetTextFont(42);
          ptxt->SetTextSize(tsiz);
          ptxt->SetNDC();
          man.add(ptxt);
          ylab -= 1.3*tsiz;
        }
        string sttl = responseTitle;
        sttl += " " + smipfac + " MIP " + sdepRange + " for #theta_{xz}=" + sthtxz + "#circ with " + slabsmr;
        string ssuf;
        ssuf += "N_{s} =" + to_string(nsam) + ", ";
        ssuf += "N_{nb} =" + to_string(nnbr);
        if ( use2dConvolution ) ssuf += string(ssuf.size() ? ", " : "") + "N_{bin} =" + to_string(nbinPerWire);
        if ( ssuf.size() ) sttl += " (" + ssuf + ")";
        if ( sttlPrefix.size() ) sttl = sttlPrefix + sttl;
        man.setTitle(sttl);
        man.setTitle(sglobalLabel);
        man.addHorizontalLine(0.0);
        man.setLabel("Run " + to_string(irun) + " event " + to_string(ievt) + chanLabel[iplt]);
      }
      // Print the waveforms for this event.
      string fnam = "wf";
      fnam += (doDeconvolution ? deconName : doConvolution ? "conv" : "noconv");
      fnam += "-sigma" + ssigl;
      fnam += "-" + strkLab;
      if ( showNsam ) {
        string snsam = to_string(nsam);
        while ( snsam.size() < 5 ) snsam = "0" + snsam;
        fnam += "-nsam" + snsam;
      }
      fnam += "-nnbr" + to_string(nnbr);
      fnam += "-nbin" + to_string(nbinPerWire);
      string sevt = to_string(ievt);
      while ( sevt.size() < 3 ) sevt = "0" + sevt;
      fnam += "-evt" + sevt;
      fnam += ".{png,tpad}";
      cout << myname << "Printing " << fnam << endl;
      pmanWfs->print(fnam);
    }
    // Notify tools of end of event.
    pplotDftBefore->endEvent(evt);
    pplotDftAfter->endEvent(evt);
  }
  cout << myname << "Closing tools." << endl;
  pplotDftBefore.reset();
  pplotDftAfter.reset();
  dmRoiview.setInt("dbg", 1);
  proiview->close(&dmRoiview);
  proiview.reset();
  cout << myname << sglobalLabel << endl;
  cout << myname << "Exiting." << endl;
  return 0;
}
