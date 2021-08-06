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
//       CVALS = NNBR:NBIN:NSIG:NOISE
//         NNBR = # neighbor cells to include (-1 for 1D simulation) [0]
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
//          deco1d:NBIN = 2D deconvolution with NBIN bins/wire
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

int run(float thtxzdeg, string sresp, string sdeco, unsigned int irun =0, unsigned int nevt =1) {
  string myname = "run: ";
  using Index = unsigned int;

  // Define some parameters.
  float dzcell = 4.79;         // Cell size in mm
  float dxdt = 0.8;            // Drift speed in mm/tick
  float dqds = 6.0;            // Rate of charge deposit in ke/mm
  Index nsam =     0;          // Number of time samples
  bool showNsam = nsam > 0;    // Flag to include nsam in file name.
  bool doRebaseline = true;    // Rebaselining flag
  bool mostDistantResponse = false;   // Only for debugging

  // Get the track slope and the nominal signal widh (ticks).
  // Use the latter to set the signal position isig to center the v-plane signal.
  float thtxz = thtxzdeg*acos(-1.0)/180.0;
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
  float wfyfac = 1.0;
  float wfymin = -3.0*wfyfac;  // Y-range for waveform plots
  float wfymax =  6.0*wfyfac;
  if ( false ) {
    wfymax = 50;
    wfymin = -wfymax;
  }
  bool doWideXRange = true;
  if ( doWideXRange ) {
    if ( nsam > 400 ) {
      wfxmin = isig - 400;
      wfxmax = isig + 400;
    } else {
      wfxmin = 0;
      wfxmax = nsam;
    }
  }

  // Decode the response nnbr-nbinPerWire-ssigl-enoise
  Index nnbr = 0;
  Index nbinPerWire = 1;
  string ssigl = "0";
  Index enoise = 0;
  bool doSig = true;
  bool doConvolution = nbinPerWire > 0;
  bool use2dConvolution = true;
  Index badUInt = 9999;
  if ( sresp.size() ) {
    vector<string> rsubs = StringManipulator(sresp).split(":", true);
    Index nsub = rsubs.size();
    string sconv = rsubs[0];
    if ( sconv == "nosig" ) {
      doSig = false;
      doConvolution = false;
      use2dConvolution = false;
    } else if ( sconv == "noconv" ) {
      doConvolution = false;
    } else if ( sconv == "conv1d" ) {
      use2dConvolution = false;
    } else if ( sconv == "conv2d" ) {
      use2dConvolution = true;
    } else {
      sconv = "UnknownConversion";
    }
    if ( nsub > 1 ) {
      StringManipulator smval(rsubs[1]);
      if ( smval.size() ) nnbr = smval.toUnsignedInt(badUInt);
    }
    if ( nsub > 2 ) {
      StringManipulator smval(rsubs[2]);
      if ( smval.size() ) nbinPerWire = smval.toUnsignedInt(badUInt);
    }
    if ( nsub > 3 ) ssigl = rsubs[3];
    if ( nsub > 4 ) {
      StringManipulator smval(rsubs[4]);
      if ( smval.size() ) enoise = smval.toUnsignedInt(badUInt);
    }
    if ( sconv == "UnknownConversion" || nnbr == badUInt || nbinPerWire == badUInt || enoise == badUInt ) {
      cout << myname << "ERROR: Invalid response configuration: " << sresp << " --> "
           << sconv << ":" << nnbr << ":" << nbinPerWire << ":" << ssigl << ":" << enoise << endl;
      return 1;
    }
  } else {
    cout << myname << "ERROR: Response configuration must be provided." << endl;
    return 1;
  }

  // Decode the deconvolution specifier.
  vector<string> dsubs = StringManipulator(sdeco).split(":", true);
  string dnam = dsubs[0];
  bool do1dDeconvolution = false;
  bool do2dDeconvolution = false;
  vector<string> deconToolNames;
  string deconTitle;
  int deconError = 0;
  if ( dnam == "nodeco" ) {
    cout << myname << "No deconvolution." << endl;
  } else if ( dnam == "deco1d" ) {
    do1dDeconvolution = true;
    string dopt = dsubs.size() == 2 ? dsubs[1] : "BAAAD";
    if ( dopt == "dft" ) {
      cout << myname << "1D DFT deconvolution." << endl;
      deconToolNames.push_back("decon");
      deconTitle = "1D DFT deconvoluted reponse to";
    } else if ( dopt == "dm" ) {
      cout << myname << "1D direct deconvolution." << endl;
      deconToolNames.push_back("dmdecon");
      deconTitle = "1D Direct matrix deconvoluted response to";
    } else if ( dopt == "fm" ) {
      cout << myname << "1D filtered matrix deconvolution." << endl;
      deconToolNames.push_back("fmdecon");
      deconTitle = "1D Filtered matrix deconvoluted reponse to";
    } else if ( dopt == "cm" ) {
      cout << myname << "1D chi-square matrix deconvolution." << endl;
      deconToolNames.push_back("cmdecon");
      deconTitle = "1D Chi-square matrix deconvoluted reponse to";
    } else deconError = 2;
  } else if ( dnam == "deco2d" ) {
    deconTitle = "2D deconvoluted reponse to";
    if ( dsubs.size() != 2 ) deconError = 3;
    cout << myname << "ERROR: 2D deconvolution is not yet supported." << endl;
    return 1;
  } else deconError = 1;
  if ( deconError ) {
    cout << "Invalid deconvolution option: " << sdeco << endl;
    return 1;
  }
  bool doDeconvolution = do1dDeconvolution || do2dDeconvolution;
  string deconName;
  for ( string dsub : dsubs ) deconName += dsub;

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
  cout << myname << "Signal tick: " << isig << endl;

  // Fetch noise scaled to 1.0 for this run.
  CeNoise cen(irun);

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
      if ( fetchTool(tnam, pcon2ds[nbin][sview], myname) ) return 2;
    }
  }
  vector<TpcDataTool*> deconTools(deconToolNames.size(), nullptr);
  cout << myname << "Fetching deconvolution tool(s)." << endl;
  auto idt = deconTools.begin();
  for ( string deconToolName : deconToolNames ) {
    if ( fetchTool(deconToolName, *(idt++), myname) ) return 3;
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
      return 2;
    }
    slabsmr = "#sigma_{smear} = " + ssigl + " ticks";
  } else if ( ssigl != "0" ) {
    cout << myname << "ERROR: Invalid ssigl: " << ssigl << endl;
    return 2;
  }
  cout << myname << "Fetching DFT before plotter." << endl;
  auto pplotDftBefore = ptm->getPrivate<TpcDataTool>("plotDftBefore");
  if ( ! pplotDftBefore ) return 3;
  cout << myname << "Fetching DFT after plotter." << endl;
  auto pplotDftAfter  = ptm->getPrivate<TpcDataTool>("plotDftAfter");
  if ( ! pplotDftAfter ) return 3;
  cout << myname << "Fetching FFT calculator." << endl;
  TpcDataTool* pfft = ptm->getShared<TpcDataTool>("adcFFT");
  if ( pfft == nullptr ) return 3;
  cout << myname << "Fetching signal finder." << endl;
  TpcDataTool* psigfind = doDeconvolution ? ptm->getShared<TpcDataTool>("deconSignalFinder")
                                             : ptm->getShared<TpcDataTool>("nodeconSignalFinder");
  if ( psigfind == nullptr ) return 3;
  cout << myname << "Fetching sample RMS calculator." << endl;
  TpcDataTool* psamrms = ptm->getShared<TpcDataTool>("adcChannelSamplelRmsPlotter");
  if ( psamrms == nullptr ) return 3;
  cout << myname << "Fetching signal RMS calculator." << endl;
  TpcDataTool* psigrms = ptm->getShared<TpcDataTool>("adcChannelSignalRmsPlotter");
  if ( psigrms == nullptr ) return 3;
  cout << myname << "Fetching rebaseline tool." << endl;
  TpcDataTool* prbl = ptm->getShared<TpcDataTool>("rebaseline");
  //TpcDataTool* prbl = ptm->getShared<TpcDataTool>("adcPedestalFit");
  if ( prbl == nullptr ) return 3;
  cout << myname << "Fetching ROI viewer." << endl;
  string rvname = "roiViewer";
  if ( thtxzdeg > 80 ) rvname += "Wide2";
  else if ( thtxzdeg > 70 ) rvname += "Wide";
  auto proiview = ptm->getPrivate<TpcDataTool>(rvname);
  if ( proiview == nullptr ) return 3;
  // End fetching tools.

  // Create charge deposition.
  cout << myname << line << endl;
  cout << myname << "Depositing charge." << endl;
  cout << myname << line << endl;
  // Charge is deposited in time-wire bins.
  // We integrate across each 0.5 us time bin starting at the left edge of bin 0.
  string sttlPrefix;   // Used to build title for waveform plots.
  map<int, vector<float>> binsamples;     // Q[ibin][isam]
  map<int, vector<float>> celsamples;     // Q[icel][isam]
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
      cout << myname << "Wire bin " << iwir << endl;
      vector<float>& wsamples = kwir.second;
      // Handle track parallel to APA
      if ( dxdz == 0.0 ) {
        float qdep = dqdw*wbWid;
        wsamples[isig] = qdep;
        cout << myname << "  Q[" << iwir << "][" << isig << "] = " << qdep << endl;
        continue;
      }
      float xwir = wbOff + wbWid*iwir;
      float tmin = ttOff + dtdw*(iwir)*wbWid;
      float tmax = ttOff + dtdw*(iwir+1)*wbWid;
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
          cout << myname << "  Q[" << iwir << "][" << it << "] = " << qdep << endl;
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
  cout << myname << "Bin sample counts, charge: " << endl;
  for ( auto& kwir : binsamples ) {
    int iwir = kwir.first;
    float qsum = 0.0;
    for ( float q : binsamples[iwir] ) qsum += q;
    cout << myname << setw(6) << iwir << ": " << binsamples[iwir].size() << ", " << qsum << endl;
  }
  cout << myname << "Cell sample counts: " << endl;
  for ( auto& kwir : celsamples ) {
    int iwir = kwir.first;
    float qsum = 0.0;
    for ( float q : celsamples[iwir] ) qsum += q;
    cout << myname << setw(6) << iwir << ": " << celsamples[iwir].size() << ", " << qsum << endl;
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

  unsigned int ncon = doConvolution ? 3 : 1;
  vector<TH1*> myhsts(ncon, nullptr);
  Index maxEvtPlot = 20;    // # events for which waveform plots are made
  unique_ptr<TPadManipulator> pman;
  // Loop over events.
  for ( Index ievt=1; ievt<nevt+1; ++ievt ) {
    cout << myname << line << endl;
    cout << myname << "######## Event " << ievt << endl;
    cout << myname << line << endl;
    string hopt = "h";
    vector<string> slabs;
    // Create event info.
    DuneEventInfo evt;
    evt.run = irun;
    evt.event = ievt;
    // Display noise state.
    cout << myname << "Noise counts: " << cen.count() << ", " << cen.next()
         << ".  Last = " << cen.last() << endl;
    // Loop over the three views unless we are only showing the deposition.
    if ( ievt < maxEvtPlot ) pman.reset(new TPadManipulator);
    else pman.reset(nullptr);
    bool doWfPlot = false;
    for ( unsigned int icon=0; icon<ncon; ++icon ) {
      Index icha = chans[icon] + thtxzdeg;
      string sview = views[icon];
      // Create ADC channel data with the energy deposition
      AdcChannelDataMap acm;
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
        //jacd.run = irun;
        //jacd.event = ievt;
        jacd.setEventInfo(irun, ievt);
        jacd.sampleUnit = "ke/tick";
        //jacd.channel = jcha;
        jacd.setChannelInfo(jcha, ifmb, ifch, chanStat);
        jacd.samples.resize(nsam);
        //jacd.channelStatus = 0;
      }
      AdcChannelData& acd = acm[icha];
      if ( ! doConvolution ) {
        // Fill samples from cells.
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
      } else if ( use2dConvolution ) {
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
              return 4;
            }
            cout << myname << "  Filling chan/bin " << jcha << "/" << jbin << " from edep bin " << ibin << endl;
            vector<float>& dest = jacd.binSamples[jbin];
            dest.resize(nsam);
            jacd.binSamples[jbin].resize(nsam);
            for ( Index isam=0; isam<nsam; ++isam ) dest[isam] = binsamples[ibin][isam];
          }  // End loop over channel bins
        }
        if ( pcon2ds.count(nbinPerWire) == 0 ) {
          cout << myname << "Tool for 2d convolution with " << nbinPerWire << " bins not found." << endl;
          return 1;
        }
        TpcDataTool* pcon2d = pcon2ds[nbinPerWire][sview];
        pcon2d->updateMap(acm).print();
      } else {
        cout << myname << "Convoluting with neighbors" << endl;
        for ( int icel : apaConCells ) {
          cout << myname << "...cell " << icel << endl;
          TpcDataTool* pncon = abs(icel) == 0 ? pconvo  :
                                  abs(icel) == 1 ? pconvo1 :
                                  abs(icel) == 2 ? pconvo2 :
                                  abs(icel) == 3 ? pconvo3 :
                                  abs(icel) == 4 ? pconvo4 : nullptr;
          if ( pncon == nullptr ) {
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
          pncon->update(acdtmp);
          for ( Index isam=0; isam<acdtmp.samples.size(); ++isam ) acd.samples[isam] += acdtmp.samples[isam];
        }  // End loop over neighbors
      }
      if ( noise > 0.0 ) {
        cout << myname << "Adding noise " << icha << endl;
        for ( float& sam : acd.samples ) sam += noise*cen.get();
        acd.dftmags.clear();
        acd.dftphases.clear();
      }
      if ( true ) {
        pfft->update(acd);
        pplotDftBefore->viewMap(acm);
        // Deconvolution.
        if ( doDeconvolution ) {
          cout << myname << "Deconvoluting channel " << icha << endl;
          for ( TpcDataTool* pdecon : deconTools ) pdecon->update(acd).print();
          if ( doRebaseline ) {
            cout << myname << "Rebaselining channel " << icha << endl;
            prbl->update(acd).print();
          }
          pplotDftAfter->viewMap(acm);
        }
      }
      // Find ROI and update ROI histograms.
      cout << myname << "Finding signal" << endl;
      psigfind->update(acd).print();
      cout << myname << "Updating ROI histograms" << endl;
      proiview->view(acd).print();
      cout << myname << "Evaluating background RMS" << endl;
      // Invert signal finder to get background nose level.
      for ( Index isam=0; isam<acd.signal.size(); ++isam ) acd.signal[isam] = acd.signal[isam] ? false : true;
      DataMap dsamrms = psamrms->view(acd);
      cout << myname << "Sample RMS: " << dsamrms.getFloat("metricValue") << endl;
      DataMap dsigrms = psigrms->view(acd);
      cout << myname << "Signal RMS: " << dsigrms.getFloat("metricValue") << endl;
      //dsamrms.print();
      // Fetch the waveform and add it to the waveform plot.
      TH1* ph = nullptr;
      if ( pman ) {
        TPadManipulator& man = *pman.get();
        DataMap dm = pwf->view(acd);
        string hnam = "prepared";
        TH1* ph = dm.getHist(hnam);
        if ( ph == nullptr ) {
          cout << myname << "Waveform hist not found: " << hnam << endl;
          return 1;
        }
        ph->Print();
        ph->SetLineColor(mycols[icon]);
        ph->SetLineWidth(mywids[icon]);
        ph->SetLineStyle(mystys[icon]);
        string slab;
        float thr = 3*noise;
        //thr = 0.5;
        if ( nevt ) {
          getResponseLabel(ph, slab, thr, 1);
          slab = mylabs[icon] + ": " + slab;
          slabs.push_back(slab);
          cout << myname << slab << endl;
        }
        man.add(ph, hopt);
        hopt = "h same";
      }
    }  // End loop over views.
    // Notify tools of end of event.
    pplotDftBefore->endEvent(evt);
    pplotDftAfter->endEvent(evt);
    if ( pman ) {
      TPadManipulator& man = *pman.get();
      // Decorate and print the waveform plot for this event.
      man.addAxis();
      man.setRangeX(wfxmin, wfxmax);
      man.setRangeY(wfymin, wfymax);
      float xlab = 0.12;
      float ylab = 0.86;
      float tsiz = 0.035;
      for ( string slab : slabs ) {
        TLatex* ptxt = new TLatex(xlab, ylab, slab.c_str());
        ptxt->SetTextFont(42);
        ptxt->SetTextSize(tsiz);
        ptxt->SetNDC();
        man.add(ptxt);
        ylab -= 1.3*tsiz;
      }
      ostringstream ssout;
      ssout << std::fixed << std::setprecision(0) << thtxzdeg;
      string stht = ssout.str();
      string sttl = responseTitle;
      sttl += " 1 MIP for #theta_{xz}=" + stht + "#circ with " + slabsmr;
      string ssuf;
      ssuf += "N_{nb} =" + to_string(nnbr);
      if ( use2dConvolution ) ssuf += string(ssuf.size() ? ", " : "") + "N_{bin} =" + to_string(nbinPerWire);
      if ( ssuf.size() ) sttl += " (" + ssuf + ")";
      if ( sttlPrefix.size() ) sttl = sttlPrefix + sttl;
      man.setTitle(sttl);
      man.addHorizontalLine(0.0);
      man.setLabel("Run " + to_string(irun) + " event " + to_string(ievt));
      string fnam = "wf";
      fnam += (doDeconvolution ? deconName : doConvolution ? "conv" : "noconv");
      fnam += "-sigma" + ssigl;
      fnam += "-thtxz" + stht;
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
      man.print(fnam);
    }
  }
  cout << myname << "Closing tools." << endl;
  pplotDftBefore.reset();
  pplotDftAfter.reset();
  proiview.reset();
  cout << myname << "Exiting." << endl;
  return 0;
}
