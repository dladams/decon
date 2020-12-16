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
// There is also the option (ssig) to add diffusion smearing of 1 or 2 ticks
// roughly corresponding respectively to 1 or 3.6 m drift.
// Noise is added using CeNoise which is random in each tick followed by
// smearing with the CE response function and then scaled to have RMS 1.0.
// This is the scaled by the noise level (noise) in ke.
//
// thtxzdeg: horizontal angle w.r.t. the APA plane
// ssig:
//      "" - no smearing
//      "1" - smear RMS = 1 ==> 1 m
//      "2" - smear RMS = 2 ==> 3.6 m
//      "nosig" - No signal (noise only)
// iproc specifies the extent of processing:
// iproc = 10*iprocRec + iprocSim
//   iprocSim: 0 - energy deposit only
//             1 - 1D response
//             2 - Neighbor cell 2D response
//             3 - Binned 2D response
//   iprocRec: 0 - no deconvolution
//             1 - 1D DFT deconvolution
//             2 - 1D directr matrix deconvolution
//             3 - 1D direct matrix deconvolution
//             4 - 1D filtered matrix deconvolution
// noise is the noise level in ke
// run is the run number and is used to seed the noise
// nevt is the number of events to process (1 track with 3 views for each event)
//
// Products of the script include:
// Waveform plots - One for each event overlaying the three views.
// Power plots - One for each view, before and after reonstruction.

int fetchTool(string tnam, AdcChannelTool*& ptoo, string prefix) {
  ptoo = ptm->getShared<AdcChannelTool>(tnam);
  if ( ptoo == nullptr ) {
    cout << prefix << "ERROR: Unable to find tool " << tnam << endl;
    return 1;
  }
  cout << prefix << "Found tool " << tnam << endl;
  return 0;
}

int run(float thtxzdeg, string ssig, int iproc =1, float noise = 0.0, unsigned int irun =0, unsigned int nevt =1) {
  string myname = "run: ";
  using Index = unsigned int;
  Index nsam =   300;     // Number of time samples
  Index nchDep = 2;       // Number of neighboring cells with signals
  int ncel = int(nchDep);
  Index nbinPerWire = 2;  // Number of energy deposit bins per wire.
  vector<int> apaConCells; // Cells contributing to central signal.
  bool mostDistantResponse =  true;   // Only for debugging
  if ( mostDistantResponse ) {
    apaConCells.push_back(-int(nchDep));
    apaConCells.push_back( int(nchDep));
  } else {
    for ( int icel=-ncel; icel<=ncel; ++icel ) apaConCells.push_back(icel);
  }
  Index isig = nsam/2;    // Position of the signal.
  float wfxmin = isig - 45;
  float wfxmax = isig + 45;
  float wfyfac = 0.1;
  float wfymin = -3.0*wfyfac;
  float wfymax =  6.0*wfyfac;
  if ( false ) {
    wfymax = 50;
    wfymin = -wfymax;
  }

  bool doWideXRange = true;
  bool doOldXRange = false;
  if ( doOldXRange ) {
    isig -= 20;
    wfxmin -= 5;
    wfxmax -= 5;
    if ( thtxzdeg > 80 ) wfxmax += 40;
    doWideXRange = false;
  }
  bool doRebaseline = true;
  if ( doWideXRange ) {
    if ( nsam > 400 ) {
      wfxmin = isig - 400;
      wfxmax = isig + 400;
    } else {
      wfxmin = 0;
      wfxmax = nsam;
    }
  }
    
  // Extract the simulation and deconvolution options.
  // The values will be checked when the correponding tools are fetched.
  int iprocSim = iproc%10;
  int iprocRec = iproc/10;

  // Fetch noise scaled to 1.0 for this run.
  CeNoise cen(irun);

  // Fetch the tools.
  // Tools that use Root objects in their dtors are created as private and deleted
  // here to ensure those objects are still present.
  cout << myname << "Fetching waveform plotter." << endl;
  AdcChannelTool* pwf = ptm->getShared<AdcChannelTool>("wf");
  cout << myname << "Fetching 1D reponse convolution tools." << endl;
  AdcChannelTool* pconvo = ptm->getShared<AdcChannelTool>("convo");
  AdcChannelTool* pconvo1 = ptm->getShared<AdcChannelTool>("convo1");
  AdcChannelTool* pconvo2 = ptm->getShared<AdcChannelTool>("convo2");
  AdcChannelTool* pconvo3 = ptm->getShared<AdcChannelTool>("convo3");
  AdcChannelTool* pconvo4 = ptm->getShared<AdcChannelTool>("convo4");
  cout << myname << "Fetching 2D reponse convolution tools." << endl;
  AdcChannelTool *pconv2dx, *pconv2du, *pconv2dv;
  if ( fetchTool("conv2dx", pconv2dx, myname) ) return 2;
  if ( fetchTool("conv2du", pconv2du, myname) ) return 2;
  if ( fetchTool("conv2dv", pconv2dv, myname) ) return 2;
  AdcChannelTool* pdecon = nullptr;
  string responseTitle;
  if ( iprocSim == 0 ) {
    cout << myname << "No convolution." << endl;
    responseTitle = "Response from";
  } else if ( iprocSim == 1 ) {
    cout << myname << "1D convolution." << endl;
    responseTitle = "1D convolution of";
  } else if ( iprocSim == 2 ) {
    cout << myname << "Neighbor convolution." << endl;
    responseTitle = "Neighbor convolution of";
  } else if ( iprocSim == 3 ) {
    cout << myname << "Binned 2D convolution." << endl;
    responseTitle = "Binned 2D convolution of";
  } else {
    cout << "Invalid simulation option: " << iprocSim << endl;
  }
  cout << myname << "Fetching deconvolution tool(s)." << endl;
  string deconName;
  if ( iprocRec == 0 ) {
    cout << myname << "No deconvolution." << endl;
  } else if ( iprocRec == 1 ) {
    cout << myname << "1D DFT deconvolution." << endl;
    deconName = "decon";
    responseTitle = "DFT deconvoluted reponse to";
  } else if ( iprocRec == 2 ) {
    cout << myname << "1D direct deconvolution." << endl;
    deconName = "dmdecon";
    responseTitle = "Direct matrix deconvoluted reponse to";
  } else if ( iprocRec == 3 ) {
    cout << myname << "1D filtered matrix deconvolution." << endl;
    deconName = "fmdecon";
    responseTitle = "Filtered matrix deconvoluted reponse to";
  } else if ( iprocRec == 4 ) {
    cout << myname << "1D chi-square matrix deconvolution." << endl;
    deconName = "cmdecon";
    responseTitle = "Chi-square matrix deconvoluted reponse to";
  } else {
    cout << "Invalid deconvolution option: " << iprocRec << endl;
  }
  if ( deconName.size() && fetchTool(deconName, pdecon, myname) ) return 3;
  AdcChannelTool* pconvg = nullptr;
  string slabsmr = "no smearing";
  bool doSig = ssig != "nosig";
  if ( ssig == "1" || ssig == "2" ) {
    string sconvg = "convg" + ssig;
    cout << myname << "Fetching filter smearing tool " << sconvg << endl;
    pconvg = ptm->getShared<AdcChannelTool>(sconvg);
    if ( pconvg == nullptr ) {
      cout << myname << "ERROR: Unable to find tool " << sconvg << endl;
      return 2;
    }
    slabsmr = "#sigma_{smear} = " + ssig + " ticks";
  } else if ( ssig == "-1" ) {
    cout << myname << "Skipping filter smearing tool." << endl;
    doSig = false;
  } else if ( ssig != "0" ) {
    cout << myname << "Invalid ssig: " << ssig << endl;
    return 2;
  }
  cout << myname << "Fetching DFT before plotter." << endl;
  auto pplotDftBefore = ptm->getPrivate<AdcChannelTool>("plotDftBefore");
  if ( ! pplotDftBefore ) return 3;
  cout << myname << "Fetching DFT after plotter." << endl;
  auto pplotDftAfter  = ptm->getPrivate<AdcChannelTool>("plotDftAfter");
  if ( ! pplotDftAfter ) return 3;
  cout << myname << "Fetching FFT calculator." << endl;
  AdcChannelTool* pfft = ptm->getShared<AdcChannelTool>("adcFFT");
  if ( pfft == nullptr ) return 3;
  cout << myname << "Fetching signal finder." << endl;
  AdcChannelTool* psigfind = ptm->getShared<AdcChannelTool>("deconSignalFinder");
  if ( psigfind == nullptr ) return 3;
  cout << myname << "Fetching sample RMS calculator." << endl;
  AdcChannelTool* psamrms = ptm->getShared<AdcChannelTool>("adcChannelSamplelRmsPlotter");
  if ( psamrms == nullptr ) return 3;
  cout << myname << "Fetching signal RMS calculator." << endl;
  AdcChannelTool* psigrms = ptm->getShared<AdcChannelTool>("adcChannelSignalRmsPlotter");
  if ( psigrms == nullptr ) return 3;
  cout << myname << "Fetching rebaseline tool." << endl;
  AdcChannelTool* prbl = ptm->getShared<AdcChannelTool>("rebaseline");
  //AdcChannelTool* prbl = ptm->getShared<AdcChannelTool>("adcPedestalFit");
  if ( prbl == nullptr ) return 3;
  cout << myname << "Fetching ROI viewer." << endl;
  string rvname = "roiViewer";
  if ( thtxzdeg > 80 ) rvname += "Wide2";
  else if ( thtxzdeg > 70 ) rvname += "Wide";
  auto proiview = ptm->getPrivate<AdcChannelTool>(rvname);
  if ( proiview == nullptr ) return 3;

  // Create charge deposition.
  // Charge is deposited in time-wire bins.
  // We integrate across each 0.5 us time bin starting at the left edge of bin 0.
  string sttlPrefix;
  map<int, vector<float>> binsamples;     // Q[ibin][isam]
  map<int, vector<float>> celsamples;     // Q[icel][isam]
  for ( int icel=-ncel; icel<=ncel; ++icel ) {
    celsamples[icel].resize(nsam, 0.0);
  }
  if ( doSig ) {
    float dzcell = 4.79;    // mm
    float dxdt = 0.8;  // mm/tick
    float dqds = 6.0;  // ke/mm
    float thtxz = thtxzdeg*acos(-1.0)/180.0;
    float dxdz = tan(thtxz);
#if 0
    if ( nbinPerWire != 1 ) {
      cout << myname << "Error: Must have 1 bin/wire for old simulation." << endl;
    }
    sttlPrefix = "Old ";
    float dxcell = dxdz*dzcell;
    float dtcell = dxcell/dxdt;
    float dscell = sqrt(dzcell*dzcell + dxcell*dxcell);
    float dqcell = dqds*dscell;
    float dt = dtcell;
    float dq = dqcell;
    cout << myname << "  dzcell = " << dzcell << " mm" << endl;
    cout << myname << "  dscell = " << dscell << " mm" << endl;
    cout << myname << "  dtcell = " << dtcell << " tick" << endl;
    cout << myname << "  dqcell = " << dqcell << " ke" << endl;
    unsigned int itick = isig;
    // Loop over wire bins.
    for ( auto& kwir : celsamples ) {
      int iwir = kwir.first;
      vector<float>& csamples = kwir.second;
      // Handle track parallel to APA
      if ( dxdz == 0.0 ) {
        csamples[isig] = dqcell;
        cout << myname << "  Q[" << iwir << "][" << isig << "] = " << dqcell << endl;
        continue;
      }
      float tmin = isig + iwir*dtcell;
      float tmax = isig + (iwir+1)*dtcell;
      int it1 = tmin - 1.0;
      int it2 = tmax + 1.0;
      float dqdt = dqcell/dtcell;
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
          csamples[it] = qdep;
          cout << myname << "  Q[" << iwir << "][" << it << "] = " << qdep << endl;
        }
      }
    }  // End loop over wire bins.
#else
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
    // Temporary: fill cells from bins.
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
#endif
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

  unsigned int ncon = 3;
  if ( iprocSim == 0 ) ncon = 1;
  vector<TH1*> myhsts(ncon, nullptr);
  Index maxEvtPlot = 20;    // # events for which waveform plots are made
  unique_ptr<TPadManipulator> pman;
  // Loop over events.
  for ( Index ievt=1; ievt<nevt+1; ++ievt ) {
    cout << myname << "################################# Event " << ievt << endl;
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
      // Create ADC channel data with the energy deposition
      AdcChannelDataMap acm;
      AdcChannelData& acd = acm[icha];
      acd.run = irun;
      acd.event = ievt;
      acd.sampleUnit = "ke/tick";
      acd.channel = icha;
      acd.samples.resize(nsam);
      acd.channelStatus = 0;
      if ( iprocSim == 0 || iprocSim == 1 ) {
        // Fill samples from cells.
        for ( Index isam=0; isam<celsamples[0].size(); ++isam ) {
          if ( isam >= nsam ) break;
          acd.samples[isam] = celsamples[0][isam];
        }
        if ( pconvg != nullptr ) {
          cout << myname << "Adding longitudinal diffusion to single channel " << icha << endl;
          pconvg->update(acd).print();
        }
      }
      if ( iprocSim > 0 ) {
        if ( iprocSim == 1 ) {
          cout << myname << "Convoluting with detector response " << icha << endl;
          pconvo->update(acd).print();
        } else if ( iprocSim == 2 ) {
          cout << myname << "Adding APA contributions from cells" << endl;
          for ( int icel : apaConCells ) {
            cout << myname << "...cell " << icel << endl;
            AdcChannelTool* pncon = abs(icel) == 0 ? pconvo  :
                                    abs(icel) == 1 ? pconvo1 :
                                    abs(icel) == 2 ? pconvo2 :
                                    abs(icel) == 3 ? pconvo3 :
                                    abs(icel) == 4 ? pconvo4 : nullptr;
            if ( pncon == nullptr ) {
              cout << myname << "WARNING: Convolution tool not found for cell " << icel << endl;
              continue;
            }
            AdcChannelData acdtmp;
            acdtmp.channel = icha + icel;
            acdtmp.samples = celsamples[icel];
            if ( pconvg != nullptr ) {
              cout << myname << "Adding longitudinal diffusion to channel " << acdtmp.channel << endl;
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
        pfft->update(acd);
        pplotDftBefore->viewMap(acm);
        // Deconvolution.
        if ( iprocRec > 0 ) {
          cout << myname << "Deconvoluting channel " << icha << endl;
          pdecon->update(acd).print();
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
        float thr = 0.0;
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
      ssout << std::fixed << std::setprecision(1) << thtxzdeg;
      string stht = ssout.str();
      string sttl = responseTitle;
      sttl += " 1 MIP for #theta_{xz}=" + stht + "#circ with " + slabsmr;
      if ( iprocSim==2 && nchDep ) sttl += " (N_{nb} =" + to_string(nchDep) + ")";
      if ( sttlPrefix.size() ) sttl = sttlPrefix + sttl;
      man.setTitle(sttl);
      man.addHorizontalLine(0.0);
      man.setLabel("Run " + to_string(irun) + " event " + to_string(ievt));
      string fnam = "wf";
      fnam += (iprocRec>0 ? deconName : iprocSim ? "conv" : "noconv");
      fnam += "-sigma" + ssig;
      fnam += "-thtxz" + stht;
      string snsam = to_string(nsam);
      while ( snsam.size() < 5 ) snsam = "0" + snsam;
      fnam += "nsam" + snsam;
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
