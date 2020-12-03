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
//    0 - Energy deposit and diffusion smearing
//    1 - Convoluted with detector response
//    2 - Deconvoluted
//    3 - Direct matrix deconvoluted
//    4 - Filtered matrix deconvoluted
//    5 - Chi-square matrix deconvoluted
// noise is the noise level in ke
// run is the run number and is used to seed the noise
// nevt is the number of events to process (1 track with 3 views for each event)
//
// Products of the script include:
// Waveform plots - One for each event overlaying the three views.
// Power plots - One for each view, before and after reonstruction.

int run(float thtxzdeg, string ssig, int iproc =1, float noise = 0.0, unsigned int irun =0, unsigned int nevt =1) {
  string myname = "run: ";
  using Index = unsigned int;
  Index nsam =    70;     // Number of time samples
  Index nchDep = 1;       // Number of neighboring cells with signals
  Index nbinPerWire = 2;  // Number of energy deposit bins per wire.
  vector<int> apaConCells = {-1, 1};      // Cells contributing to central signal.
  Index isig = nsam/2;    // Position of the signal.
  float wfxmin = isig - 45;
  float wfxmax = isig + 45;
  float wfymin = -2.0;
  float wfymax = 7.0;
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
    

  // Check arguments.
  if ( iproc > 5 ) {
    cout << myname << "Invalid conv option: " << iproc << endl;
    return 1;
  }

  // Fetch noise scaled to 1.0 for this run.
  CeNoise cen(irun);

  // Fetch the tools.
  // Tools that use Root objects in their dtors are created as private and deleted
  // here to ensure those objects are still present.
  cout << myname << "Fetching waveform plotter." << endl;
  AdcChannelTool* pwf = ptm->getShared<AdcChannelTool>("wf");
  cout << myname << "Fetching reponse convolution tool." << endl;
  AdcChannelTool* pconvo = ptm->getShared<AdcChannelTool>("convo");
  AdcChannelTool* pconvo1 = ptm->getShared<AdcChannelTool>("convo1");
  AdcChannelTool* pdecon = nullptr;
  cout << myname << "Fetching deconvolution tool." << endl;
  string deconName;
  if ( iproc > 1 ) {
    deconName = iproc==5 ? "cmdecon" : iproc==4 ? "fmdecon" : iproc==3 ? "dmdecon" : "decon";
    pdecon = ptm->getShared<AdcChannelTool>(deconName);
    if ( pdecon == nullptr ) {
      cout << myname << "ERROR: Unable to find tool " << deconName << endl;
      return 2;
    }
  }
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
  // Charge is deposited in time-wire pixels with
  // We integrate across each 0.5 us time bin starting at the left edge of bin 0.
  string sttlPrefix;
  map<int, vector<float>> pixsamples;
  map<int, vector<float>> celsamples;
  int ncel = int(nchDep);
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
    // Loop over wire pixels.
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
    }  // End loop over wire pixels.
#else
  //   w is wire coordinate in wire spacing units
  //   t is drift coordinate in ticks (not time or distance)
  //   wpWid = wire pixel width [wire spacings]
  //   tpWid = time pixel width [ticks]
  //   Wire pixel 0 is at the left edge of the central cell.
  //   Time pixel 0 is tick 0
  //   (wtoff, ttoff) is a point on the track
  //   dtdw derived from dxdz is the track angle
    int iwir1 = -ncel*nbinPerWire;
    int iwir2 = (ncel + 1)*nbinPerWire;
    for ( int iwir=iwir1; iwir<iwir2; ++iwir ) {
      pixsamples[iwir].resize(nsam, 0.0);
    }
    sttlPrefix = "New ";
    float wpWid = 1.0/nbinPerWire;
    float tpWid = 1.0;
    float wpOff = 0.0;
    float dtdw = dxdz*dzcell/dxdt;    // dimensionless
    float dzds = cos(thtxz);
    float dxds = sin(thtxz);
    float dqdw = dqds*dzcell/dzds;
    float dqdt = dqds*dxdt/dxds;
    float wtOff = 0.0;
    float ttOff = isig;
    // Loop over wire pixels.
    for ( auto& kwir : pixsamples ) {
      int iwir = kwir.first;
      cout << myname << "Wire pixel " << iwir << endl;
      vector<float>& wsamples = kwir.second;
      // Handle track parallel to APA
      if ( dxdz == 0.0 ) {
        float qdep = dqdw*wpWid;
        wsamples[isig] = qdep;
        cout << myname << "  Q[" << iwir << "][" << isig << "] = " << qdep << endl;
        continue;
      }
      float xwir = wpOff + wpWid*iwir;
      float tmin = ttOff + dtdw*(iwir)*wpWid;
      float tmax = ttOff + dtdw*(iwir+1)*wpWid;
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
    }  // End loop over wire pixels.
    // Temporary: fill cells from pixels.
    int ipix = -ncel*nbinPerWire;
    if ( pixsamples.size() ) {
      for ( int icel=-ncel; icel<=ncel; ++icel ) {
        vector<float>& dest = celsamples[icel];
        for ( int ibin=0; ibin<nbinPerWire; ++ibin, ++ipix ) {
          std::transform(dest.begin(), dest.end(), pixsamples[ipix].begin(), dest.begin(), std::plus<float>());
        }
      }
    }
#endif
  }
  cout << myname << "Pixel sample counts: " << endl;
  for ( auto& kwir : pixsamples ) {
    int iwir = kwir.first;
    cout << myname << setw(6) << iwir << ": " << pixsamples[iwir].size() << endl;
  }
  cout << myname << "Cell sample counts: " << endl;
  for ( auto& kwir : celsamples ) {
    int iwir = kwir.first;
    cout << myname << setw(6) << iwir << ": " << celsamples[iwir].size() << endl;
  }

  // Assign chanel numbers, colors, styles, etc. for the three views.
  // Note channel number will be offset by the track angle.
  LineColors lc;
  vector<int> chans = {2100, 100, 1100};
  vector<int> mycols = {lc.blue(), lc.red(), lc.green()};
  vector<int> mywids = {2, 2, 4};
  vector<int> mystys = {1, 2, 3};
  vector<string> mylabs = {"X", "U", "V"};

  unsigned int ncon = 3;
  if ( iproc == 0 ) ncon = 1;
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
      // Fill samples from cells.
      for ( Index isam=0; isam<celsamples[0].size(); ++isam ) {
        if ( isam >= nsam ) break;
        acd.samples[isam] = celsamples[0][isam];
      }
      if ( pconvg != nullptr ) {
        cout << myname << "Adding longitudinal diffusion " << icha << endl;
        pconvg->update(acd).print();
      }
      if ( iproc > 0 ) {
        cout << myname << "Convoluting with detector response " << icha << endl;
        pconvo->update(acd).print();
        if ( apaConCells.size() ) {
          cout << myname << "Adding APA contributions from neighbors" << endl;
          for ( int icel : apaConCells ) {
            cout << myname << "...cell " << icel << endl;
            AdcChannelTool* pncon = abs(icel) == 1 ? pconvo1 : nullptr;
            if ( pncon == nullptr ) {
              cout << myname << "WARNING: Convolution tool not found for cell " << icel << endl;
              continue;
            }
            AdcChannelData acdtmp;
            acdtmp.channel = icha + icel;
            acdtmp.samples = celsamples[icel];
            pncon->update(acdtmp);
            for ( Index isam=0; isam<acdtmp.samples.size(); ++isam ) acd.samples[isam] += acdtmp.samples[isam];
          }
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
        if ( iproc > 1 ) {
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
      string sttl = iproc == 5 ? "Chi-square matrix deconvoluted reponse" :
                    iproc == 4 ? "Filtered matrix deconvoluted reponse" :
                    iproc == 3 ? "Direct matrix deconvoluted reponse" :
                    iproc == 2 ? "DFT deconvoluted reponse" :
                    iproc == 1 ? "Response to" : "Deposit from";
      sttl += " to 1 MIP for #theta_{xz}=" + stht + "#circ with " + slabsmr;
      if ( nchDep ) sttl += " (N_{neighbor} =" + to_string(nchDep) + ")";
      if ( sttlPrefix.size() ) sttl = sttlPrefix + sttl;
      man.setTitle(sttl);
      man.addHorizontalLine(0.0);
      man.setLabel("Run " + to_string(irun) + " event " + to_string(ievt));
      string fnam = "wf";
      fnam += (iproc>1 ? deconName : iproc ? "conv" : "noconv");
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
