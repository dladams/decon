// Signal extends over ntick tiks smeared with sconvg.
// act = 0: charge deposit (no
// ssig:
//      "" - no smearing
//      "1" - smear RMS = 1 ==> 1 m
//      "2" - smear RMS = 2 ==> 3.6 m
//  "nosig" - No signal (noise only)
// iconv:
//   0 - energy deposit
//   1 - convoluted with response
//   2 - deconvoluted
int run(float thtxzdeg, string ssig, int iconv =1, float noise = 0.0, unsigned int seed =0, unsigned int nevt =1) {
  using Index = unsigned int;
  // Scale noise for amplifier suppression.
  CeNoise cen(seed);
  string myname = "run: ";
  if ( iconv > 2 ) {
    cout << myname << "Invalid conv option: " << iconv << endl;
    return 1;
  }
  cout << myname << "Fetching waveform plotter." << endl;
  AdcChannelTool* pwf = ptm->getShared<AdcChannelTool>("wf");
  cout << myname << "Fetching reponse convolution tool." << endl;
  AdcChannelTool* pconvo = ptm->getShared<AdcChannelTool>("convo");
  cout << myname << "Fetching deconvolution tool." << endl;
  AdcChannelTool* pdecon = ptm->getShared<AdcChannelTool>("decon");
  AdcChannelTool* pconvg = nullptr;
  string slabsmr = "no smearing";
  bool doSig = ssig != "nosig";
  // Fetch tools.
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
  AdcChannelTool* psigfind = ptm->getShared<AdcChannelTool>("adcThresholdSignalFinder");
  if ( psigfind == nullptr ) return 3;
  cout << myname << "Fetching sample RMS calculator." << endl;
  AdcChannelTool* psamrms = ptm->getShared<AdcChannelTool>("adcChannelSamplelRmsPlotter");
  if ( psamrms == nullptr ) return 3;
  cout << myname << "Fetching signal RMS calculator." << endl;
  AdcChannelTool* psigrms = ptm->getShared<AdcChannelTool>("adcChannelSignalRmsPlotter");
  if ( psigrms == nullptr ) return 3;
  vector<float> samples(100, 0.0);
  if ( doSig ) {
    float dz = 4.79;    // mm
    float dxdt = 0.8;  // mm/tick
    float dqds = 6.0;  // ke/mm
    float thtxz = thtxzdeg*acos(-1.0)/180.0;
    float dxdz = tan(thtxz);
    float dx = dxdz*dz;
    float dt = dx/dxdt;
    float ds = sqrt(dz*dz + dx*dx);
    float dq = dqds*ds;
    cout << myname << "  ds = " << ds << " mm" << endl;
    cout << myname << "  dt = " << dt << " tick" << endl;
    cout << myname << "  dq = " << dq << " ke" << endl;
    unsigned int itick = 0;
    float dtrem = dt;
    float dqrem = dq;
    if ( dt > 1.0 ) {
      float dqdt = dq/dt;
      while ( dtrem >= 1.0 ) {
        samples[itick] = dqdt;
        cout << myname << "  Q[" << itick << "] = " << dqdt << endl;
        ++itick;
        dtrem -= 1.0;
        dqrem -= dqdt;
      }
    }
    if ( dqrem > 0.0 ) {
      samples[itick] = dqrem;
      cout << myname << "  Q[" << itick << "] = " << dqrem << endl;
    }
  }
  cout << myname << "Sample count: " << samples.size() << endl;
  TPadManipulator man;
  string hopt = "h";
  LineColors lc;
  vector<int> chans = {2000, 0, 1000};
  vector<int> mycols = {lc.blue(), lc.red(), lc.green()};
  vector<int> mywids = {2, 2, 4};
  vector<int> mystys = {1, 2, 3};
  vector<string> mylabs = {"X", "U", "V"};
  unsigned int ncon = 3;
  vector<TH1*> myhsts(3, nullptr);
  if ( iconv == 0 ) ncon = 1;
  vector<string> slabs;
  Index nsam = 1000;
  Index isig = 500;
  for ( Index ievt=1; ievt<nevt+1; ++ievt ) {
    for ( unsigned int icon=0; icon<ncon; ++icon ) {
      Index icha = chans[icon] + thtxzdeg;
      Index irun = 100;
      DuneEventInfo evt;
      evt.run = irun;
      evt.event = ievt;
      AdcChannelDataMap acm;
      AdcChannelData& acd = acm[icha];
      acd.run = irun;
      acd.event = ievt;
      acd.sampleUnit = "ke/tick";
      acd.channel = icha;
      acd.samples.resize(nsam);
      Index isam = isig;
      for ( Index jsam=0; jsam<samples.size(); ++jsam ) {
        if ( isam >= nsam ) break;
        acd.samples[isam] = samples[jsam];
        ++isam;
      }
      if ( iconv > 0 ) {
        if ( pconvg != nullptr ) {
          cout << myname << "Smearing channel " << icha << endl;
          pconvg->update(acd).print();
        }
        cout << myname << "Convoluting channel " << icha << endl;
        pconvo->update(acd).print();
        if ( noise > 0.0 ) {
          for ( float& sam : acd.samples ) sam += noise*cen.get();
          acd.dftmags.clear();
          acd.dftphases.clear();
        }
        pfft->update(acd);
        pplotDftBefore->viewMap(acm);
        pplotDftBefore->endEvent(evt);
        if ( iconv > 1 ) {
          cout << myname << "Deconvoluting channel " << icha << endl;
          pdecon->update(acd).print();
          pplotDftAfter->viewMap(acm);
          pplotDftAfter->endEvent(evt);
        }
      }
      psigfind->update(acd);
      for ( Index isam=0; isam<acd.signal.size(); ++isam ) acd.signal[isam] = acd.signal[isam] ? false : true;
      DataMap dsamrms = psamrms->view(acd);
      cout << myname << "Sample RMS: " << dsamrms.getFloat("metricValue") << endl;
      DataMap dsigrms = psigrms->view(acd);
      cout << myname << "Signal RMS: " << dsigrms.getFloat("metricValue") << endl;
      //dsamrms.print();
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
      if ( nevt < 2 ) {
        getResponseLabel(ph, slab, thr, 1);
        slab = mylabs[icon] + ": " + slab;
        slabs.push_back(slab);
        cout << myname << slab << endl;
      }
      man.add(ph, hopt);
      hopt = "h same";
    }
  }
  man.addAxis();
  man.setRangeX(isig- 50, isig+150);
  float ymin = -5;
  float ymax = 35.0;
  if ( iconv > 0 ) {
    ymin = -1.5;
    ymax = 6.0;
  }
  if ( nevt > 1 ) {
    ymin = -5;
    ymax = 10;
  }
  man.setRangeY(ymin, ymax);
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
  string sttl = iconv == 2 ? "Deconvoluted reponse" :
                iconv == 1 ? "Response to" : "Deposit from";
  sttl += " to 1 MIP for #theta_{xz}=" + stht + "#circ with " + slabsmr;
  man.setTitle(sttl);
  man.addHorizontalLine(0.0);
  string fnam = "wf";
  fnam += (iconv==2 ? "deconv" : iconv ? "conv" : "noconv");
  fnam += "-sigma" + ssig;
  fnam += "-thtxz" + stht;
  string snsam = to_string(nsam);
  while ( snsam.size() < 5 ) snsam = "0" + snsam;
  fnam += "nsam" + snsam;
  fnam += ".{png,tpad}";
  cout << myname << "Printing " << fnam << endl;
  man.print(fnam);
  cout << myname << "Closing tools." << endl;
  pplotDftBefore.reset();
  pplotDftAfter.reset();
  cout << myname << "Exiting." << endl;
  return 0;
}
