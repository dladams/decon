// Signal extends over ntick tiks smeared with sconvg.
// act = 0: charge deposit (no
int run(float thtxzdeg, string ssig, unsigned int iconv =1, float noise = 0.0, unsigned int seed =0) {
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
  if ( ssig.size() && ssig != "0" ) {
    string sconvg = "convg" + ssig;
    cout << myname << "Fetching filter smearing tool " << sconvg << endl;
    pconvg = ptm->getShared<AdcChannelTool>(sconvg);
    if ( pconvg == nullptr ) {
      cout << myname << "ERROR: Unable to find tool " << sconvg << endl;
      return 2;
    }
    slabsmr = "#sigma_{smear} = " + ssig + " ticks";
  }
  vector<float> samples(100, 0.0);
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
  cout << myname << "Sample count: " << samples.size() << endl;
  TPadManipulator man;
  string hopt = "h";
  LineColors lc;
  vector<int> chans = {1600, 0, 800};
  vector<int> mycols = {lc.blue(), lc.red(), lc.green()};
  vector<int> mywids = {2, 2, 4};
  vector<int> mystys = {1, 2, 3};
  unsigned int ncon = 3;
  if ( iconv == 0 ) ncon = 1;
  vector<string> slabs;
  Index nsam = 1000;
  Index isig = 500;
  for ( unsigned int icon=0; icon<ncon; ++icon ) {
    Index icha = chans[icon];
    AdcChannelData acd;
    acd.run = 100;
    acd.event = 1;
    acd.sampleUnit = "ke/tick";
    acd.channel = icha;
    acd.samples.resize(nsam);
    Index isam = isig;
    for ( Index jsam=0; jsam<samples.size(); ++jsam ) {
      if ( isam >= nsam ) break;
      acd.samples[isam] = samples[jsam];
      ++isam;
cout << "CCCCCCCCCC " << isam << ": " << acd.samples[isam] << endl;
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
      }
      if ( iconv > 1 ) {
        cout << myname << "Deconvoluting channel " << icha << endl;
        pdecon->update(acd).print();
      }
    }
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
    getResponseLabel(ph, slab, 0.0, 1);
    slabs.push_back(slab);
    cout << myname << slab << endl;
    man.add(ph, hopt);
    hopt = "h same";
  }
  man.addAxis();
  if ( iconv != 2 ) man.setRangeX(isig-210, isig+250);
  else man.setRangeX(isig-300, isig+300);
  float ymin = -5;
  float ymax = 35.0;
  if ( iconv > 0 ) {
    ymin = -1.5;
    ymax = 6.0;
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
  fnam += ".png";
  cout << myname << "Printing " << fnam << endl;
  man.print(fnam);
  return 0;
}
