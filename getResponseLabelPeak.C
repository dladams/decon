// getReponseLabelPeak.C
//
// Script that builds a label that includes the signed signal area and asymmetry.
// Signal lobes are found by looking for the first bin below the negative threshold.
//
void getResponseLabelPeak(TH1* phk, string& slab, double ythr =0.0, int prec =3) {
  const string myname = "getResponseLabel: ";
  int nbin = phk->GetNbinsX();
  int ineg=1;
  int imax = 0;
  double ymax = 0.0;
  while ( ineg<=nbin ) {
    double yval = phk->GetBinContent(ineg);
    if ( yval < -1.e-5 ) break;
    if ( yval > ymax ) {
      ymax = yval;
      imax = ineg;
    }
    ++ineg;
  }
  ostringstream sslab;
  int ithr1 = 1;
  int ithr2 = nbin;
  if ( ythr > 0.0 ) {
    ithr1 = ineg;
    for ( int ibin=ineg-1; ibin>=1; --ibin ) {
      if ( phk->GetBinContent(ibin) > ythr ) ithr1 = ibin;
    }
    if ( ineg <= nbin ) {
      ithr2 = ineg - 1;
      for ( int ibin=ineg; ibin<=nbin; ++ibin ) {
        if ( phk->GetBinContent(ibin) < -ythr ) ithr2 = ibin;
      }
    } else {
      ithr2 = ithr1;
      for ( int ibin=ithr1+1; ibin<=nbin; ++ibin ) {
        if ( phk->GetBinContent(ibin) > ythr ) ithr2 = ibin;
      }
    }
    sslab << "t=" << ythr << "(" << ithr1 << ", " << ineg << ", " << ithr2 << "): ";
  }
  //cout << myname << "ineg: " << ineg << endl;
  double aplu = phk->Integral(ithr1,ineg-1);
  double amin = phk->Integral(ineg, ithr2);
  double area = phk->Integral(ithr1, ithr2);
  double asum = aplu - amin;
  double asym = (aplu + amin)/asum;
  sslab << "A = " << std::fixed << std::setprecision(prec) << area;
  sslab << ", #Sigma|A| = " << std::fixed << std::setprecision(prec) << asum;
  sslab << ", #Phi = " << std::fixed << std::setprecision(3) << asym;
  slab += sslab.str();
}

