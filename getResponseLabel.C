// getReponseLabelPeak.C
//
// Script that builds a label that includes the signed signal area and asymmetry.
// Signal is anythin above threhold.
//
void getResponseLabel(TH1* phk, string& slab, double ythr =0.0, int prec =3) {
  const string myname = "getResponseLabel: ";
  if ( ythr < 0.0 ) ythr = 0.0;
  int nbin = phk->GetNbinsX();
  ostringstream sslab;
  double aplu = 0.0;
  double amin = 0.0;
  for ( int ibin=1; ibin<=nbin; ++ibin ) {
    double val = phk->GetBinContent(ibin);
    if ( val > ythr ) aplu += val;
    if ( val < ythr ) amin += val;
  }
  //cout << myname << "ineg: " << ineg << endl;
  double area = aplu + amin;
  double asum = aplu - amin;
  double asym = asum > 0 ? (aplu + amin)/asum : 0.0;
  sslab << "A = " << std::fixed << std::setprecision(prec) << area;
  sslab << ", #Sigma|A| = " << std::fixed << std::setprecision(prec) << asum;
  sslab << ", #Phi = " << std::fixed << std::setprecision(3) << asym;
  slab += sslab.str();
}

