// Plot2dDft.h
//
// Class to plot power for a 2D DFT.
//

class Plot2dDft {

public:

  int makeHisto(const RealDftData& dft);

  TH2* getHisto() { return m_ph; }

private:

  TH2* m_ph = nullptr;

};

Plot2dDft::makeHisto(const RealDftData& dft) {
  Index n0 = dft.size(0);
  Index n1 = dft.size(1);
  m_ph = new TH2("hdft", "", n0, 0, n0, n1, 0, n1);
  float pwrsum = 0.0;
  typename C::IndexArray isams;
  Index& i0 = isams[0];
  Index& i1 = isams[1];
  for ( i1=0; i1<n1; ++i1 ) {
    cout << setw(4) << i1 << ":";
    for ( i0=0; i0<n0; ++i0 ) {
      float pwr = std::norm(dft.value(isams));
      pwrsum += pwr;
      cout << fixed << setprecision(5) << setw(10) << pwr;
      m_ph->SetBinContent(i0, i1, val);
    }
    cout << endl;
  }
}
