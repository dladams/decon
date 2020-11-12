int noise(float sig, unsigned int nval =1000) {
  vector<float> vals(nval,0);
  for ( float& val : vals ) val = sig*gRandom->Gaus();
  float sum0 = 0.0;
  float sum1 = 0.0;
  float sum2 = 0.0;
  for ( float val : vals ) {
    ++sum0;
    sum1 += val;
    sum2 += val*val;
  }
  float mean = sum1/sum0;
  float rms = sqrt(sum2/sum0 - mean*mean);
  cout << " input RMS: " << rms << endl;
  TH1* ph = new TH1F("h", "h", nval, 0, nval);
  int ibin = 0;
  for ( float val : vals ) {
    ph->SetBinContent(++ibin, val);
  }
  cesmearTH1(ph);
  ibin = 0;
  sum0 = 0;
  sum1 = 0;
  sum2 = 0;
  while ( ibin < nval ) {
    float val = ph->GetBinContent(++ibin);
    ++sum0;
    sum1 += val;
    sum2 += val*val;
  }
  mean = sum1/sum0;
  float oldrms = rms;
  rms = sqrt(sum2/sum0 - mean*mean);
  cout << "output RMS: " << rms << " (" << oldrms/rms << ")" << endl;
  return 0;
}
