// Display gaus values.
void dogaus(double sig) {
  double xmax = int(4*sig);
  for ( double x=-xmax; x<=xmax; x+=1.0 ) cout << setw(4) << TMath::Gaus(x, 0, sig, true) << ",";
  cout << endl;
}
