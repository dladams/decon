class CeNoise {
  using Index = unsigned int;
  TRandom ran;
  Index nran = 10000;
  Index iran = nran;
  float scale = 0.0;
  vector<double> cache;
public:
  CeNoise(Index seed) : ran(seed), cache(nran) { }
  double get() {
    if ( iran >= nran ) {
      for ( double& val : cache ) val = ran.Gaus();
      cesmearVector(cache);
      iran = 0;
      if ( scale == 0.0 ) {
        float sum0 = 0.0;
        float sum1 = 0.0;
        float sum2 = 0.0;
        for ( float val : cache ) {
          ++sum0;
          sum1 += val;
          sum2 += val*val;
        }
        float mean = sum1/sum0;
        float rms = sqrt(sum2/sum0 - mean*mean);
        scale = 1/rms;
      }
    }
    return scale*cache[iran++];
  }
};

