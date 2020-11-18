class CeNoise {
public:
  using Index = unsigned int;
  using Float = double;
private:
  TRandom ran;
  Index nran;
  Index iran;
  Index icnt;
  Float m_scale;
  vector<Float> cache;
  Float m_value;
  Float* address;
public:
  CeNoise(Index seed) :
  nran(100000), iran(0), icnt(0), m_scale(0.0),
  ran(seed), cache(nran), m_value(0.0), address(&cache.front()) { }
  Float get() {
    if ( address != &cache.front() ) cout << "ERROR: Address changed." << endl;
    if ( iran%nran == 0 ) {
      for ( Index ival=0; ival<cache.size(); ++ival ) {
        Float& val = cache[ival];
        val = ran.Gaus();
        if ( fabs(val) > 11 ) cout << "ERROR: Invalid raw entry cache[ " << ival << "] = " << cache[ival] << endl;
      }
      cesmearVector(cache);
      for ( Index ival=0; ival<cache.size(); ++ival ) {
        Float val = cache[ival];
        if ( fabs(val) > 11 ) {
          cout << "ERROR: Invalid smeared entry cache[ " << ival << "] = " << cache[ival] << endl;
          return 0.0;
        }
      }
      if ( address != &cache.front() ) cout << "ERROR: Address changed after smear." << endl;
      iran = 0;
      if ( scale() == 0.0 ) {
        Float sum0 = 0.0;
        Float sum1 = 0.0;
        Float sum2 = 0.0;
        for ( Float val : cache ) {
          ++sum0;
          sum1 += val;
          sum2 += val*val;
        }
        Float mean = sum1/sum0;
        Float rms = sqrt(sum2/sum0 - mean*mean);
        m_scale = 1/rms;
      }
    }
    ++icnt;
    m_value = scale()*cache[iran++];
    return m_value;
  }
  Index cacheSize() const { return nran; }
  double scale() const { return m_scale; }
  Index next() const { return iran; }
  Index count() const { return icnt; }
  double last() const { return m_value; }
};

int test_CeNoise(unsigned int seed =1) {
  const string myname = "test_CeNoise: ";
  CeNoise cen(seed);
  using Index = CeNoise::Index;
  for ( Index i=0; i<5*cen.cacheSize(); ++i ) {
    if ( i%100 == 0 ) cout << myname << setw(5) << i << ": " << cen.get() << endl;
  }
  cout << myname << "Scale: " << cen.scale() << endl;
  return 0;
}
