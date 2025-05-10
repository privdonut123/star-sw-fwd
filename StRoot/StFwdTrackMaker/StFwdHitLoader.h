
class StFwdHitLoader {
 public:
    StFwdHitLoader() {}
    ~StFwdHitLoader() {}
    void clear() {
        mFwdHitsFtt.clear();
        mFwdHitsFst.clear();
        mFwdHitsEpd.clear();
    }

    void loadFttHits( FwdDataSource::McTrackMap_t &mcTrackMap, FwdDataSource::HitMap_t &hitMap, int count ){
    void loadFttHitsFromStEvent( FwdDataSource::McTrackMap_t &mcTrackMap, FwdDataSource::HitMap_t &hitMap, int count ){
    void loadFttHitsFromGEANT( FwdDataSource::McTrackMap_t &mcTrackMap, FwdDataSource::HitMap_t &hitMap, int count ){


    int loadFstHits( FwdDataSource::McTrackMap_t &mcTrackMap, FwdDataSource::HitMap_t &hitMap );
    int loadFstHitsFromMuDst( FwdDataSource::McTrackMap_t &mcTrackMap, FwdDataSource::HitMap_t &hitMap );
    int loadFstHitsFromStEvent( FwdDataSource::McTrackMap_t &mcTrackMap, FwdDataSource::HitMap_t &hitMap );
    int loadFstHitsFromStRnDHits( FwdDataSource::McTrackMap_t &mcTrackMap, FwdDataSource::HitMap_t &hitMap );
    int loadFstHitsFromGEANT( FwdDataSource::McTrackMap_t &mcTrackMap, FwdDataSource::HitMap_t &hitMap );

    enum DataSource { GEANT, StEvent, MuDst };
    void setFttDataSource( DataSource ds ) { mFttDataSource = ds; }
    void setFstDataSource( DataSource ds ) { mFstDataSource = ds; }
    void setEpdDataSource( DataSource ds ) { mEpdDataSource = ds; }
    void setDataSource( DataSource ds ) { mFttDataSource = ds; mFstDataSource = ds; mEpdDataSource = ds; }
  protected:
    DataSource mFttDataSource;
    DataSource mFstDataSource;
    DataSource mEpdDataSource;

    // Pointers to these are used by StFwdTrackMaker, clear the vectors after each event
    vector<FwdHit> mFwdHitsFtt;
    vector<FwdHit> mFwdHitsFst;
    vector<FwdHit> mFwdHitsEpd;

};


