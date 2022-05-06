
#include "GenFit/FitStatus.h"
#include "GenFit/GFRaveVertexFactory.h"
#include <math.h>

typedef struct{
  float x;
  float y;
  float z;
}  Point;

float DEGS_TO_RAD = 3.14159f/180.0f;
float SCALE = 0.1;
// int numVertices = 0;    // Tallies the number of vertex points added.

//------------------------
//-- Prints a sphere as a "standard sphere" triangular mesh with the specified
//-- number of latitude (nLatitude) and longitude (nLongitude) lines and
//-- writes results to the specified output file (fout).





class ObjExporter {
public:

    void vert( ofstream &of, float x, float y, float z ){
        of << "v " << x << " " << y << " " << z << endl;
        numVertices++;
    }

    void sphere(TVector3 pt, float radius, int nLatitude, int nLongitude, ofstream &fout){
        int p, s, i, j;
        float x, y, z, out;
        int nPitch = nLongitude + 1;

        float pitchInc = (180. / (float)nPitch) * DEGS_TO_RAD;
        float rotInc   = (360. / (float)nLatitude) * DEGS_TO_RAD;

        //## PRINT VERTICES:

        vert(fout, pt.X(), pt.Y()+radius, pt.Z());    // Top vertex.
        vert(fout, pt.X(), pt.Y()-radius, pt.Z());    // Bottom vertex.
        // numVertices = numVertices+2;

        int fVert = numVertices;    // Record the first vertex index for intermediate vertices.
        for(p=1; p<nPitch; p++) {    // Generate all "intermediate vertices":
            out = radius * sin((float)p * pitchInc);
            if(out < 0) out = -out;    // abs() command won't work with all compilers
            y   = radius * cos(p * pitchInc);
            // printf("OUT = %g\n", out);    // bottom vertex
            // printf("nPitch = %d\n", nPitch);    // bottom vertex
            for(s=0; s<nLatitude; s++) {
                x = out * cos(s * rotInc);
                z = out * sin(s * rotInc);
                fout << TString::Format("v %g %g %g\n", x+pt.X(), y+pt.Y(), z+pt.Z()).Data();
                numVertices++;
            }
        }
      
      //## PRINT SQUARE FACES BETWEEN INTERMEDIATE POINTS:
      
      for(p=1; p<nPitch-1; p++) {
        for(s=0; s<nLatitude; s++) {
          i = p*nLatitude + s;
          j = (s==nLatitude-1) ? i-nLatitude : i;
          fout << TString::Format("f %d %d %d %d\n", 
                  (i+1-nLatitude)+fVert, (j+2-nLatitude)+fVert, (j+2)+fVert, (i+1)+fVert);
        }
      }
      
      //## PRINT TRIANGLE FACES CONNECTING TO TOP AND BOTTOM VERTEX:
      
      int offLastVerts  = fVert + (nLatitude * (nLongitude-1));
      for(s=0; s<nLatitude; s++)
      {
        j = (s==nLatitude-1) ? -1 : s;
        fout << TString::Format("f %d %d %d\n", fVert-1, (j+2)+fVert,        (s+1)+fVert       ).Data();
        fout << TString::Format("f %d %d %d\n", fVert,   (s+1)+offLastVerts, (j+2)+offLastVerts).Data();
      }
    }

    static TVector3 trackPosition( genfit::Track * t, float z, float * cov = 0 ){

        try {
            auto plane = genfit::SharedPlanePtr(
                // these normals make the planes face along z-axis
                new genfit::DetPlane(TVector3(0, 0, z), TVector3(1, 0, 0), TVector3(0, 1, 0) )
            );

            genfit::MeasuredStateOnPlane tst = t->getFittedState(1);
            auto TCM = t->getCardinalRep()->get6DCov(tst);
            //  can get the track length if needed
            // double len = t->getCardinalRep()->extrapolateToPlane(tst, detSi, false, true);
            double len = t->getCardinalRep()->extrapolateToPlane(tst, plane, false, true);

            TCM = t->getCardinalRep()->get6DCov(tst);

            // can get the projected positions if needed
            float x = tst.getPos().X();
            float y = tst.getPos().Y();
            float _z = tst.getPos().Z();

            if ( cov ){
                cov[0] = TCM(0,0); cov[1] = TCM(1,0); cov[2] = TCM(2,0);
                cov[3] = TCM(0,1); cov[4] = TCM(1,1); cov[5] = TCM(2,1);
                cov[6] = TCM(0,2); cov[7] = TCM(1,2); cov[8] = TCM(2,2);
            }


            return TVector3( x, y, _z );
        } catch ( genfit::Exception e ){
            LOG_INFO << "E: " << e.what() << endm;
            return TVector3( -990, -990, -990 );
        }

        return TVector3( -99, -99, -99 );
    }
    
    void output( std::string filename, 
            std::vector< Seed_t> seeds, 
            std::vector< genfit::Track *> tracks, 
            const std::vector< genfit::GFRaveVertex *> &vertices, 
            std::vector<TVector3> &fttHits,
            std::vector<TVector3> &fstHits,
            std::vector<TVector3> &fcsPreHits, // EPD = preshower
            std::vector<TVector3> &fcsClusters ){
		
		LOG_INFO << "Writing: " << filename << endm;
        numVertices = 0;
        // OPEN output
        ofstream ofile( (filename + ".obj" ).c_str() );

        ofile << "\nmtllib materials.mtl\n\n" << endl;


        TVector3 startPos;
        if ( vertices.size() > 0 ){
            ofile << "o FwdVertices" << endl;
            size_t ivert = 0;
            for ( auto v : vertices ) {
                // ofile << "v " << v->getPos().X() << " " << v->getPos().Y() << " " << v->getPos().Z() << endl;
                startPos.SetXYZ( v->getPos().X(), v->getPos().Y(), v->getPos().Z() );
                // ivert++;
                sphere( TVector3( v->getPos().X() * SCALE, v->getPos().Y() * SCALE, -v->getPos().Z() * SCALE ), 0.5, 10, 10, ofile );
            }
        }

        if ( fttHits.size() > 0 ){
            ofile << "\n" << endl;
            ofile << "o fttHits" << endl;
            ofile << "usemtl stgc_hits\n" << endl;
            for ( auto p : fttHits ){
                sphere( TVector3( p.X() * SCALE, p.Y() * SCALE, -p.Z() * SCALE ), 0.15, 12, 12, ofile );
            }
        }


        if ( fstHits.size() > 0 ) {
            ofile << "\n" << endl;
            ofile << "o fstHits" << endl;
            ofile << "usemtl fst_hits\n" << endl;
            for ( auto p : fstHits ){
             
                sphere( TVector3( p.X() * SCALE, p.Y() * SCALE, -p.Z() * SCALE ), 0.15, 10, 10, ofile );
            }
        }

        if ( fcsPreHits.size() > 0 || fcsClusters.size() > 0 ){ 
            ofile << "\n" << endl;
            ofile << "o fcs" << endl;
            ofile << "usemtl fcs_hits\n" << endl;
            for ( auto p : fcsPreHits ){
             
                sphere( TVector3( p.X() * SCALE, p.Y() * SCALE, -p.Z() * SCALE ), 0.25, 10, 10, ofile );
            }
            ofile << "\n\n";
            for ( auto p : fcsClusters ){
             
                sphere( TVector3( p.X() * SCALE, p.Y() * SCALE, -p.Z() * SCALE ), 0.75, 10, 10, ofile );
            }
        }
        

        if ( seeds.size() > 0 ){
            ofile << "\n\no FwdSeeds\n" << endl;
            ofile << "usemtl seeds\n" << endl;
            // numVertices = 0;
            for ( auto s : seeds ) {
                size_t vStart = numVertices;
                for ( auto h : s ){

                    vert( ofile, h->getX() * SCALE, h->getY() * SCALE, -h->getZ() * SCALE );

                }

                ofile << "l ";
                for ( size_t i = vStart; i < numVertices; i++){
                    ofile << i+1 << " "; 
                }
                ofile << endl;
            }
        }

        // ofile.open( (filename + "_tracks.obj" ).c_str()  );
        // numVertices = 0;
        // EXPORT TRACKS

        if ( tracks.size() > 0 ){
            ofile << "\n\no FwdTracks\n" << endl;
            ofile << "usemtl tracks\n" << endl;
            float zStep = 4.0; // cm
            for ( auto t : tracks ) {
                size_t vStart = numVertices;
    
                if ( false ) {
    
                    TVector3 point;
                    point = startPos;
                    vert( ofile, point.X() * SCALE, point.Y() * SCALE, -point.Z() * SCALE );
    
                    point = trackPosition( t, 281.0 );
                    if ( point.X() > -90 && point.Y() > -90 ) 
                        vert( ofile, point.X() * SCALE, point.Y() * SCALE, -point.Z() * SCALE );
    
                    point = trackPosition( t, 304.0 );
                    if ( point.X() > -90 && point.Y() > -90 )
                        vert( ofile, point.X() * SCALE, point.Y() * SCALE, -point.Z() * SCALE );
    
                    point = trackPosition( t, 327.0 );
                    if ( point.X() > -90 && point.Y() > -90 )
                        vert( ofile, point.X() * SCALE, point.Y() * SCALE, -point.Z() * SCALE );
    
                    point = trackPosition( t, 349.5 );
                    if ( point.X() > -90 && point.Y() > -90 )
                        vert( ofile, point.X() * SCALE, point.Y() * SCALE, -point.Z() * SCALE );
    
                    point = trackPosition( t, 800.0 );
                    if ( point.X() > -90 && point.Y() > -90 )
                        vert( ofile, point.X() * SCALE, point.Y() * SCALE, -point.Z() * SCALE );
    
                } else {
                    TVector3 lpoint;
                    for ( float z = startPos.Z(); z < 800; z += zStep ){
                        TVector3 point = trackPosition( t, z );
                        // if ( point.X() == -99 && point.Y() == -99 && point.Z() == -99 ){
                        //     point = lpoint;
                        // }
    
                        if ( point.X() < -900 && point.Y() < -900 ) break;
                        if ( point.X() < -90 && point.Y() < -90 ) { z+= 50; continue;}
    
                        vert( ofile, point.X() * SCALE, point.Y() * SCALE, -point.Z() * SCALE );
                        // ofile << "v " << point.X() << " " << point.Y() << " " << point.Z() << endl;
                        // ivert ++;
                        lpoint = point;
                    }
                }
                
                ofile << "l ";
                for ( size_t i = vStart; i < numVertices; i++){
                    ofile << i+1 << " "; 
                }
                ofile << endl;
            }
        }


        ofile.close();	
    }

    size_t numVertices;

};
