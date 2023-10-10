#include <string>
#include <iostream>

#include "TString.h"
#include <TApplication.h>
#include "TSystem.h"

// Objects
#include "CommandLineParser.h"
#include "SEventHandle.h"
#include "SParserConfig.h"
#include "STPCAnalyzerTreeReader.h"
#include "TPCSimpleEvents.h"
#include "TPCSimpleHits.h"
#include "TPCLinesParameters.h"
#include "TPCLinesAlgo.h"
#include "ChargeDensity.h"
#include "ChargeDensityPset.h"
#include "FRANSTTreeHandle.h"


std::vector<SHit> GetFRANSHitsView(
    std::string view,
    std::vector<int> *_X,
    std::vector<double> *_Y,
    std::vector<double> *_Int,
    std::vector<double> *_Wi,
    std::vector<double> *_ST,
    std::vector<double> *_ET)
{

    // Channel boundaries
    std::map<std::string, std::vector<int>> fChB =  {
        {"U0", {0, 1983}},
        {"V0", {1984, 3967}},
        {"C0", {3968, 5631}},
        {"U1", {5632, 7631}},
        {"V1", {7632, 9599}},
        {"C1", {9600, 11263}}
    };

    // reset variables
    std::vector<SHit> hitList;

    int nTotalHits = _X->size();

    for (int i = 0; i < nTotalHits; i++) {

        // filter channels for the view
        int x = _X->at(i);        
        
        if ( x > fChB[view][0] && x <= fChB[view][1]) {
            SHit hit(-1, _X->at(i), _Y->at(i), _Wi->at(i), _Int->at(i), _ST->at(i), _ET->at(i));
            hitList.push_back(hit);
        }
    }

    return hitList;

}


void RunAlgoFRANS(const CommandLineParser& parser)
{

    int Debug = parser.getDebug();
    std::string ConfPsetPath = parser.getPsetPath();
    if(ConfPsetPath==""){
        ConfPsetPath="config.txt";
    }
    int DebugMode = parser.getDebugMode();
    int n = parser.getN();
    int nskip = parser.getNskip();
    int event = parser.getEvent();
    int sr = parser.getSr();
    std::string file_name = parser.getFileName();
    std::string directory_path = parser.getDirectoryPath();
    std::string ext = parser.getExtension();
    int vertexOption = parser.getVertexOption();

    double fFRANSScoreCut = 0.15;

    // Input tree
    std::string fTreeName = "ana/AnaTPCTree"; 

    // Set batch mode
    if(Debug==0) gROOT->SetBatch(true);
    
    // Get input files
    std::vector<TString> fFilePaths = GetInputFileList(file_name, ext, directory_path);

    // ----------- ALGORITHM PARAMETERS --------------------------------- 

    // ---- FRAMS parameters ----------------------------------------
    FRAMSPsetType fPsetFRANS = ReadFRANSPset(ConfPsetPath);
    fPsetFRANS.Verbose = Debug;
    // ---- TPCLines parameter ----------------------------------------
    // Parameter sets
    TrackFinderAlgorithmPsetType fPsetTrackFinder = ReadTrackFinderAlgorithmPset(ConfPsetPath);
    fPsetTrackFinder.Verbose = Debug;
    HoughAlgorithmPsetType fPsetHough = ReadHoughAlgorithmPset(ConfPsetPath);
    fPsetHough.Verbose = Debug;
    VertexFinderAlgorithmPsetType fPsetVertexFinder = ReadVertexFinderAlgorithmPset(ConfPsetPath);
    fPsetVertexFinder.Verbose = Debug;
    TPCLinesAlgoPsetType fPsetAnaView = ReadTPCLinesAlgoPset(ConfPsetPath);
    fPsetAnaView.Verbose = Debug;
    fPsetAnaView.DebugMode = DebugMode;
    fPsetAnaView.HoughAlgorithmPset = fPsetHough;
    fPsetAnaView.TrackFinderAlgorithmPset = fPsetTrackFinder;
    fPsetAnaView.VertexFinderAlgorithmPset = fPsetVertexFinder;

    // View to use
    std::string fView = fPsetAnaView.View;

    
    // Define the program control variables
    int fNEv = n;
    int fEv = event;
    int fSubRun = sr;
    int fNEvSkip = nskip;
    int nEvents=0;


    // Output directory for display
    gSystem->Exec(("rm -rf "+fPsetFRANS.OutputPath).c_str());
    gSystem->Exec( ("mkdir "+fPsetFRANS.OutputPath).c_str());
    gSystem->Exec( ("mkdir "+fPsetFRANS.OutputPath+"/rootfiles").c_str());
    

    // Output directory for analysis results
    std::string fAnaReultsPath = "anaResults";
    gSystem->Exec(("rm -rf "+fAnaReultsPath).c_str());
    gSystem->Exec(("mkdir "+fAnaReultsPath).c_str());
    std::string tree_dirname = "framsReco/";
    if(vertexOption==1) tree_dirname = "framsTrue/";
    if(vertexOption==2) tree_dirname = "framsMine/";
    std::string tree_name = "FRAMSTree";
    TTree* fTree = new TTree((tree_name).c_str(), "FRAMS Output Tree");
    FRANSTTree myTree(fTree);


    // Define FRANS LINES ALGORITHM
    ChargeDensity _FRAMSAlgo(fPsetFRANS);
    ChargeDensity _FRAMSAlgoPANDORA(fPsetFRANS);
    // Define TPC LINES ALGORITHM
    TPCLinesAlgo _TPCLinesAlgo(fPsetAnaView);
    // Effiency status
    EfficiencyCalculator _EfficiencyCalculator(fAnaReultsPath);

    // TTree loop
    int nEntries = 0;
    for (const auto& filepath : fFilePaths) {

        std::cout<<" ANALYZING THE FILE: "<<filepath<<std::endl;
        MyTPCTreeReader treeReader(filepath, fTreeName);
        

        for (int entry = 0; entry < treeReader.NEntries(); entry++) {
            std::cout<<"Entry "<<entry<<std::endl;
            // get the entry
            treeReader.GetEntry(entry);

            // check program control variables
            SEventId ev(treeReader.runID, treeReader.subrunID, treeReader.eventID);
            if (fNEv > 0 && nEvents >= fNEv) continue;
            if (treeReader.eventID != fEv && fEv != -1) continue;
            if (treeReader.subrunID != fSubRun && fSubRun != -1) continue;
            nEntries++;
            if (nEntries <= fNEvSkip) continue;
            
            // Check if it's signal
            bool isSignal=false;
            std::cout<<"Perl\n";
            for(size_t k=0; k<treeReader.truePrimeriesPDG->size();k++){
                if(treeReader.truePrimeriesPDG->at(k)==3122) isSignal=true;
            }

            nEvents++;
            std::cout << "\n\n ************** Analyzing: " << ev;
            
            // True vertex
            std::vector<double> VertexXYZ = {treeReader.nuvX, treeReader.nuvY, treeReader.nuvZ};
            std::vector<int> VertexUVYT = {treeReader.nuvU, treeReader.nuvV, treeReader.nuvC, treeReader.nuvTimeTick};
            std::cout << "   - True vertex (X, Y, Z, T) " << VertexXYZ[0] << " " << VertexXYZ[1] << " " << VertexXYZ[2] << " " << treeReader.nuvT << " (U, V, C, TT): " << VertexUVYT[0] << " " << VertexUVYT[1] << " " << VertexUVYT[2] << " " << VertexUVYT[3] << std::endl;

            // Reco vertex
            std::vector<double> RecoVertexXYZ = {treeReader.recnuvX, treeReader.recnuvY, treeReader.recnuvZ};
            std::vector<int> RecoVertexUVYT = {treeReader.recnuvU, treeReader.recnuvV, treeReader.recnuvC, treeReader.recnuvTimeTick};
            std::cout << "  - Reco vertex (X, Y, Z) " << RecoVertexXYZ[0] << " " << RecoVertexXYZ[1] << " " << RecoVertexXYZ[2] << " (U, V, C, TT): " << RecoVertexUVYT[0] << " " << RecoVertexUVYT[1] << " " << RecoVertexUVYT[2] << " "<< RecoVertexUVYT[3]  << std::endl;

            // Set TPC using the drift coordinate of the reco vertex
            //TPC = (nuvX > 0) ? 1 : 0;
            int TPC = (RecoVertexXYZ[0] > 0) ? 1 : 0;
            
            // We need a minimum number of hits to run the track finder
            size_t nhits = treeReader.hitsChannel->size();
            std::cout << "  - NHits: " << nhits << std::endl;
            if(nhits<=3){
                std::cout<<"   SKIPPED NHits \n";
                _EfficiencyCalculator.UpdateSkipped(ev);
                continue;
            }
            if( (vertexOption==0 || vertexOption==2) and treeReader.recnuvC==-1){
                std::cout<<"   SKIPPED RecoVertex \n";
                _EfficiencyCalculator.UpdateSkipped(ev);
                continue;
            }
            if(vertexOption==1 and treeReader.nuvC==-1){
                std::cout<<"   SKIPPED TrueVertex \n";
                _EfficiencyCalculator.UpdateSkipped(ev);
                continue;
            }

            // Assing the view and TPC
            std::string view;
            if(fView=="Best"){
                std::string bestView=_TPCLinesAlgo.GetBestView(treeReader.hitsView, treeReader.hitsChi2);
                std::cout<<"  Using best view: "<<bestView+std::to_string(TPC)<<std::endl;
                view = bestView+std::to_string(TPC);
            }    
            else{
                view = fView+std::to_string(TPC);
            }
            
            // Set the hits
            std::vector<SHit> hitList = GetFRANSHitsView(view,
                                                        treeReader.hitsChannel,
                                                        treeReader.hitsPeakTime,
                                                        treeReader.hitsIntegral, 
                                                        treeReader.hitsRMS,
                                                        treeReader.hitsStartT, 
                                                        treeReader.hitsEndT);
                                        
            // Set the vertex
            // true
            double vertexXTrue = VertexUVYT[2];
            if (view == "U0" || view == "U1") vertexXTrue = VertexUVYT[0];
            if (view == "V0" || view == "V1") vertexXTrue = VertexUVYT[1];
            double vertexYTrue = VertexUVYT[3];
            SVertex fVertexTrue(SPoint(vertexXTrue, vertexYTrue), view);
            
            // reco
            double vertexXReco = RecoVertexUVYT[2];
            if (view == "U0" || view == "U1") vertexXReco = RecoVertexUVYT[0];
            if (view == "V0" || view == "V1") vertexXReco = RecoVertexUVYT[1];
            double vertexYReco = RecoVertexUVYT[3];
            SVertex fVertexReco(SPoint(vertexXReco, vertexYReco), view);

            SVertex fVertexMine;
            
            SEvent recoEvent;
            if(vertexOption==2){
                // Set the hits
                bool filled = _TPCLinesAlgo.SetHitList(view, RecoVertexUVYT, VertexUVYT, 
                                        treeReader.hitsChannel,
                                        treeReader.hitsPeakTime,
                                        treeReader.hitsIntegral, 
                                        treeReader.hitsRMS,
                                        treeReader.hitsStartT, 
                                        treeReader.hitsEndT,
                                        treeReader.hitsChi2,
                                        "");

                if(filled){
                    std::cout<<" Analyzing TPCLines\n";
                    // Analyze
                    _TPCLinesAlgo.AnaView(ev.Label());
                    recoEvent = _TPCLinesAlgo.GetRecoEvent();
                    fVertexMine = SVertex( SPoint((double)_TPCLinesAlgo.GetMainVertex().X()+_TPCLinesAlgo.ShiftX(),
                                                (double) _TPCLinesAlgo.GetMainVertex().Y()+_TPCLinesAlgo.ShiftY())
                                            , "");
                }
                else{
                    fVertexMine = fVertexReco;
                }
            }

            // set it
            SVertex fVertex = fVertexReco;
            if(vertexOption==1) fVertex = fVertexTrue;
            else if(vertexOption==2) fVertex = fVertexMine;

            
            _FRAMSAlgoPANDORA.Fill(hitList, fVertexReco);
            _FRAMSAlgo.Fill(hitList, fVertex);

            myTree.FillData(2, treeReader.runID, treeReader.subrunID, treeReader.eventID, isSignal,
                            _FRAMSAlgo.Delta(), _FRAMSAlgo.Eta(), _FRAMSAlgo.FitScore(), _FRAMSAlgo.Alpha(),
                            _FRAMSAlgo.Omega(), _FRAMSAlgo.Tau(), _FRAMSAlgo.Iota(),
                            recoEvent.GetNOrigins(), recoEvent.GetNOriginsMult(1), recoEvent.GetNOriginsMult(2), recoEvent.GetNOriginsMultGt(3), recoEvent.HitDensity() );
            myTree.FillTree();
  
            std::string outputLabel="";
            if(_FRAMSAlgo.Score()>=fFRANSScoreCut){
                _EfficiencyCalculator.UpdateSelected(ev);
                outputLabel="Accepted";   
            }
            else{
                _EfficiencyCalculator.UpdateNotSelected(ev);
                outputLabel="Rejected";
            }
            std::string outputLabel2 = (_FRAMSAlgoPANDORA.Score()>=fFRANSScoreCut)? "Accepted":"Rejected";
            _FRAMSAlgo.Display(ev.Label()+"_FRANS_"+outputLabel);
            _FRAMSAlgoPANDORA.Display(ev.Label()+"_FRANSPANDORAVx_"+outputLabel2);
        }

    }
    
    // Print final status
    std::cout<<_EfficiencyCalculator;
    _EfficiencyCalculator.DrawHistograms();
        

    TFile outputFile((fAnaReultsPath+"/FRAMSTree.root").c_str(), "RECREATE");
    fTree->Write();
    delete fTree;
    outputFile.Close();

    return;
}




int main(int argc, char* argv[]){

    CommandLineParser parser(argc, argv);
    
    // Create a ROOT application object
    TApplication *myApp = new TApplication("myApp", &argc, argv);

    RunAlgoFRANS(parser);

    // Run the ROOT event loop
    myApp->Run();

    return 0;
}



/*_FRAMSAlgo.Fill(hitList, fVertexReco);
double scoreReco = _FRAMSAlgo.Score();
_FRAMSAlgo.Fill(hitList, fVertexTrue);
double scoreTrue = _FRAMSAlgo.Score();
_FRAMSAlgo.Fill(hitList, fVertexMine);
double scoreMine = _FRAMSAlgo.Score();
std::cout<<"JUJU "<<_TPCLinesAlgo.GetMainVertex().X()<<" "<<_TPCLinesAlgo.ShiftX()<<std::endl;
std::cout<<" FRANS Reco vertex (PANDORA) Score="<<scoreReco<<" Vertex="<<fVertexReco;            
std::cout<<" FRANS True vertex (PANDORA) "<<scoreTrue<<" Vertex="<<fVertexTrue;
std::cout<<" FRANS Using my fabolous vertex "<<scoreMine<<" Vertex="<<fVertexMine;*/