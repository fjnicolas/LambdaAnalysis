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

#include "TImage.h"


std::vector<SHit> GetFRANSHitsView(
    int view,
    int tpc,
    std::vector<int> *_X,
    std::vector<double> *_Y,
    std::vector<double> *_Int,
    std::vector<double> *_Wi,
    std::vector<double> *_ST,
    std::vector<double> *_ET,
    std::vector<int> *_View,
    std::vector<double> *_Chi2)
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

    // set variables
    std::vector<SHit> hitList;
    int nTotalHits = _X->size();

    // loop over the hits
    for (int i = 0; i < nTotalHits; i++) {

        // filter channels for the view        
        if ( _View->at(i)==view ) {
            SHit hit(-1, _X->at(i), _Y->at(i), _Wi->at(i), _Int->at(i), _ST->at(i), _ET->at(i), _Chi2->at(i));
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
    FRAMSPsetType fPsetFRANS = ReadFRANSPset( FindFile("chargedensityalg_config.fcl"), "ChargeDensityAlg:");
    fPsetFRANS.Verbose = Debug;
    std::cout<<"  njvisfnvioanvownvr "<<fPsetFRANS.TMVAFilename<<" "<<fPsetFRANS.OutputPath<<std::endl;
    fPsetFRANS.TMVAFilename = FindFile(fPsetFRANS.TMVAFilename);
    
    // ---- TPCLines parameters ----------------------------------------
    TrackFinderAlgorithmPsetType fPsetTrackFinder = ReadTrackFinderAlgorithmPset( FindFile("trackfinderalg_config.fcl"), "TrackFinderAlg:");
    fPsetTrackFinder.Verbose = Debug;
    HoughAlgorithmPsetType fPsetHough = ReadHoughAlgorithmPset( FindFile("houghalg_config.fcl"), "HoughAlg:");
    fPsetHough.Verbose = Debug;
    VertexFinderAlgorithmPsetType fPsetVertexFinder = ReadVertexFinderAlgorithmPset( FindFile("vertexfinderalg_config.fcl"), "VertexFinderAlg:");
    fPsetVertexFinder.Verbose = Debug;
    TPCLinesAlgoPsetType fPsetAnaView = ReadTPCLinesAlgoPset( FindFile("tpclinesalg_config.fcl"), "TPCLinesAlg:");
    fPsetAnaView.Verbose = Debug;
    fPsetAnaView.DebugMode = DebugMode;
    fPsetAnaView.HoughAlgorithmPset = fPsetHough;
    fPsetAnaView.TrackFinderAlgorithmPset = fPsetTrackFinder;
    fPsetAnaView.VertexFinderAlgorithmPset = fPsetVertexFinder;

    // Define the program control variables
    int fNEv = n;
    int fEv = event;
    int fSubRun = sr;
    int fNEvSkip = nskip;
    int nEvents=0;


    // output ROOT files with analysis results
    TFile* anaOutputFile = new TFile("LambdaAnaOutput.root", "RECREATE");    
    TDirectory *fransPlotsDirectory = anaOutputFile->mkdir("FRANSPlots");

    // Output directory for analysis results
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
    EfficiencyCalculator _EfficiencyCalculator;

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
            int view = fPsetAnaView.View;
            if(fPsetAnaView.View==-1){ // use best view
                int bestView=_TPCLinesAlgo.GetBestView(treeReader.hitsView, treeReader.hitsChi2);
                std::cout<<"  Using best view: "<<bestView<<" TPC="<<std::to_string(TPC)<<std::endl;
                view = bestView;
            }
            
            // Set the hits
            std::vector<SHit> hitList = GetFRANSHitsView(view,
                                                        TPC,
                                                        treeReader.hitsChannel,
                                                        treeReader.hitsPeakTime,
                                                        treeReader.hitsIntegral, 
                                                        treeReader.hitsRMS,
                                                        treeReader.hitsStartT, 
                                                        treeReader.hitsEndT,
                                                        treeReader.hitsView,
                                                        treeReader.hitsChi2);
                                        
            // Set the vertex
            // true
            double vertexXTrue = VertexUVYT[2];
            if ( view == 0 ) vertexXTrue = VertexUVYT[0];
            if ( view == 1 ) vertexXTrue = VertexUVYT[1];
            double vertexYTrue = VertexUVYT[3];
            SVertex fVertexTrue(SPoint(vertexXTrue, vertexYTrue), std::to_string(view));
            
            // reco
            double vertexXReco = RecoVertexUVYT[2];
            if ( view == 0 ) vertexXReco = RecoVertexUVYT[0];
            if ( view == 1 ) vertexXReco = RecoVertexUVYT[1];
            double vertexYReco = RecoVertexUVYT[3];
            SVertex fVertexReco(SPoint(vertexXReco, vertexYReco), std::to_string(view));

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
                                        treeReader.hitsView,
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
  


            _FRAMSAlgo.Score()>=fFRANSScoreCut ? _EfficiencyCalculator.UpdateSelected(ev):_EfficiencyCalculator.UpdateNotSelected(ev);
            std::string outputLabel = (_FRAMSAlgo.Score()>=fFRANSScoreCut)? "plot_Accepted":"plot_Rejected";
            std::string outputLabel2 = (_FRAMSAlgoPANDORA.Score()>=fFRANSScoreCut)? "Accepted":"Rejected";

                        
            TCanvas *cDisplay = new TCanvas( (outputLabel+"_"+ev.Label()+"_vw"+std::to_string(view)).c_str(), outputLabel.c_str(), 600, 0, 800, 1200);
            _FRAMSAlgo.Display(cDisplay);

            // Save TCanvas and pdf
            TImage *img = TImage::Create();
            img->FromPad(cDisplay);
            fransPlotsDirectory->cd();
            img->Write( ("image/"+outputLabel+"_"+ev.Label()+"_vw"+std::to_string(view)+".pdf").c_str() );
            cDisplay->Write();

            delete cDisplay;
            delete img;
        }

    }
    
    // Origins Ana Results
    std::cout<<_EfficiencyCalculator;
    TDirectory *originsAnaDirectory = anaOutputFile->mkdir("originsAnaDirectory");
    TCanvas *cOriginsAna = new TCanvas("cOriginsAna", "cOriginsAna", 0, 0, 1400,900);
    originsAnaDirectory->cd();
    _EfficiencyCalculator.DrawHistograms(cOriginsAna);
    cOriginsAna->Write();
        
    // Write the FRANS TTree
    anaOutputFile->cd();
    fTree->Write();
    delete fTree;

    anaOutputFile->Write();
    anaOutputFile->Close();

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