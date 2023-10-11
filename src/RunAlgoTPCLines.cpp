#include <string>
#include <iostream>

#include "TString.h"
#include <TApplication.h>

#include "CommandLineParser.h"
#include "SEventHandle.h"
#include "SParserConfig.h"
#include "STPCAnalyzerTreeReader.h"
#include "TPCSimpleEvents.h"
#include "TPCLinesParameters.h"
#include "TPCLinesAlgo.h"

void RunAlgoTPCLines(const CommandLineParser& parser)
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


    // Get input files
    std::vector<TString> fFilePaths = GetInputFileList(file_name, ext, directory_path);

    // Set batch mode
    if(Debug==0) gROOT->SetBatch(true);

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

    fPsetAnaView.Print();

    // Check if the directory exists, create it if not
    gSystem->Exec(("rm -rf "+fPsetAnaView.OutputPath).c_str());
    if (!gSystem->OpenDirectory(fPsetAnaView.OutputPath.c_str())) {    
        gSystem->Exec( ("mkdir "+fPsetAnaView.OutputPath).c_str());
        gSystem->Exec( ("mkdir "+fPsetAnaView.OutputPath+"/rootfiles").c_str());
    }  

    // output ROOT files with analysis results
    TFile* anaOutputFile = new TFile("LambdaAnaOutput.root", "RECREATE");

    // View to use
    std::string fView = fPsetAnaView.View;


    // ------------------------------------------------------------------ 
    // Define the program control variables
    int fNEv = n;
    int fEv = event;
    int fSubRun = sr;
    int fNEvSkip = nskip;
    int nEvents=0;


    // Define TPC LINES ALGORITHM
    TPCLinesAlgo _TPCLinesAlgo(fPsetAnaView);
    // Effiency status
    EfficiencyCalculator _EfficiencyCalculator;

    // TTree loop
    std::vector<int> fNVertex;
    int nEntries = 0;
    for (const auto& filepath : fFilePaths) {

        std::cout<<" ANALYZING THE FILE: "<<filepath<<std::endl;
        MyTPCTreeReader treeReader(filepath, "ana/AnaTPCTree");
        
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
            if(nhits<=3 || treeReader.recnuvU==-1){
                std::cout<<"   SKIPPED NHits or RecoVertex \n";
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
            bool filled = _TPCLinesAlgo.SetHitList(view, RecoVertexUVYT, VertexUVYT, 
                                    treeReader.hitsChannel,
                                    treeReader.hitsPeakTime,
                                    treeReader.hitsIntegral, 
                                    treeReader.hitsRMS,
                                    treeReader.hitsStartT, 
                                    treeReader.hitsEndT,
                                    treeReader.hitsChi2,
                                    "");
            
            // Analyze
            SEvent recoEvent;
            if(filled){
                _TPCLinesAlgo.AnaView(ev.Label());
                recoEvent = _TPCLinesAlgo.GetRecoEvent();
            }
            
            int nAngles = recoEvent.GetNAngles();
            // Update the efficiency calculator
            if(nAngles>0){
                _EfficiencyCalculator.UpdateSelected(ev);

                if(Debug==-12){
                    gROOT->SetBatch(false);
                    _TPCLinesAlgo.Display(ev.Label(), "Misselected");
                    gROOT->SetBatch(true);
                }

            }
            else
                _EfficiencyCalculator.UpdateNotSelected(ev);

            if(recoEvent.GetNOrigins()>0){
                _EfficiencyCalculator.UpdateHistograms(recoEvent);

                if(Debug==-13){

                    std::vector<SOrigin> origins = recoEvent.GetOrigins();

                    if(origins.size()==2){
                        int maxMult = std::max(origins[0].Multiplicity(), origins[1].Multiplicity());
                        int minMult = std::min(origins[0].Multiplicity(), origins[1].Multiplicity());
                        if(maxMult==2 && minMult==1){
                            gROOT->SetBatch(false);
                            _TPCLinesAlgo.Display(ev.Label(), "Misselected");
                            gROOT->SetBatch(true);
                        }
                    }
                }
            }
            std::cout<<_EfficiencyCalculator;

        }
    }

    // Origins Ana Results
    std::cout<<_EfficiencyCalculator;
    TDirectory *originsAnaDirectory = anaOutputFile->mkdir("originsAnaDirectory");
    TCanvas *cOriginsAna = new TCanvas("cOriginsAna", "cOriginsAna", 0, 0, 1400,900);
    originsAnaDirectory->cd();
    _EfficiencyCalculator.DrawHistograms(cOriginsAna);
    cOriginsAna->Write();


    anaOutputFile->Write();
    anaOutputFile->Close();

    return;

}

int main(int argc, char* argv[]){

    CommandLineParser parser(argc, argv);

    std::cout << "Debug: " << parser.getDebug() << std::endl;
    std::cout << "DebugMode: " << parser.getDebugMode() << std::endl;
    std::cout << "n: " << parser.getN() << std::endl;
    std::cout << "nskip: " << parser.getNskip() << std::endl;
    std::cout << "event: " << parser.getEvent() << std::endl;
    std::cout << "sr: " << parser.getSr() << std::endl;
    std::cout << "file_name: " << parser.getFileName() << std::endl;
    std::cout << "directory_path: " << parser.getDirectoryPath() << std::endl;
    std::cout << "ext: " << parser.getExtension() << std::endl;

    // Create a ROOT application object
    TApplication *myApp = new TApplication("myApp", &argc, argv);


    RunAlgoTPCLines(parser);


    // Run the ROOT event loop
    myApp->Run();    

    return 0;
}