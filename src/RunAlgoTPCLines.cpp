#include <string>
#include <iostream>

#include "TString.h"
#include <TApplication.h>
#include "TImage.h"

#include "CommandLineParser.h"
#include "SEventHandle.h"
#include "LambdaAnaTree.h"
#include "SParserConfig.h"
#include "STPCAnalyzerTreeReader.h"
#include "TPCSimpleEvents.h"
#include "TPCLinesParameters.h"
#include "TPCLinesAlgo.h"
#include "ChargeDensity.h"
#include "ChargeDensityPset.h"

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

    // Output ana tree
    TTree *fTreeAna = new TTree("LambdafAnaTreeHandle", "LambdafAnaTreeHandle");
    LambdaAnaTree fAnaTreeHandle(fTreeAna);

    // Get input files
    std::vector<TString> fFilePaths = GetInputFileList(file_name, ext, directory_path);

    // Set batch mode
    if(Debug==0) gROOT->SetBatch(true);

    // ---- FRAMS parameters ----------------------------------------
    double fFRANSScoreCut = 0.1;

    FRAMSPsetType fPsetFRANS = ReadFRANSPset( FindFile("chargedensityalg_config.fcl"), "ChargeDensityAlg:");
    fPsetFRANS.Verbose = Debug;
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


    fPsetAnaView.Print();
    fPsetFRANS.Print();

    // Check if the directory exists, create it if not
    gSystem->Exec(("rm -rf "+fPsetAnaView.OutputPath).c_str());
    if (!gSystem->OpenDirectory(fPsetAnaView.OutputPath.c_str())) {    
        gSystem->Exec( ("mkdir "+fPsetAnaView.OutputPath).c_str());
        gSystem->Exec( ("mkdir "+fPsetAnaView.OutputPath+"/rootfiles").c_str());
    }  

    // output ROOT files with analysis results
    TFile* anaOutputFile = new TFile("LambdaAnaOutput.root", "RECREATE");


    // ------------------------------------------------------------------ 
    // Define the program control variables
    int fNEv = n;
    int fEv = event;
    int fSubRun = sr;
    int fNEvSkip = nskip;
    int nEvents=0;



    // Define TPC LINES ALGORITHM
    TPCLinesAlgo _TPCLinesAlgo(fPsetAnaView);
    // Define FRANS LINES ALGORITHM
    ChargeDensity _FRAMSAlgo(fPsetFRANS);
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
            std::cout << "   Interaction mode: "<<treeReader.intMode<<" NLambda: "<<treeReader.intNLambda<<std::endl;
            
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


            int view = fPsetAnaView.View;
            if(fPsetAnaView.View==-1){ // use best view
                int bestView=_TPCLinesAlgo.GetBestView(treeReader.hitsView, treeReader.hitsChi2);
                std::cout<<"  Using best view: "<<bestView<<" TPC="<<std::to_string(TPC)<<std::endl;
                view = bestView;
            }
            
            // Get the hits in the view
            std::vector<SHit> hitList = GetHitsInView(view,
                                                        treeReader.hitsChannel,
                                                        treeReader.hitsPeakTime,
                                                        treeReader.hitsIntegral, 
                                                        treeReader.hitsRMS,
                                                        treeReader.hitsStartT, 
                                                        treeReader.hitsEndT,
                                                        treeReader.hitsView,
                                                        treeReader.hitsChi2);

            // Set the hits
            bool filled = _TPCLinesAlgo.SetHitList(view, RecoVertexUVYT, VertexUVYT, hitList);
            
            // Analyze
            SEvent recoEvent;
            if(filled){
                _TPCLinesAlgo.AnaView(ev.Label());
                recoEvent = _TPCLinesAlgo.GetRecoEvent();
            }

            // FRANS part
            std::vector<STriangle> angleList = recoEvent.GetAngleList();
            std::vector<SOrigin> associatedOrigins = recoEvent.GetAssociatedOrigins();
            int bestTriangleIx = -1;
            double bestFRANSScore = -1000;
            TCanvas *cDisplayPANDORA = new TCanvas( "FinalRecoFRANSPANDORA", "FinalRecoFRANS", 600, 0, 800, 1200);
            TCanvas *cDisplay = new TCanvas( "FinalRecoFRANS", "FinalRecoFRANS", 700, 0, 900, 1200);
            

            SVertex fVertexReco = SVertex( SPoint( RecoVertexUVYT[2], RecoVertexUVYT[3]), "");
            SVertex fVertexTrue = SVertex( SPoint( VertexUVYT[2], VertexUVYT[3]), "");

            _FRAMSAlgo.Fill(hitList, fVertexReco);
            _FRAMSAlgo.Display(cDisplayPANDORA);
            double FRANSScorePANDORA = _FRAMSAlgo.Score();

            for(size_t orix=0; orix<angleList.size(); orix++){
                SVertex fVertex(associatedOrigins[orix].GetPoint(), std::to_string(view));
                std::cout<<" LAMBDA CANDIDATE "<<fVertex<<std::endl;

                std::cout<<" COVERED AREA: "<<angleList[orix].ComputeCoveredArea();

                SVertex fVertexMine = SVertex( SPoint( fVertex.X()+_TPCLinesAlgo.ShiftX(),
                                                fVertex.Y()+_TPCLinesAlgo.ShiftY())
                                            , "");

                                        

                std::cout<<" VertexTrue "<<fVertexTrue;
                std::cout<<" VertexMine "<<fVertexMine;
                
                _FRAMSAlgo.Fill(hitList, fVertexMine);

                double score = _FRAMSAlgo.Score();
                if(score>bestFRANSScore){
                    bestFRANSScore = score;
                    bestTriangleIx = orix;
                    _FRAMSAlgo.Display(cDisplay);
                }
            }
            
            int nAngles = recoEvent.GetNAngles();

            int nOrigins = recoEvent.GetNOrigins();
            int nOriginsMultGT3;
            if(bestTriangleIx!=-1) nOriginsMultGT3 = recoEvent.GetNOriginsMultGt(3, angleList[bestTriangleIx].GetTrack1().GetId(), angleList[bestTriangleIx].GetTrack2().GetId());
            else nOriginsMultGT3 = recoEvent.GetNOriginsMultGt(3);
            nOriginsMultGT3 = recoEvent.GetNOriginsMultGt(3);

            // Print the track associations
            recoEvent.PrintTrackConnections();

    
            SEvent notAssociatedRecoEvent ({}, {}, {}, {}, 0, {});
            if(bestTriangleIx!=-1){
                notAssociatedRecoEvent = recoEvent.UnassociatedOrigins(angleList[bestTriangleIx]);
            }

            std::cout<<"  - Not associated reco event: \n";
            std::cout<<" NOrigins: "<<notAssociatedRecoEvent.GetNOrigins()<<std::endl;
            std::cout<<" NOrigins mult 2: "<<notAssociatedRecoEvent.GetNOriginsMult(2)<<std::endl;
            std::cout<<" NOrigins mult 1: "<<notAssociatedRecoEvent.GetNOriginsMult(1)<<std::endl;
            std::cout<<" NOrigins mult GT 3: "<<notAssociatedRecoEvent.GetNOriginsMultGt(3)<<std::endl;


            int nDirtHitsInTriangle = 0;
            double nFractionDirtHitsInTriangle = 0;
            int nDirtHitsInTriangleWires = 0;
            double nFractionDirtHitsInTriangleWires = 0;
            if(bestTriangleIx!=-1){
                recoEvent.FreeHitsAroundTriangle(angleList[bestTriangleIx],
                                                nDirtHitsInTriangle,
                                                nFractionDirtHitsInTriangle,
                                                nDirtHitsInTriangleWires,
                                                nFractionDirtHitsInTriangleWires);
            }

            
            int nFreeHits = 0;
            int fNUnassociatedHits = 0;
            if(bestTriangleIx!=-1){
                recoEvent.GetUnassociatedHits(angleList[bestTriangleIx], nFreeHits, fNUnassociatedHits);
            }

            //bool accepted = nAngles>0 && bestFRANSScore>fFRANSScoreCut;
            bool accepted = nAngles>0 && bestFRANSScore>fFRANSScoreCut && nOrigins<=6 && nOriginsMultGT3==0;
    
            
            // Update the efficiency calculator
            if(accepted){
                _EfficiencyCalculator.UpdateSelected(ev);

                if(Debug==-12){
                    gROOT->SetBatch(false);
                    _TPCLinesAlgo.Display( "Misselected"+ev.Label());
                    gROOT->SetBatch(true);
                }

            }
            else
                _EfficiencyCalculator.UpdateNotSelected(ev);

            
            std::string outNamePreffix = accepted? "Accepted":"Rejected";
            std::string outputLabel = "FinalReco" + outNamePreffix + ev.Label() + "_" + std::to_string(_EfficiencyCalculator.NEvents()); 
            
            if(recoEvent.GetNOrigins()>0){
                _EfficiencyCalculator.UpdateHistograms(recoEvent);

                if(Debug==-13){

                    std::vector<SOrigin> origins = recoEvent.GetOrigins();

                    if(origins.size()==2){
                        int maxMult = std::max(origins[0].Multiplicity(), origins[1].Multiplicity());
                        int minMult = std::min(origins[0].Multiplicity(), origins[1].Multiplicity());
                        if(maxMult==2 && minMult==1){
                            gROOT->SetBatch(false);
                            _TPCLinesAlgo.Display("Misselected"+ev.Label());
                            gROOT->SetBatch(true);
                        }
                    }
                }
            }
            std::cout<<_EfficiencyCalculator;

            fAnaTreeHandle.ResetVars();
            fAnaTreeHandle.fEventID = treeReader.eventID;
            fAnaTreeHandle.fSubrunID = treeReader.subrunID;
            fAnaTreeHandle.fRunID = treeReader.runID;
            fAnaTreeHandle.fSliceID = 1;

            fAnaTreeHandle.fIntMode = treeReader.intMode;
            fAnaTreeHandle.fIntNLambda = treeReader.intNLambda;

            fAnaTreeHandle.fFRANSScorePANDORA = FRANSScorePANDORA;

            // Number of origins variables
            fAnaTreeHandle.fNOrigins = recoEvent.GetNOrigins();
            fAnaTreeHandle.fNOriginsMult1 = recoEvent.GetNOriginsMult(1);
            fAnaTreeHandle.fNOriginsMult2 = recoEvent.GetNOriginsMult(2);
            fAnaTreeHandle.fNOriginsMultGT3 = recoEvent.GetNOriginsMultGt(3);
            fAnaTreeHandle.fNOriginsPairOneTwo = recoEvent.GetNOriginsMult(2) * recoEvent.GetNOriginsMult(1);

            // Angle variables
            fAnaTreeHandle.fNAngles = recoEvent.GetNAngles(); 
            if(bestTriangleIx!=-1){
                STriangle bestTriangle = angleList[bestTriangleIx];

                fAnaTreeHandle.fAngleFRANSScore = bestFRANSScore;
                fAnaTreeHandle.fAngleGap = bestTriangle.GetGap();
                fAnaTreeHandle.fAngleDecayContainedDiff = bestTriangle.GetDecayAngleDifference();
                fAnaTreeHandle.fAngleNHits = bestTriangle.GetNHitsTriangle();
                fAnaTreeHandle.fAngleNHitsTrack1 = bestTriangle.GetNHitsTrack1();
                fAnaTreeHandle.fAngleNHitsTrack2 = bestTriangle.GetNHitsTrack2();
                fAnaTreeHandle.fAngleNHitsMainTrack = bestTriangle.GetNHitsMainTrack();
                fAnaTreeHandle.fAngleLengthTrack1 = bestTriangle.GetLengthTrack1();
                fAnaTreeHandle.fAngleLengthTrack2 = bestTriangle.GetLengthTrack2();
                fAnaTreeHandle.fAngleLengthMainTrack = bestTriangle.GetLengthMainTrack();

                std::cout<<"  - Best angle: "<<bestTriangleIx<<" FRANS: "<<bestFRANSScore<<" Gap: "<<bestTriangle.GetGap()<<" DecayContainedDiff: "<<bestTriangle.GetDecayAngleDifference()<<std::endl;
                std::cout<<"  - NHits: "<<bestTriangle.GetNHitsTriangle()<<" NHitsTrack1: "<<bestTriangle.GetNHitsTrack1()<<" NHitsTrack2: "<<bestTriangle.GetNHitsTrack2()<<" NHitsMainTrack: "<<bestTriangle.GetNHitsMainTrack()<<std::endl;
                std::cout<<"  - LengthTrack1: "<<bestTriangle.GetLengthTrack1()<<" LengthTrack2: "<<bestTriangle.GetLengthTrack2()<<" LengthMainTrack: "<<bestTriangle.GetLengthMainTrack()<<std::endl;
                // cout minimum hits
                std::cout<<"  - Minimum hits: "<<std::min(bestTriangle.GetNHitsTrack1(), bestTriangle.GetNHitsTrack2())<<std::endl;
                //cout opening angle
                std::cout<<"  - Opening angle: "<<bestTriangle.GetOpeningAngle()<<std::endl;
            }

            int associatedHits = recoEvent.NHits();
            int totalHits = _TPCLinesAlgo.GetNInputHits();
            int unassociatedHits = totalHits - associatedHits;
            fAnaTreeHandle.fNUnassociatedHits = unassociatedHits;
            std::cout<<"  - Associated hits: "<<associatedHits<<" Total hits: "<<totalHits<<" Unassociated hits: "<<unassociatedHits<<std::endl;

            fAnaTreeHandle.fTruthIsFiducial = true;
            fAnaTreeHandle.fRecoIsFiducial = true;
            
            //fAnaTreeHandle.FillTree();

            TCanvas *cTPCDisplay = new TCanvas( ("FinalReco"+ev.Label()).c_str(), "FinalRecoTPCLines", 0, 0, 1000, 800);
            _TPCLinesAlgo.Display("", cTPCDisplay);
            cTPCDisplay->SaveAs( (fPsetAnaView.OutputPath+"/"+outputLabel+".pdf").c_str() );
            delete cTPCDisplay;
            
            outputLabel+="FRANS";
            if(bestFRANSScore!=-1000){
                cDisplay->SaveAs( (fPsetAnaView.OutputPath+"/"+outputLabel+".pdf").c_str() );
            }

            outputLabel+="FRANSPANDORA";
            cDisplayPANDORA->SaveAs( (fPsetAnaView.OutputPath+"/"+outputLabel+".pdf").c_str() );
            delete cDisplay;

        }
    }

    // Origins Ana Results
    std::cout<<_EfficiencyCalculator;
    TDirectory *originsAnaDirectory = anaOutputFile->mkdir("originsAnaDirectory");
    TCanvas *cOriginsAna = new TCanvas("cOriginsAna", "cOriginsAna", 0, 0, 1400,900);
    originsAnaDirectory->cd();
    _EfficiencyCalculator.DrawHistograms(cOriginsAna);
    cOriginsAna->Write();
    TDirectory *fAnaTreeHandleDirectory = anaOutputFile->mkdir("originsAna");
    fAnaTreeHandleDirectory->cd();
    fAnaTreeHandle.WriteTree();
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