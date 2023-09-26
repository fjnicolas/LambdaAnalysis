// Objects
#include "src/FRAMS/ChargeDensity/ChargeDensity.cc"


void RunAlgoFRANS(int Debug=0, int DebugMode=-1, int n=1e6, int nskip=-1, int event=-1, int sr=-1, std::string file_name="", const char *directory_path=".", const char *ext=".root")
{

    // ----------- ALGORITHM PARAMETERS --------------------------------- 
    // View to use
    std::string fView="C";

    // ---- FRAMS parameters ----------------------------------------

    FRAMSPsetType fPsetFRAMS(
        true,                          // ApplyRawSmoothing
        false,                         // ApplySmoothing
        false,                         // ApplyCumulativeSmoothing
        4,                             // NDriftPack
        1,                             // NWirePack
        0.3f,                          // ExpoAvSmoothPar
        1,                             // UnAvNeighbours
        0.8,                           // CumulativeCut
        3,                             // SlidingWindowN
        3,                             // NSamplesBeginSlope
        70,                            // MaxRadius
        false,                         // DebugMode
        true,                          // CalculateScore
        "FRAMSSelectionTMVA_BDT.weights.xml"  // TMVAFilename
    );
    
    // ------------------------------------------------------------------ 


    // Get the candidate Files
    std::vector<TString> fFilePaths;
    TSystemDirectory dir(".", ".");
    TList *files = dir.GetListOfFiles();
    TString targetFileName(file_name);
    std::cout<<" Target file name "<<targetFileName<<std::endl;
    gSystem->Exec("ls");
    if (files){
        TSystemFile *file;
        TString fname;
        TIter next(files);
        while ((file=(TSystemFile*)next())) {
            fname = file->GetName();
            
            if (!file->IsDirectory() && fname.EndsWith(ext)){
                if(fname.Contains(targetFileName)){
                    fFilePaths.push_back(fname);
                }
                std::cout << fname << " Target:" << targetFileName <<std::endl;
                std::cout<<" Contains: "<<fname.Contains(targetFileName)<<std::endl;
                
                
            }
        }
    }

    // Define the program control variables
    int fNEv = n;
    int fEv = event;
    int fSubRun = sr;
    int fNEvSkip = nskip;
    int nEvents=0;
    // Output paths for displays
    std::string fAppDisplayPath = "plots/";
    if(DebugMode==0) fAppDisplayPath = "plotsbg";
    else if(DebugMode==1) fAppDisplayPath = "plotssignal";
    // SBND readout parameters
    int fStampTime = -200;
    double fSamplingFrequency = 500;
    double fSamplingTime = 0.5;

    // Define TPC LINES ALGORITHM
    FRAMSAlgo _FRAMSAlgo(fPsetAnaView);
    // Effiency status
    EfficiencyCalculator _EfficiencyCalculator;

    // TTree loop
    std::vector<int> fNVertex;
    int nEntries = 0;
    for (const auto& filepath : fFilePaths) {

        std::cout<<" ANALYZING THE FILE: "<<filepath<<std::endl;

        TFile* f= new TFile(filepath);
	    TTree* tree = (TTree*)f->Get("ana/AnaTPCTree");
        
        // event ID
        int eventID;
        int subrunID;
        int runID;
        tree->SetBranchAddress("EventID",&eventID);
        tree->SetBranchAddress("RunID",&runID);
        tree->SetBranchAddress("SubRunID",&subrunID);
        // true Vertex
        double nuvE, nuvT, nuvX, nuvY, nuvZ;
        int nuvU, nuvV, nuvC, nuvTimeTick;
        int TPC;
        double nuvDriftTime;
        tree->SetBranchAddress("TrueVEnergy",&nuvE);
        tree->SetBranchAddress("TrueVt",&nuvT);
        tree->SetBranchAddress("TrueVx",&nuvX);
        tree->SetBranchAddress("TrueVy",&nuvY);
        tree->SetBranchAddress("TrueVz",&nuvZ);
        tree->SetBranchAddress("TrueVU",&nuvU);
        tree->SetBranchAddress("TrueVV",&nuvV);
        tree->SetBranchAddress("TrueVC",&nuvC);
        tree->SetBranchAddress("TrueVTimeTick",&nuvTimeTick);
        // reco Vertex
        double recnuvX = -1, recnuvY = -1, recnuvZ = -1;
        int recnuvU = -1, recnuvV = -1, recnuvC = -1, recnuvTimeTick = -1;
        tree->SetBranchAddress("RecoVx",&recnuvX);
        tree->SetBranchAddress("RecoVy",&recnuvY);
        tree->SetBranchAddress("RecoVz",&recnuvZ);
        tree->SetBranchAddress("RecoVU",&recnuvU);
        tree->SetBranchAddress("RecoVV",&recnuvV);
        tree->SetBranchAddress("RecoVC",&recnuvC);
        tree->SetBranchAddress("RecoVTimeTick",&recnuvTimeTick);

        std::vector<double> * hitsIntegral=new std::vector<double>;
        std::vector<double> * hitsPeakTime=new std::vector<double>;
        std::vector<int> * hitsChannel=new std::vector<int>;
        std::vector<double> * hitsRMS=new std::vector<double>;
        //std::vector<double> * hitsWidth=new std::vector<double>;
        std::vector<double> * hitsStartT=new std::vector<double>;
        std::vector<double> * hitsEndT=new std::vector<double>;
        std::vector<double> * hitsChi2=new std::vector<double>;
        std::vector<double> * hitsNDF=new std::vector<double>;
        std::vector<int> * hitsClusterID=new std::vector<int>;        
        tree->SetBranchAddress("HitsIntegral",&hitsIntegral);
        tree->SetBranchAddress("HitsPeakTime",&hitsPeakTime);
        tree->SetBranchAddress("HitsChannel",&hitsChannel);
        tree->SetBranchAddress("HitsRMS",&hitsRMS);
        //tree->SetBranchAddress("HitsWidth",&hitsWidth);
        tree->SetBranchAddress("HitsStartT",&hitsStartT);
        tree->SetBranchAddress("HitsEndT",&hitsEndT);
        tree->SetBranchAddress("HitsChi2",&hitsChi2);
        tree->SetBranchAddress("HitsNDF",&hitsNDF);
        tree->SetBranchAddress("HitsClusterID",&hitsClusterID);


        for (int entry = 0; entry < tree->GetEntries(); entry++) {
            std::cout<<"Entry "<<entry<<std::endl;
            // get the entry
            tree->GetEntry(entry);

            // check program control variables
            SEventId ev(runID, subrunID, eventID);
            if (fNEv > 0 && nEvents >= fNEv) continue;
            if (eventID != fEv && fEv != -1) continue;
            if (subrunID != fSubRun && fSubRun != -1) continue;
            nEntries++;
            if (nEntries <= fNEvSkip) continue;
            
            nEvents++;
            std::cout << "\n\n ************** Analyzing: " << ev;
            
            // True vertex
            nuvDriftTime = nuvTimeTick * fSamplingTime + fStampTime;
            std::vector<double> VertexXYZ = {nuvX, nuvY, nuvZ};
            std::vector<int> VertexUVYT = {nuvU, nuvV, nuvC, nuvTimeTick};
            std::cout << "   - True vertex (X, Y, Z, T) " << VertexXYZ[0] << " " << VertexXYZ[1] << " " << VertexXYZ[2] << " " << nuvT << " (U, V, C, TT): " << VertexUVYT[0] << " " << VertexUVYT[1] << " " << VertexUVYT[2] << " " << VertexUVYT[3] << std::endl;

            // Reco vertex
            double recnuvDriftTime = recnuvTimeTick * fSamplingTime + fStampTime;
            std::vector<double> RecoVertexXYZ = {recnuvX, recnuvY, recnuvZ};
            std::vector<int> RecoVertexUVYT = {recnuvU, recnuvV, recnuvC, recnuvTimeTick};
            std::cout << "  - Reco vertex (X, Y, Z) " << RecoVertexXYZ[0] << " " << RecoVertexXYZ[1] << " " << RecoVertexXYZ[2] << " (U, V, C, TT): " << RecoVertexUVYT[0] << " " << RecoVertexUVYT[1] << " " << RecoVertexUVYT[2] << " "<< RecoVertexUVYT[3]  << std::endl;

            // Set TPC using the drift coordinate of the reco vertex
            //TPC = (nuvX > 0) ? 1 : 0;
            TPC = (RecoVertexXYZ[0] > 0) ? 1 : 0;
            
            // We need a minimum number of hits to run the track finder
            size_t nhits = hitsChannel->size();
            std::cout << "  - NHits: " << nhits << std::endl;
            if(nhits<=3 || recnuvU==-1){
                std::cout<<"   SKIPPED NHits or RecoVertex \n";
                _EfficiencyCalculator.UpdateSkipped(ev);
                continue;
            }

            // Assing the view and TPC
            std::string view = fView+std::to_string(TPC);

            // Set the hits
            //_FRAMSAlgo.Fill(view)
        }
    }
    
    // Print final status
    std::cout<<_EfficiencyCalculator;
    _EfficiencyCalculator.DrawHistograms();


    return 0;

}


