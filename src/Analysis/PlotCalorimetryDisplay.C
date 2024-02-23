#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TString.h>
#include <TMath.h>

//---------  Load function
void PlotCalorimetryDisplay(){
    std::cout<<"LambdaAnalysis function loaded"<<std::endl;
    return;
}

//---------  Color list
std::vector<int> fColorList = {kRed, kBlue, kGreen, kMagenta, kCyan, kYellow, kOrange, kViolet, kTeal, kAzure, kGray, kPink, kSpring, kWhite, kBlack};

//--------- Define a function to add points to the graph and create a TGraph2D
void addPointsAndCreateGraph(TGraph2D *graph, std::vector<double> x, std::vector<double> y, std::vector<double> z){

    for(size_t k = 0; k < x.size(); k++) {
        graph->SetPoint(k, x.at(k), y.at(k), z.at(k));
    }

    return;
}


// Define a function to draw the graph
void drawGraph(TGraph2D *graph, const char* title, Color_t color){
    graph->SetMarkerColor(color);
    graph->SetMarkerStyle(20);
    graph->SetTitle(title);
    graph->Draw("P same");
}


//---------  Main function
void RunCalorimetryDisplay(std::string fInputFileName="", int event=-1, std::string fTreeDirName = "lambdaPidPandoraAna/", std::string fTreeName = "CalorimetryTree")
{

    //--------- Input TTree
    TFile *fFile = new TFile(fInputFileName.c_str(),"READ");
    TTree *fTree = (TTree *)fFile->Get((fTreeDirName+fTreeName).c_str());

    //--------- Binning
    Int_t nBins = 200;
    Double_t xMin = 0;
    Double_t xMax = 20;
    Double_t yMin = 0;
    Double_t yMax = 40;
    Double_t yMaxCharge = 10000;;
    // Stats to 0
    gStyle->SetOptStat(0);
    
    //--------- Set Branches 
    Int_t RunID;
    Int_t SubRunID;
    Int_t EventID;
    fTree->SetBranchAddress("RunID", &RunID);
    fTree->SetBranchAddress("SubRunID", &SubRunID);
    fTree->SetBranchAddress("EventID", &EventID);
    Bool_t RecoStatus;
    fTree->SetBranchAddress("RecoStatus", &RecoStatus);
    Double_t InvariantMass, KELI, KEHI, OpeningAngle;
    fTree->SetBranchAddress("InvariantMass", &InvariantMass);
    fTree->SetBranchAddress("KELI", &KELI);
    fTree->SetBranchAddress("KEHI", &KEHI);
    fTree->SetBranchAddress("OpeningAngle", &OpeningAngle);
    // True variables
    Double_t ProtonKE, PionKE, DecayCosOpeningAngle;
    fTree->SetBranchAddress("ProtonKE", &ProtonKE);
    fTree->SetBranchAddress("PionKE", &PionKE);
    fTree->SetBranchAddress("DecayCosOpeningAngle", &DecayCosOpeningAngle);
    // Space points
    std::vector<std::vector<double>> *SpacePointsX = new std::vector<std::vector<double>>;
    std::vector<std::vector<double>> *SpacePointsY = new std::vector<std::vector<double>>;
    std::vector<std::vector<double>> *SpacePointsZ = new std::vector<std::vector<double>>;
    fTree->SetBranchAddress("SpacePointsX", &SpacePointsX);
    fTree->SetBranchAddress("SpacePointsY", &SpacePointsY);
    fTree->SetBranchAddress("SpacePointsZ", &SpacePointsZ);
    // Residual range LI
    std::vector<double> *ResidualRangeLI = new std::vector<double>;
    std::vector<double> *DepositedEnergyLI = new std::vector<double>;
    std::vector<double> *DepositedChargeLI = new std::vector<double>;
    fTree->SetBranchAddress("ResidualRangeLI", &ResidualRangeLI);
    fTree->SetBranchAddress("DepositedEnergyLI", &DepositedEnergyLI);
    fTree->SetBranchAddress("DepositedChargeLI", &DepositedChargeLI);
    // Residual range HI
    std::vector<double> *ResidualRangeHI = new std::vector<double>;
    std::vector<double> *DepositedEnergyHI = new std::vector<double>;
    std::vector<double> *DepositedChargeHI = new std::vector<double>;
    fTree->SetBranchAddress("ResidualRangeHI", &ResidualRangeHI);
    fTree->SetBranchAddress("DepositedEnergyHI", &DepositedEnergyHI);
    fTree->SetBranchAddress("DepositedChargeHI", &DepositedChargeHI);

    
    // deDx (all tracks)
    std::vector<std::vector<double>> *ResidualRange = new std::vector<std::vector<double>>;
    std::vector<std::vector<double>> *DepositedEnergy = new std::vector<std::vector<double>>;
    std::vector<std::vector<double>> *DepositedCharge = new std::vector<std::vector<double>>;
    std::vector<int> *BestPIDAll = new std::vector<int>;
    fTree->SetBranchAddress("ResidualRangeAll", &ResidualRange);
    fTree->SetBranchAddress("DepositedEnergyAll", &DepositedEnergy);
    fTree->SetBranchAddress("DepositedChargeAll", &DepositedCharge);
    fTree->SetBranchAddress("BestPIDAll", &BestPIDAll);


    //--------- Tree loop
    for(int i=0; i<fTree->GetEntries(); i++){
        
        fTree->GetEntry(i);
        std::cout<<i<<" RunID: "<<RunID<<", SubRunID: "<<SubRunID<<", EventID: "<<EventID<<std::endl;

        if(event!=-1 && EventID!=event) continue;
        if( event==-2 && RecoStatus!=1) continue;

        std::cout<<"InvariantMass: "<<InvariantMass<<std::endl;
        std::cout<<"KELI: "<<KELI<<" KEPion: "<<1000*PionKE<<std::endl;
        std::cout<<"KEHI: "<<KEHI<<" KEProton: "<<1000*ProtonKE<<std::endl;
        std::cout<<"CosOpeningAngle: "<<OpeningAngle<<" DecayCosOpeningAngle: "<<std::acos(DecayCosOpeningAngle)*TMath::RadToDeg()<<std::endl;


        //-------- Canvas
        std::string canvasTitle = "RunID: "+std::to_string(RunID)+", SubRunID: "+std::to_string(SubRunID)+", EventID: "+std::to_string(EventID);
        TCanvas *c1 = new TCanvas("c1",canvasTitle.c_str(),800,600);
        TPad *pad1 = new TPad("pad1","This is pad1",0.02,0.52,0.48,0.98);
        TPad *pad2 = new TPad("pad2","This is pad2",0.02,0.02,0.48,0.48);
        TPad *pad3 = new TPad("pad3","This is pad3",0.52,0.52,0.98,0.98);
        TPad *pad4 = new TPad("pad4","This is pad4",0.52,0.02,0.98,0.48);
        pad1->Draw();
        pad2->Draw();
        pad3->Draw();
        pad4->Draw();


        //-------- Plot dEdx from TTree
        std::cout<<" ResidualRangeLI Size: "<<ResidualRangeLI->size()<<std::endl;
        if(ResidualRangeLI->size()>0 && ResidualRangeHI->size()>0){
            pad1->cd();
            TH2F *h_dEdxLI = new TH2F("h_dEdxLI","h_dEdxLI",nBins,xMin,xMax, nBins,yMin,yMax);
            fTree->Draw("DepositedEnergyLI:ResidualRangeLI>>h_dEdxLI","","colz", 1, i);
            TH2F *h_dEdxHI = new TH2F("h_dEdxHI","h_dEdxHI",nBins,xMin,xMax, nBins,yMin,yMax);
            fTree->Draw("DepositedEnergyHI:ResidualRangeHI>>h_dEdxHI","","colz", 1, i);

            // Colors and markers
            h_dEdxLI->SetMarkerColor(kRed);
            h_dEdxHI->SetMarkerColor(kBlue);
            h_dEdxLI->SetMarkerStyle(20);
            h_dEdxHI->SetMarkerStyle(20);
            
            //Draw histograms
            h_dEdxLI->GetYaxis()->SetTitle("dE/dx [MeV]");
            h_dEdxLI->GetXaxis()->SetTitle("Residual Range [cm]");
            h_dEdxLI->Draw("");
            h_dEdxHI->Draw("same");
            
            // Legend
            TLegend *leg = new TLegend(0.1,0.7,0.3,0.9);
            leg->AddEntry(h_dEdxLI,"Low Ionization","p");
            leg->AddEntry(h_dEdxHI,"High Ionization","p");
            leg->Draw("same");
            pad1->Update();

            //-------- Plot dQdx from TTree
            pad2->cd();
            TH2F *h_dQdxLI = new TH2F("h_dQdxLI","h_dQdxLI",nBins,xMin,xMax, nBins,yMin,yMaxCharge);
            fTree->Draw("DepositedChargeLI:ResidualRangeLI>>h_dQdxLI","","colz", 1, i);
            TH2F *h_dQdxHI = new TH2F("h_dQdxHI","h_dQdxHI",nBins,xMin,xMax, nBins,yMin,yMaxCharge);
            fTree->Draw("DepositedChargeHI:ResidualRangeHI>>h_dQdxHI","","colz", 1, i);

            // Colors and markers
            h_dQdxLI->SetMarkerColor(kRed);
            h_dQdxHI->SetMarkerColor(kBlue);
            h_dQdxLI->SetMarkerStyle(20);
            h_dQdxHI->SetMarkerStyle(20);

            //Draw histograms
            h_dQdxLI->GetYaxis()->SetTitle("dQ/dx [ADC/cm]");
            h_dQdxLI->GetXaxis()->SetTitle("Residual Range [cm]");
            h_dQdxLI->Draw("");
            h_dQdxHI->Draw("same");

            // Legend
            leg->Draw("same");

        }

        //-------- Plot dEdx (all tracks)
        pad3->cd();
        std::vector<TH1F*> h_dEdxAll;
        for(int j=0; j<ResidualRange->size(); j++){
            std::string histName = "h_dEdx_"+std::to_string(j);
            TH1F *h_dEdx = new TH1F(histName.c_str(),histName.c_str(),nBins,xMin,xMax);
            for(int k=0; k<ResidualRange->at(j).size(); k++){
                h_dEdx->Fill(ResidualRange->at(j).at(k), DepositedEnergy->at(j).at(k));
            }
            h_dEdx->SetMarkerColor(fColorList.at(j));
            h_dEdx->SetMarkerStyle(20);
            h_dEdxAll.push_back(h_dEdx);
        }
        
        if(h_dEdxAll.size()>0){
            h_dEdxAll.at(0)->Draw("p");
            for(int j=1; j<h_dEdxAll.size(); j++){
                h_dEdxAll.at(j)->Draw("p same");
            }

            TLegend *leg = new TLegend(0.1,0.7,0.3,0.9);
            for(int j=0; j<h_dEdxAll.size(); j++){
                std::string pidLabel = "Track "+std::to_string(j)+" PDG="+std::to_string(BestPIDAll->at(j));
                leg->AddEntry(h_dEdxAll.at(j), pidLabel.c_str(),"p");
            }
            leg->Draw("same");
            h_dEdxAll.at(0)->GetYaxis()->SetTitle("dE/dx [MeV]");
            h_dEdxAll.at(0)->GetXaxis()->SetTitle("Residual Range [cm]");
        }
        

        //-------- Draw SpacePoints
        pad4->cd();
        std::cout<<" SpacePoints Size: "<<SpacePointsX->size()<<std::endl;
        if(SpacePointsX->size()>0){
            // Function to get minimum and maximum XYZ values of the vector of vectors
            auto getMinMax = [](std::vector<double> v){
                double min = *std::min_element(v.begin(), v.end());
                double max = *std::max_element(v.begin(), v.end());
                return std::make_pair(min, max);
            };
        

            // Global min and max values
            std::pair<double, double> xMinMaxGlobal = {200, -200};
            std::pair<double, double> yMinMaxGlobal = {200, -200};
            std::pair<double, double> zMinMaxGlobal = {500, 0};

        
            std::vector<TGraph2D*> graphList;
            int graphCounter = 0;
            for(int j=0; j<SpacePointsX->size(); j++){
                std::cout<<"SpacePointsX size: "<<SpacePointsX->at(j).size()<<std::endl;
                if(SpacePointsX->at(j).size()==0) continue;
                
                // Get min and max values of X, Y and Z
                std::pair<double, double> xMinMax = getMinMax(SpacePointsX->at(j));
                std::pair<double, double> yMinMax = getMinMax(SpacePointsY->at(j));
                std::pair<double, double> zMinMax = getMinMax(SpacePointsZ->at(j));
                // Update global min and max values
                xMinMaxGlobal = std::make_pair(std::min(xMinMax.first, xMinMaxGlobal.first), std::max(xMinMax.second, xMinMaxGlobal.second));
                yMinMaxGlobal = std::make_pair(std::min(yMinMax.first, yMinMaxGlobal.first), std::max(yMinMax.second, yMinMaxGlobal.second));
                zMinMaxGlobal = std::make_pair(std::min(zMinMax.first, zMinMaxGlobal.first), std::max(zMinMax.second, zMinMaxGlobal.second));
            
                TGraph2D *graph = new TGraph2D();
                addPointsAndCreateGraph(graph, SpacePointsZ->at(j), SpacePointsX->at(j), SpacePointsY->at(j));    
                graphList.push_back(graph);
            
            }

            TH3F *hFrame = new TH3F("frame",";Z [cm]; X [cm]; Y [cm]", nBins,zMinMaxGlobal.first,zMinMaxGlobal.second, nBins, xMinMaxGlobal.first,xMinMaxGlobal.second, nBins,yMinMaxGlobal.first,yMinMaxGlobal.second);
            hFrame->Draw();
            for(int j=0; j<graphList.size(); j++){
                drawGraph(graphList.at(j), "SpacePoints", fColorList.at(j));
            }
            pad4->Update();     
        }

        c1->cd();
        c1->Update();
        c1->WaitPrimitive();
        
    }


    

    
    
    
    

}