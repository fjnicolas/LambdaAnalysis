#include <cstdlib>
#include <iostream>
#include <string>
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TEfficiency.h"
#include "TMVA/Reader.h"

int MacroFRAMSEvaluate(
    const std::string& weightsPath = "FRAMSSelectionTMVA_BDT.weights.xml",
    const std::string& fileName = "LambdaAnaOutput.root",
    const std::string& treeName = "FRAMSTree"
) {
    const double scoreCut = 0.15;

    bool evaluteSignalOnly = true;

    TFile inputFile(fileName.c_str(), "READ");
    TTree* tree = static_cast<TTree*>(inputFile.Get(treeName.c_str()));

    TMVA::Reader reader;
    Float_t alpha_C, eta_C, delta_C, fitScore_C, gap, protonKE, pionKE;
    Int_t isSignal;

    //reader.AddVariable("Alpha_C", &alpha_C);
    reader.AddVariable("Eta_C", &eta_C);
    reader.AddVariable("Delta_C", &delta_C);
    reader.AddVariable("FitScore_C", &fitScore_C);
    reader.AddSpectator("Gap", &gap);
    reader.AddSpectator("ProtonKE", &protonKE);
    reader.AddSpectator("PionKE", &pionKE);
    reader.BookMVA("FRAMS BDT", weightsPath.c_str());

    double Alpha_C, Eta_C, Delta_C, FitScore_C, Gap, ProtonKE, PionKE;
    tree->SetBranchAddress("Alpha_C", &Alpha_C);
    tree->SetBranchAddress("Eta_C", &Eta_C);
    tree->SetBranchAddress("Delta_C", &Delta_C);
    tree->SetBranchAddress("FitScore_C", &FitScore_C);
    tree->SetBranchAddress("Gap", &Gap);
    tree->SetBranchAddress("ProtonKE", &ProtonKE);
    tree->SetBranchAddress("PionKE", &PionKE);
    tree->SetBranchAddress("IsSignal", &isSignal);

    TEfficiency* efficiency = new TEfficiency("effS", "Selection Efficiency;Gap [cm];#epsilon", 21, -1, 20);
    TH2F hC("hC", ";Gap [cm]; Score", 21, -1, 20, 40, -1., 1.);

    for (Long64_t ievt = 0; ievt < tree->GetEntries(); ievt++) {
        tree->GetEntry(ievt);

        alpha_C = Alpha_C;
        eta_C = Eta_C;
        delta_C = Delta_C;
        fitScore_C = FitScore_C;
        gap = Gap;
        protonKE = ProtonKE;
        pionKE = PionKE;

        if(evaluteSignalOnly && !isSignal) continue;

        // Retrieve the corresponding MVA output
        double score = reader.EvaluateMVA("FRAMS BDT");
        std::cout << ievt << " " << isSignal << " " << eta_C << " Gap=" << gap << " " << score << std::endl;

        bool passCut = score > scoreCut;
        efficiency->Fill(passCut, gap);
        hC.Fill(gap, score);
    }



    TCanvas* c1 = new TCanvas("example", "", 800, 400);
    c1->Divide(2, 1);
    

    c1->cd(1);
    gPad->SetBottomMargin(0.15);
    gPad->SetLeftMargin(0.15);
    efficiency->Draw("AP");
    gPad->Update();
    efficiency->GetPaintedGraph()->GetXaxis()->SetTitleOffset(1.2);
    efficiency->GetPaintedGraph()->GetYaxis()->SetTitleOffset(1.2);

    c1->cd(2);
    gPad->SetBottomMargin(0.15);
    gPad->SetLeftMargin(0.15);
    hC.GetXaxis()->SetTitleOffset(1.);
    hC.GetYaxis()->SetTitleOffset(1.);
    hC.Draw();

    c1->cd();
    c1->Update();
    c1->WaitPrimitive();

    delete efficiency;
    return 0;
}
