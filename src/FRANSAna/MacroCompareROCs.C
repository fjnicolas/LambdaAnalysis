#include "TFile.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include <string>
#include <vector>

std::vector<int> colors = {kBlue, kRed, kGreen, kOrange, kMagenta, kCyan, kYellow, kBlack};

//std::vector<std::string> myFiles = {"TMVAResultsFinal/TMTMVAResultsNew2.root", "TMVAResultsFinal/TMTMVAResultsNew4.root", "TMVAResultsFinal/TMTMVAResultsOld.root", "TMVAResultsFinal/TMTMVAResultsNew4NoAlpha.root"}; std::vector<std::string> myLabels = {"New2", "New4", "Old", "New4 no alpha"};

//std::vector<std::string> myFiles = {"TMVAResultsFinal/TMTMVAResultsResNew4Alpha.root", "TMVAResultsFinal/TMTMVAResultsResNew4NoAlpha.root", "TMVAResultsFinal/TMTMVAResultsResOld.root"};
//std::vector<std::string> myLabels = {"Alpha", "NoAlpha", "Old"};

//std::vector<std::string> myFiles = {"TMVAResultsFinal/TMTMVAResultsFullNew2.root", "TMVAResultsFinal/TMTMVAResultsFullNew4.root", "TMVAResultsFinal/TMTMVAResultsFullOld.root"}; std::vector<std::string> myLabels = {"New2", "New4", "Old"};

std::vector<std::string> myFiles = {"TMVAResultsFinal/TMVAResultsStd/TMVAResults_reco.root"};
std::vector<std::string> myLabels = {"Old"};

void CompareROCCurves(const std::vector<std::string>& fileNames = myFiles, const std::vector<std::string>& labels = myLabels) {
    TCanvas* c = new TCanvas("c", "ROC Comparison", 800, 600);
    c->Divide(1,1);


    c->cd(1);
    gPad->SetBottomMargin(0.15);
    gPad->SetLeftMargin(0.15);


    std::vector<TH1F*> rocs;
    for (size_t i = 0; i < fileNames.size(); i += 1) {
        std::string path = "/Users/franciscojaviernicolas/Work/HyperonsAna/SampleLines/AnaNewFRANS/";
        TFile* file1 = TFile::Open((path + fileNames[i]).c_str());
       
        std::string fPlotName = "dataset/Method_BDT/BDT/MVA_BDT_rejBvsS";

        TH1F* roc = (TH1F*)file1->Get(fPlotName.c_str());
        rocs.push_back(roc);
        std::cout<<i<<std::endl;
    }

    TH2F hFrame("hFrame", "ROC Comparison;Signal efficiency; BG rejection", 500, 0, 1, 500, 0, 1);
    hFrame.SetStats(0);
    hFrame.Draw();
    hFrame.GetXaxis()->SetTitleOffset(1.2);
    hFrame.GetYaxis()->SetTitleOffset(1.2);

    TLegend* legend = new TLegend(0.2, 0.2, 0.5, 0.4);
    
    for (size_t i = 0; i < fileNames.size(); i += 1) {

        std::cout<<i<<std::endl;

        rocs.at(i)->SetLineColor(colors.at(i));
        rocs.at(i)->SetLineWidth(2);
        rocs.at(i)->Draw("L SAME");        

        legend->AddEntry(rocs.at(i), labels[i].c_str(), "l");
        legend->Draw();
    }

    c->cd();
    c->Update();
    c->WaitPrimitive();
    c->SaveAs("ROC_Comparison.pdf");
}