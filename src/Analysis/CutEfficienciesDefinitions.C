#ifndef CUTEFFICIENCIES_DEFINITIONS_H
#define CUTEFFICIENCIES_DEFINITIONS_H

#include "CutEfficienciesStyle.C"

class SampleDef {
public:
    SampleDef(const std::string& var = "", const std::string& label = "", bool isSignal = false, const std::string& weight = "1")
        : fVar(var), fLabel(label), fIsSignal(isSignal), fWeight(weight), fNEvents(0) {
    }

    // Getter methods
    TString GetVar() const {
        return fVar;
    }

    TString GetLabel() const {
        return fLabel;
    }

    std::string GetVarS() const {
        return fVar.Data();
    }

    std::string GetLabelS() const {
        return fLabel.Data();
    }

    TString GetWeight() const {
        return fWeight;
    }

    std::string GetWeightS() const {
        return fWeight.Data();
    }

    bool IsSignal() const {
        return fIsSignal;
    }

    // Setter methods
    void SetVar(const TString& var) {
        fVar = var;
    }

    void SetLabel(const TString& label) {
        fLabel = label;
    }

    void SetIsSignal(bool isSignal) {
        fIsSignal = isSignal;
    }

    void SetWeight(const TString& weight) {
        fWeight = weight;
    }

    void SetNEvents(int nEvents) {
        fNEvents = nEvents;
    }

    int GetNEvents() const {
        return fNEvents;
    }

private:
    TString fVar;
    TString fLabel;
    bool fIsSignal;
    TString fWeight;
    int fNEvents;
};


class Binning {
public:
    Binning(double x1 = 0, double x2 = 1, int nBins = 2)
        : fX1(x1), fX2(x2), fNBins(nBins) {
    }

    // Getter methods
    double GetX1() const {
        return fX1;
    }

    double GetX2() const {
        return fX2;
    }

    int GetNBins() const {
        return fNBins;
    }

    // Setter methods
    void SetX1(double x1) {
        fX1 = x1;
    }

    void SetX2(double x2) {
        fX2 = x2;
    }

    void SetNBins(int nBins) {
        fNBins = nBins;
    }

private:
    double fX1;
    double fX2;
    int fNBins;
};

enum class CutType {
    kNone=-1,
    kRight=0,
    kLeft,
    kCenter,
    k2D
};

class PlotDef {
public:
    PlotDef(const TString& var = "",
            const TString& cut = "",
            const CutType cutType = CutType::kRight,
            double cutValue = 0,
            const Binning& bins = Binning(0, 1, 2),
            bool accumulateCut = true,
            const TString& varLabel = "",
            const TString& cutLabel = "",
            bool log = false,
            const TString& suffix = "")
        : fVar(var),
          fCut(cut),
          fCutType(cutType),
          fCutValue(cutValue),
          fBins(bins),
          fAccumulateCut(accumulateCut),
          fVarLabel(varLabel),
          fCutLabel(cutLabel),
          fLog(log),
          fSuffix(suffix)
    {
        if(cutType == CutType::kNone){
            fCut = "0==0";
        }
        else if(cutType == CutType::kCenter){
            fCut =  fCut + ("==" + std::to_string(cutValue)).c_str();
        }
        else if(cutType == CutType::kRight){
            fCut =  fCut + (">=" + std::to_string(cutValue)).c_str();
            std::ostringstream streamObj;
            streamObj <<  std::fixed << std::setprecision(2) << cutValue;
            fCutLabel = fCutLabel + ("\\ >= " + streamObj.str()).c_str();
        }
        else if(cutType == CutType::kLeft){
            fCut =  fCut + ("<=" + std::to_string(cutValue)).c_str();
            std::ostringstream streamObj;
            streamObj <<  std::fixed << std::setprecision(2) << cutValue;
            fCutLabel = fCutLabel + ("\\ <=" + streamObj.str()).c_str();
        }
    }

    // Getter methods
    TString GetSuffix() const {
        return fSuffix;
    }

    TString GetVar() const {
        return fVar;
    }

    std::string GetVarS() const {
        return fVar.Data();
    }
    
    TString GetCut() const {
        return fCut;
    }

    double GetCutValue() const {
        return fCutValue;
    }

    TString GetVarLabel(bool useStdString=0) const {
        return fVarLabel;
    }

    std::string GetVarLabelS(bool useStdString=0) const {
        return fVarLabel.Data();
    }

    TString GetCutLabel() const {
        return fCutLabel;
    }

    CutType GetCutType() const {
        return fCutType;
    }

    Binning GetBins() const {
        return fBins;
    }

    bool GetAccumulateCut() const {
        return fAccumulateCut;
    }

    bool GetLog() const {
        return fLog;
    }

    // Setter methods
    void SetSuffix(const TString& suffix) {
        fSuffix = suffix;
    }

    void SetVar(const TString& var) {
        fVar = var;
    }

    void SetCut(const TString& cut) {
        fCut = cut;
    }

    void SetVarLabel(const TString& varLabel) {
        fVarLabel = varLabel;
    }

    void SetCutLabel(const TString& cutLabel) {
        fCutLabel = cutLabel;
    }

    void SetBins(const Binning& bins) {
        fBins = bins;
    }

    void SetAccumulateCut(bool accumulateCut) {
        fAccumulateCut = accumulateCut;
    }

    void SetLog(bool log) {
        fLog = log;
    }

private:
    TString fVar;
    TString fCut;
    CutType fCutType;
    double fCutValue;
    Binning fBins;
    bool fAccumulateCut;
    TString fVarLabel;
    TString fCutLabel;
    bool fLog;
    TString fSuffix;
};



class AnaPlot{
    public:
        AnaPlot(int plotIndex, PlotDef plotDef, std::vector<SampleDef> intTypes);

        void DrawHistograms(TTree* fTree, TCut currentCut, bool afterCut=0);

        PlotDef GetPlotDef() const {
            return fPlotDef;
        };

        std::map<std::string, int> GetCountsV() const {
            return fCountsV;
        };



    private:
        
        int fPlotIndex;
        PlotDef fPlotDef;
        
        int fNTypes;
        std::vector<SampleDef> fIntTypes;

        std::map<std::string, TH1F*> fHistV;
        std::map<std::string, TH1F*> fHistCumulativeV;
        std::map<std::string, int> fCountsV;

        std::map<std::string, TH2F*> fHist2DV;
        std::map<std::string, TH2F*> fHist2DCumulativeV;

        TCanvas *fCanvas;
        CutStyler *fStyler;
};


AnaPlot::AnaPlot(int cutIndex, PlotDef plotDef, std::vector<SampleDef> intTypes)
    : fPlotIndex(cutIndex),
    fPlotDef(plotDef),
    fIntTypes(intTypes),
    fNTypes(intTypes.size()),
    fCanvas(new TCanvas( ("c_"+std::to_string(fPlotIndex)).c_str(), "Selection", 800,600)),
    fStyler(new CutStyler(0))
{

    //--------- Hisogram initializations
    for (size_t j = 0; j < fIntTypes.size(); ++j){
        
        std::string intTypeLabel = fIntTypes[j].GetLabelS();
        std::string plotLabel = std::to_string(fPlotIndex)+"_"+std::to_string(j);

        TH1F *hAux = new TH1F(plotLabel.c_str(), plotLabel.c_str(), fPlotDef.GetBins().GetNBins(), fPlotDef.GetBins().GetX1(), fPlotDef.GetBins().GetX2());
        TH1F *hAuxCumulative = new TH1F( (plotLabel+"C").c_str(), plotLabel.c_str(), fPlotDef.GetBins().GetNBins(), fPlotDef.GetBins().GetX1(), fPlotDef.GetBins().GetX2());

        fCountsV[intTypeLabel] = 0;
        fHistV[intTypeLabel] = hAux;
        fHistCumulativeV[intTypeLabel] = hAuxCumulative;

        // if 2D Cut
        if(fPlotDef.GetCutType()==CutType::k2D){
            TH2F *hAux2D = new TH2F(plotLabel.c_str(), plotLabel.c_str(), fPlotDef.GetBins().GetNBins(), fPlotDef.GetBins().GetX1(), fPlotDef.GetBins().GetX2(), fPlotDef.GetBins().GetNBins(), fPlotDef.GetBins().GetX1(), fPlotDef.GetBins().GetX2());
            TH2F *hAux2DCumulative = new TH2F((plotLabel+"C").c_str(), plotLabel.c_str(), fPlotDef.GetBins().GetNBins(), fPlotDef.GetBins().GetX1(), fPlotDef.GetBins().GetX2(), fPlotDef.GetBins().GetNBins(), fPlotDef.GetBins().GetX1(), fPlotDef.GetBins().GetX2());
            
            fHist2DV[intTypeLabel] = hAux2D;
            fHist2DCumulativeV[intTypeLabel] = hAux2DCumulative;
        }

        
    }

}


void DrawVerticalLineWithArrow(double x, double x1, double x2, double y1, double y2, CutType cutType) {

    if(cutType == CutType::kNone) return;
    // Get Y range value
    double minY = std::min(y1, y2);
    double maxY = std::max(y1, y2);

    // Get X range value
    double x_min = std::min(x1, x2);
    double x_max = std::max(x1, x2);    

    double minX = x;
    double maxX = x;
    if(cutType == CutType::kRight ){
        maxX = x_max;
    }
    else if(cutType == CutType::kLeft ){
        minX = x_min;
    }
    else{
        maxX = x+1;
    }
    
   
    if(cutType == CutType::kCenter ){
        // Create a TLine object for the vertical line
        TLine *line1 = new TLine(minX, minY, minX, maxY);
        line1->SetLineColor(40);
        line1->SetLineWidth(2);
        line1->Draw("same");
        TLine *line2 = new TLine(maxX, minY, maxX, maxY);
        line2->SetLineColor(40);
        line2->SetLineWidth(2);
        line2->Draw("same");

        // Calculate the position for the arrowhead (in the middle of the line)
        double yMid = (line1->GetY1() + line1->GetY2()) / 2;
        std::string marker = "|>";
        TArrow *arrow1 = new TArrow(minX - 0.02, yMid, minX + 0.02, yMid, 0.02, marker.c_str());
        arrow1->SetLineColor(40);
        arrow1->SetFillColor(40);
        arrow1->Draw();

        // Calculate the position for the arrowhead (in the middle of the line)
        yMid = (line2->GetY1() + line2->GetY2()) / 2;
        marker = "<|";
        TArrow *arrow2 = new TArrow(maxX - 0.02, yMid, maxX + 0.02, yMid, 0.02, marker.c_str());
        arrow2->SetLineColor(40);
        arrow2->SetFillColor(40);
        arrow2->Draw();
    }
    else{
        // Create a TLine object for the vertical line
        if(cutType==CutType::kLeft) x+=1;
        TLine *line = new TLine(x, minY, x, maxY);
        line->SetLineColor(40);
        line->SetLineWidth(2);
        line->Draw("same");

        // Calculate the position for the arrowhead (in the middle of the line)
        double yMid = (line->GetY1() + line->GetY2()) / 2;
        std::string marker = (cutType == CutType::kRight)? "|>":"<|";
        TArrow *arrow = new TArrow(x - 0.02, yMid, x + 0.02, yMid, 0.02, marker.c_str());
        arrow->SetLineColor(40);
        arrow->SetFillColor(40);
        arrow->Draw();
    }
    
    // Create a TBox object for the highlight box
    if(cutType==CutType::kLeft) maxX+=1;
    TBox *highlightBox = new TBox(minX, minY, maxX, maxY);
    highlightBox->SetFillColorAlpha(40, 0.2); // Set transparent green color
    highlightBox->SetFillStyle(1001); // Set the fill style
    highlightBox->Draw("same");

    return;
}

void AnaPlot::DrawHistograms(TTree* fTree, TCut currentCut, bool afterCut){

    // --- Frame histograms and legend
    std::string plotIndex = std::to_string(fPlotIndex);
    
    if(fPlotDef.GetCutType()==CutType::k2D){
        std::string plotLabel = "hAux2D"+plotIndex;
        std::cout<<"2D plots not implemented yet"<<std::endl;
        TCanvas *c3 = new TCanvas(("c3"+plotIndex).c_str(),"2DSelection",0, 0, 1000, 750);
        c3->cd();

        TPad *pad1 = new TPad("pad1","pad1",0,0.5,0.5,1);
        pad1->Draw();
        TPad *pad2 = new TPad("pad2","pad2",0.5,0.5,1,1);
        pad2->Draw();
        TPad *pad3 = new TPad("pad3","pad3",0,0,0.5,0.5);
        pad3->Draw();
        TPad *pad4 = new TPad("pad4","pad4",0.5,0,1,0.5);
        pad4->Draw();
        
         
        std::string plotAxisLabels="1;# hits track 1;# hits track 2";

        plotLabel = "hAux2D";
        TH2F *hAux2D = new TH2F(plotLabel.c_str(), plotAxisLabels.c_str(), fPlotDef.GetBins().GetNBins(), fPlotDef.GetBins().GetX1(), fPlotDef.GetBins().GetX2(), fPlotDef.GetBins().GetNBins(), fPlotDef.GetBins().GetX1(), fPlotDef.GetBins().GetX2());

        std::cout<<" Cut name " << fPlotDef.GetVarS() << std::endl; 
        for (size_t j = 0; j < fIntTypes.size(); ++j) {
        
           
            
            TCut sampelCut(fIntTypes[j].GetVar());
            std::string plotLabel = "hAux2D";

             // --- Interaction type
            std::cout<<j<<" "<<sampelCut<<std::endl;
            // --- Draw distribution
            pad1->cd();
            fTree->Draw( (fPlotDef.GetVarS()+">>"+std::to_string(fPlotIndex)+"_"+std::to_string(j)).c_str(), sampelCut && currentCut );
        }

        hAux2D->Draw();
        
        for (size_t j = 0; j < fIntTypes.size(); ++j) {
            std::string intTypeLabel = fIntTypes[j].GetLabelS();
            fHist2DV[intTypeLabel]->Draw("same");
            // set marker styles and colors
            fHist2DV[intTypeLabel]->SetLineColor(fStyler->GetColor(j));
            fHist2DV[intTypeLabel]->SetLineWidth(2);
            fHist2DV[intTypeLabel]->SetMarkerStyle(fStyler->GetMarkerStyle(j));
            fHist2DV[intTypeLabel]->SetMarkerSize(1);
            fHist2DV[intTypeLabel]->SetMarkerColorAlpha(fStyler->GetColor(j), 0.6);
        }

        c3->Update();
        c3->SaveAs( ("OutputPlots/test2D_"+std::to_string(fPlotIndex)+".pdf").c_str());
        return;
    }

    std::cout<<"Making plot "<<fPlotIndex<<std::endl;


    TCanvas *c2 = new TCanvas(("c2"+plotIndex).c_str(),"Selection",0, 0, 1000, 750);
    c2->cd();

    TPad *pad1 = new TPad("pad1","pad1",0,0.5,0.5,1);
    pad1->Draw();
    TPad *pad2 = new TPad("pad2","pad2",0.5,0.5,1,1);
    pad2->Draw();
    TPad *pad3 = new TPad("pad3","pad3",0,0,0.5,0.5);
    pad3->Draw();
    TPad *pad4 = new TPad("pad4","pad4",0.5,0,1,0.5);
    pad4->Draw();
    
    std::string plotAxisLabels=";"+fPlotDef.GetVarLabelS()+"; # events";

    std::string plotLabel = "hAux"+std::to_string(fPlotIndex);
    TH1F *hAux = new TH1F( plotLabel.c_str(), plotAxisLabels.c_str(), fPlotDef.GetBins().GetNBins(), fPlotDef.GetBins().GetX1(), fPlotDef.GetBins().GetX2());

    plotLabel = "hAuxEff"+std::to_string(fPlotIndex);
    TH1F *hAuxEff = new TH1F(plotLabel.c_str(), plotAxisLabels.c_str(), fPlotDef.GetBins().GetNBins(), fPlotDef.GetBins().GetX1(), fPlotDef.GetBins().GetX2());

    plotLabel = "hAuxSignificance"+std::to_string(fPlotIndex);
    TH1F *hAuxSignificance = new TH1F(plotLabel.c_str(), plotAxisLabels.c_str(), fPlotDef.GetBins().GetNBins(), fPlotDef.GetBins().GetX1(), fPlotDef.GetBins().GetX2());

    plotLabel = "hCounter"+std::to_string(fPlotIndex);
    TH1F *hCounter = new TH1F( plotLabel.c_str(), plotAxisLabels.c_str(), 1, 0, 1);

    TLegend *legend = new TLegend(0.75, 0.75, 0.9, 0.9);
       
    pad1->cd();
    hAux->Draw();
    int maxVal=0;

    for (size_t j = 0; j < fIntTypes.size(); ++j) {
        
        // --- Interaction type
        fHistV[fIntTypes[j].GetLabelS()]->Draw();
        TCut sampelCut(fIntTypes[j].GetVar());
       
        std::string intTypeLabel = fIntTypes[j].GetLabelS();
        std::string plotLabel = std::to_string(fPlotIndex)+"_"+std::to_string(j);

        // --- Draw distribution
        pad1->cd();
        fTree->Draw( (fPlotDef.GetVarS()+">>"+plotLabel).c_str(), sampelCut && currentCut );
        
        
        // --- Draw cumulative distribution
        pad3->cd();

        if(fPlotDef.GetCutType()==CutType::kCenter){
            fHistCumulativeV[intTypeLabel] = (TH1F*)fHistV[intTypeLabel]->Clone();
        }
        else if(fPlotDef.GetCutType()==CutType::kRight){
            if(fIntTypes[j].IsSignal()){
                fHistCumulativeV[intTypeLabel] = (TH1F*)fHistV[intTypeLabel]->GetCumulative(false);
            }
            else{
                fHistCumulativeV[intTypeLabel] = (TH1F*)fHistV[intTypeLabel]->GetCumulative();
            }
        }
        else{
            if(fIntTypes[j].IsSignal()){
                fHistCumulativeV[intTypeLabel] = (TH1F*)fHistV[intTypeLabel]->GetCumulative();
            }
            else{
                fHistCumulativeV[intTypeLabel] = (TH1F*)fHistV[intTypeLabel]->GetCumulative(false);
            }
        }



        pad1->cd();
        // --- Get counts
        plotLabel = "hCounter"+std::to_string(fPlotIndex);
        fTree->Draw( ("0>>"+plotLabel).c_str(), fIntTypes[j].GetWeight()*TCut(sampelCut && currentCut && fPlotDef.GetCut()), "hist");
        std::cout<<hCounter->Integral()<<std::endl;
        fCountsV[intTypeLabel] = hCounter->Integral();
        //fCountsV[intTypeLabel] = fTree->Draw( "", "4"*TCut(sampelCut && currentCut && fPlotDef.GetCut()), "hist").GetIntegral();
        //fCountsV[intTypeLabel] = fTree->Draw( "", "4*RecoIsFiducial", "goff");


        // --- Set styles
        fHistV[intTypeLabel]->SetLineColor(fStyler->GetColor(j)); 
        fHistV[intTypeLabel]->SetLineWidth(2);
        fHistV[intTypeLabel]->SetLineStyle(fStyler->GetLineStyle(j));
        fHistV[intTypeLabel]->SetMarkerStyle(fStyler->GetMarkerStyle(j));
        fHistV[intTypeLabel]->SetMarkerColor(fStyler->GetColor(j));
        fHistV[intTypeLabel]->SetMarkerSize(1);
        fHistV[intTypeLabel]->SetStats(0);
        legend->AddEntry(fHistV[intTypeLabel], fIntTypes[j].GetLabelS().c_str(), "lp");
        

        fHistCumulativeV[intTypeLabel]->SetLineColor(fStyler->GetColor(j));
        fHistCumulativeV[intTypeLabel]->SetLineWidth(2);
        fHistCumulativeV[intTypeLabel]->SetMarkerStyle(fStyler->GetMarkerStyle(j));
        fHistCumulativeV[intTypeLabel]->SetMarkerColor(fStyler->GetColor(j));
        fHistCumulativeV[intTypeLabel]->SetMarkerSize(1);
        fHistCumulativeV[intTypeLabel]->Scale(100./fHistV[intTypeLabel]->Integral());
        fHistCumulativeV[intTypeLabel]->SetLineStyle(fStyler->GetLineStyle(j));
        fHistCumulativeV[intTypeLabel]->SetStats(0);
        

        // --- Set maximum
        if(fHistV[intTypeLabel]->GetMaximum()>maxVal) maxVal = fHistV[intTypeLabel]->GetMaximum();

    }

    hAux->SetStats(0);
    // check hAux number of bins is the same as the range of the plot
    if( std::abs( fPlotDef.GetBins().GetX1()-fPlotDef.GetBins().GetX2() ) == fPlotDef.GetBins().GetNBins()){
        hAux->SetNdivisions(hAux->GetXaxis()->GetNbins());
    }
    legend->SetFillColorAlpha(0, 0);
    
    pad1->cd();
    hAux->Draw();
    hAux->GetYaxis()->SetRangeUser(0., maxVal*1.1);
    for (size_t j = 0; j < fIntTypes.size(); ++j) {
        fHistV[fIntTypes[j].GetLabelS()]->Draw("same hist");
    }
    
    


    DrawVerticalLineWithArrow( fPlotDef.GetCutValue(), hAux->GetXaxis()->GetXmin(),  hAux->GetXaxis()->GetXmax(), 0, maxVal*1.1, fPlotDef.GetCutType() );
    legend->Draw("same");
   
    pad3->cd();
    hAux->Draw();
    hAux->GetYaxis()->SetRangeUser(0.5, maxVal*1.1);
    pad3->SetLogy();
    for (size_t j = 0; j < fIntTypes.size(); ++j) {
        fHistV[fIntTypes[j].GetLabelS()]->Draw("same hist");
    }
    DrawVerticalLineWithArrow( fPlotDef.GetCutValue(), hAux->GetXaxis()->GetXmin(),  hAux->GetXaxis()->GetXmax(), 0, maxVal*1.1, fPlotDef.GetCutType() );
    legend->Draw("same"); 
    

    pad2->cd();
    hAuxEff->GetYaxis()->SetRangeUser(-0.1, 101);
    hAuxEff->GetYaxis()->SetTitle("Efficiency [%]");
    hAuxEff->SetStats(0);
    hAuxEff->Draw();

    gStyle->SetOptLogy(0);
    for (size_t j = 0; j < fIntTypes.size(); ++j) {
        fHistCumulativeV[fIntTypes[j].GetLabelS()]->Draw("same p lhist");
    }
    legend->Draw("same");
    pad2->Update();

    // pad4 with signal over background
    pad4->cd();
    hAuxSignificance->SetStats(0);
    hAuxSignificance->GetYaxis()->SetTitle("S/#sqrt{BG}");
    hAuxSignificance->Draw();
    // new Th1 to sotre the signal
    TH1F *hAuxSigSignal = new TH1F("hAuxSigSignal", plotAxisLabels.c_str(), fPlotDef.GetBins().GetNBins(), fPlotDef.GetBins().GetX1(), fPlotDef.GetBins().GetX2());
    for (size_t j = 0; j < fIntTypes.size(); ++j) {
        if(fIntTypes[j].IsSignal()){
            hAuxSigSignal->Add(fHistCumulativeV[fIntTypes[j].GetLabelS()]);
        }
    }

    
    // now draw all the other backgrounds
    // vector to store TGrapsh
    std::vector<TGraph*> grV;
    // minmium and maximum Y values
    double minY = 1000000;
    double maxY = 0;
    for (size_t j = 0; j < fIntTypes.size(); ++j) {
        if(!fIntTypes[j].IsSignal()){
           
            //create x and y vectors and loop over bins, add signal/sqrt(bg)
            std::vector<double> x;
            std::vector<double> y;
            for (int i = 1; i <= fPlotDef.GetBins().GetNBins(); i++) {
                double binValueBG = (100-fHistCumulativeV[fIntTypes[j].GetLabelS()]->GetBinContent(i));
                double binValueSignal = hAuxSigSignal->GetBinContent(i);
                double binCenter = fHistV[fIntTypes[j].GetLabelS()]->GetBinCenter(i);
                if(binValueBG>0){
                    x.push_back(binCenter);
                    y.push_back(binValueSignal/std::sqrt(binValueBG));
                }
            }

            // create TGraph
            TGraph *gr = new TGraph(x.size(), &x[0], &y[0]);
            gr->SetLineColor(fStyler->GetColor(j));
            gr->SetLineWidth(2);
            gr->SetMarkerStyle(fStyler->GetMarkerStyle(j));
            gr->SetMarkerColor(fStyler->GetColor(j));
            gr->SetMarkerSize(1);
            gr->SetLineStyle(fStyler->GetLineStyle(j));

            grV.push_back(gr);

            double min = gr->GetYaxis()->GetXmin();
            double max = gr->GetYaxis()->GetXmax();
            if(min<minY) minY = min;
            if(max>maxY) maxY = max;
        }
    }

    // now draw all the other backgrounds

    hAuxSignificance->GetYaxis()->SetRangeUser(minY*0.9, maxY*1.1);
    for (size_t j = 0; j < grV.size(); ++j) {
        grV[j]->Draw("same lp");
    }
    
    legend->Draw("same");
    pad4->Update();

    // --- Save canvas
    c2->Update();
    std::string stageLabel = afterCut ? "1after" : "0before";
    c2->SaveAs(("OutputPlots/plot"+std::to_string(fPlotIndex)+stageLabel+fPlotDef.GetVarS()+".pdf").c_str());
    TFile *fFile = new TFile("OutputPlots/OutputPlots.root","UPDATE");
    fFile->cd();
    c2->Write();
    fFile->Close();
    
    return;
    
}

#endif