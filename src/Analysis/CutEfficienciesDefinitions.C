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
    kRightInt,
    kLeftInt,
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
            bool norm = false,
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
          fNormalize(norm),
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
        else if(cutType == CutType::kLeft || cutType == CutType::kLeftInt){
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

    bool GetNormalize() const {
        return fNormalize;
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
    bool fNormalize;
    bool fLog;
    TString fSuffix;
};



class AnaPlot{
    public:
        AnaPlot(int plotIndex, PlotDef plotDef, std::vector<SampleDef> intTypes, std::vector<PlotDef> phaseSpaceVaribles={});

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
        
        std::vector<PlotDef> fPhaseSpaceVars;
        std::vector<TH1F*> fHistVPhaseSpace0;
        std::vector<TH1F*> fHistVPhaseSpace1;

        TCanvas *fCanvas;
        CutStyler *fStyler;

        bool fDrawPhaseSpace;
};


AnaPlot::AnaPlot(int cutIndex, PlotDef plotDef, std::vector<SampleDef> intTypes,  std::vector<PlotDef> phaseSpaceVaribles)
    : fPlotIndex(cutIndex),
    fPlotDef(plotDef),
    fIntTypes(intTypes),
    fNTypes(intTypes.size()),
    fPhaseSpaceVars(phaseSpaceVaribles),
    fCanvas(new TCanvas( ("c_"+std::to_string(fPlotIndex)).c_str(), "Selection", 800,600)),
    fStyler(new CutStyler(0))
{

    //--------- Hisogram initializations
    for (size_t j = 0; j < fIntTypes.size(); ++j){
        
        fDrawPhaseSpace = (fPhaseSpaceVars.size()>0);

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

    //--------- Phase space initializations
    for (size_t j = 0; j < fPhaseSpaceVars.size(); ++j){
        std::string plotLabel0 = std::to_string(fPlotIndex)+"_"+std::to_string(j)+"_PS0";
        TH1F *hAux0 = new TH1F(plotLabel0.c_str(), (";"+fPhaseSpaceVars[j].GetVarLabelS()+";# events").c_str(), fPhaseSpaceVars[j].GetBins().GetNBins(), fPhaseSpaceVars[j].GetBins().GetX1(), fPhaseSpaceVars[j].GetBins().GetX2());
        
        std::string plotLabel1 = std::to_string(fPlotIndex)+"_"+std::to_string(j)+"_PS1";
        TH1F *hAux1 = new TH1F(plotLabel1.c_str(), (";"+fPhaseSpaceVars[j].GetVarLabelS()+";# events").c_str(), fPhaseSpaceVars[j].GetBins().GetNBins(), fPhaseSpaceVars[j].GetBins().GetX1(), fPhaseSpaceVars[j].GetBins().GetX2());
        
        hAux0->SetStats(0);
        hAux1->SetStats(0);
        hAux0->SetLineColor(kOrange-3);
        hAux1->SetLineColor(kAzure+2);
        hAux0->SetLineWidth(2);
        hAux1->SetLineWidth(2);
        hAux0->SetLineStyle(kSolid);
        hAux1->SetLineStyle(kDashed);

        fHistVPhaseSpace0.push_back(hAux0);
        fHistVPhaseSpace1.push_back(hAux1);    
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
    else if(cutType == CutType::kLeft || cutType == CutType::kLeftInt ){
        minX = x_min;
    }
    
    if(cutType == CutType::kLeftInt || cutType == CutType::kCenter ){
        maxX = maxX+1;
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
        if(cutType==CutType::kLeftInt) x+=1;
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
    std::cout<<"minX "<<minX<<" maxX "<<maxX<<" minY "<<minY<<" maxY "<<maxY<<std::endl;
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


    TCanvas *c2 = new TCanvas(("c2"+plotIndex).c_str(),"Selection",0, 0, 1200, 750);
    c2->cd();

    TPad *pad1 = new TPad("pad1","pad1", 0, 0.5, 0.33, 1);
    pad1->Draw();
    TPad *pad2 = new TPad("pad2","pad2", 0, 0, 0.33,0.5);
    pad2->Draw();
    TPad *pad3 = new TPad("pad3","pad3", 0.33, 0.5, 0.66, 1);
    pad3->Draw();
    TPad *pad4 = new TPad("pad4","pad4", 0.33, 0, 0.66, 0.5);
    pad4->Draw();
    TPad *pad5 = new TPad("pad5","pad5", 0.66, 0.5, 1, 1);
    pad5->Draw();
    TPad *pad6 = new TPad("pad6","pad6", 0.66, 0, 1, 0.5);
    pad6->Draw();

    std::string plotAxisLabels=";"+fPlotDef.GetVarLabelS()+"; # events";

    std::string plotLabel = "hAux"+std::to_string(fPlotIndex);
    TH1F *hAux = new TH1F( plotLabel.c_str(), plotAxisLabels.c_str(), fPlotDef.GetBins().GetNBins(), fPlotDef.GetBins().GetX1(), fPlotDef.GetBins().GetX2());

    plotLabel = "hAuxEff"+std::to_string(fPlotIndex);
    TH1F *hAuxEff = new TH1F(plotLabel.c_str(), plotAxisLabels.c_str(), fPlotDef.GetBins().GetNBins(), fPlotDef.GetBins().GetX1(), fPlotDef.GetBins().GetX2());

    plotLabel = "hAuxSignificance"+std::to_string(fPlotIndex);
    TH1F *hAuxSignificance = new TH1F(plotLabel.c_str(), plotAxisLabels.c_str(), fPlotDef.GetBins().GetNBins(), fPlotDef.GetBins().GetX1(), fPlotDef.GetBins().GetX2());

    plotLabel = "hCounter"+std::to_string(fPlotIndex);
    TH1F *hCounter = new TH1F( plotLabel.c_str(), plotAxisLabels.c_str(), 1, 0, 1);

    TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9);
       
    pad1->cd();
    hAux->Draw();
    int maxVal=0;
    int maxValInt=0;

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
        if(fHistV[intTypeLabel]->GetMaximum()>maxVal){
            maxVal = fHistV[intTypeLabel]->GetMaximum();
            maxValInt = fHistV[intTypeLabel]->Integral();
        }

        if(fDrawPhaseSpace){
            // --- Draw phase space
            if(fIntTypes[j].IsSignal()){


                for(int i=0; i<fPhaseSpaceVars.size(); i++){
                    std::string plotLabel0 = std::to_string(fPlotIndex)+"_"+std::to_string(i)+"_PS0";
                    fTree->Draw( (fPhaseSpaceVars[i].GetVarS()+">>"+plotLabel0).c_str(), sampelCut, "hist");
                    std::string plotLabel1 = std::to_string(fPlotIndex)+"_"+std::to_string(i)+"_PS1";
                    fTree->Draw( (fPhaseSpaceVars[i].GetVarS()+">>"+plotLabel1).c_str(), sampelCut && currentCut && fPlotDef.GetCut(), "hist");
                }
            }
        }

    }

    // ------ Draw histograms ------ 
    pad1->cd();

    hAux->SetStats(0);
    // check hAux number of bins is the same as the range of the plot
    if( std::abs( fPlotDef.GetBins().GetX1()-fPlotDef.GetBins().GetX2() ) == fPlotDef.GetBins().GetNBins()){
        hAux->SetNdivisions(hAux->GetXaxis()->GetNbins());
    }
    legend->SetFillColorAlpha(0, 0);
    
    bool fNormalize = fPlotDef.GetNormalize();
    if(fNormalize)
        hAux->GetYaxis()->SetRangeUser(0.,  1.1*(1.*maxVal/maxValInt));
    else
        hAux->GetYaxis()->SetRangeUser(0., maxVal*1.1);
    hAux->Draw();

    for (size_t j = 0; j < fIntTypes.size(); ++j) {
        if(fNormalize){
            fHistV[fIntTypes[j].GetLabelS()]->Scale(1./fHistV[fIntTypes[j].GetLabelS()]->Integral());
        }
        std::cout<<"Integral "<<fHistV[fIntTypes[j].GetLabelS()]->Integral()<<std::endl;
        fHistV[fIntTypes[j].GetLabelS()]->Draw("same hist");
    }
    

    DrawVerticalLineWithArrow( fPlotDef.GetCutValue(), hAux->GetXaxis()->GetXmin(),  hAux->GetXaxis()->GetXmax(), 0, maxVal*1.1, fPlotDef.GetCutType() );
    legend->Draw("same");
   
    pad2->cd();
    hAux->Draw();
    if(fNormalize)
        hAux->GetYaxis()->SetRangeUser(0.005, 1.1*(1.*maxVal/maxValInt));
    else
        hAux->GetYaxis()->SetRangeUser(0.5, maxVal*1.1);
    pad2->SetLogy();
    for (size_t j = 0; j < fIntTypes.size(); ++j) {
        fHistV[fIntTypes[j].GetLabelS()]->Draw("same hist");
    }
    DrawVerticalLineWithArrow( fPlotDef.GetCutValue(), hAux->GetXaxis()->GetXmin(),  hAux->GetXaxis()->GetXmax(), 0, maxVal*1.1, fPlotDef.GetCutType() );
    legend->Draw("same"); 
    

    // ------ Draw efficiency curves ------ 
    pad3->cd();
    hAuxEff->GetYaxis()->SetRangeUser(-0.1, 101);
    hAuxEff->GetYaxis()->SetTitle("Efficiency [%]");
    hAuxEff->SetStats(0);
    hAuxEff->Draw();


    gStyle->SetOptLogy(0);
    for (size_t j = 0; j < fIntTypes.size(); ++j) {
        fHistCumulativeV[fIntTypes[j].GetLabelS()]->Draw("same p lhist");
    }

    legend->Draw("same");
    pad3->Update();


    // ------ Draw Significance ------ 
    pad4->cd();
    hAuxSignificance->SetStats(0);
    

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

    // Legend for significance
    TLegend *legendSig = new TLegend(0.7, 0.7, 0.9, 0.9);

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

            // add to legend
            legendSig->AddEntry(gr, fIntTypes[j].GetLabelS().c_str(), "lp");

            double min = gr->GetYaxis()->GetXmin();
            double max = gr->GetYaxis()->GetXmax();
            if(min<minY) minY = min;
            if(max>maxY) maxY = max;
        }
    }

    hAuxSignificance->GetYaxis()->SetRangeUser(minY*0.9, maxY*1.1);
    for (size_t j = 0; j < grV.size(); ++j) {
        grV[j]->Draw("same lp");
    }
    
    if(grV.size()>1){
        legendSig->Draw("same");
    }
    pad4->Update();

    // ------ Draw ROC curve ------ 
    pad5->cd();

    // TH1 to sotre the signal
    TH1F *hAuxEffSignal = new TH1F("hAuxEffSignal", plotAxisLabels.c_str(), fPlotDef.GetBins().GetNBins(), fPlotDef.GetBins().GetX1(), fPlotDef.GetBins().GetX2());
    for (size_t j = 0; j < fIntTypes.size(); ++j) {
        if(fIntTypes[j].IsSignal()){
            hAuxEffSignal->Add(fHistCumulativeV[fIntTypes[j].GetLabelS()]);
        }
    }

    TH1F *hAuxROC = new TH1F("hAuxROC", "ROC;;", 200, 0, 100);
    hAuxROC->SetStats(0);
    hAuxROC->GetXaxis()->SetTitle("Signal efficiency [%]");
    hAuxROC->GetYaxis()->SetTitle("BG rejection [%]");
    hAuxROC->GetYaxis()->SetRangeUser(-0.1, 101);
    hAuxROC->Draw();

    std::vector<TGraph*> grVROC;
    for (size_t j = 0; j < fIntTypes.size(); ++j) {

        if(!fIntTypes[j].IsSignal()){
            //create x and y vectors and loop over bins, add signal/sqrt(bg)
            std::vector<double> x;
            std::vector<double> y;
            for (int i = 1; i <= fPlotDef.GetBins().GetNBins(); i++) {
                double binValueBG = (fHistCumulativeV[fIntTypes[j].GetLabelS()]->GetBinContent(i));
                double binValueSignal = hAuxEffSignal->GetBinContent(i);
                
                if(binValueBG>0){
                    x.push_back(binValueSignal);
                    y.push_back(binValueBG);
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

            grVROC.push_back(gr);   
        }
    }
    hAuxROC->Draw();
    legend->Draw("same");

    for (size_t j = 0; j < grVROC.size(); ++j) {
        grVROC[j]->Draw("same l");
    }
    pad5->Update();


    // ------ Draw phase space ------ 
    // Canvas
    TCanvas *cPhaseSpace = new TCanvas(("cPhaseSpace"+plotIndex).c_str(),"Phase Space",0, 0, 1200, 850);
    cPhaseSpace->cd();

    // vector of Pads
    std::vector<TPad*> padV = buildpadcanvas( (int) fPhaseSpaceVars.size()/2, (int) fPhaseSpaceVars.size()/2);

    TLegend *legendRP = new TLegend(0.7,0.7,0.9,0.9);

    // Draw each plot
    for(int i=0; i<fPhaseSpaceVars.size(); i++){
        padV[i+1]->cd();
        //fHistVPhaseSpace1[i]->GetXaxis()->SetTitle("x");
        fHistVPhaseSpace1[i]->GetYaxis()->SetTitle("y");
        fHistVPhaseSpace0[i]->GetYaxis()->SetTitle("y");
        TRatioPlot * rp = new TRatioPlot(fHistVPhaseSpace1[i], fHistVPhaseSpace0[i]);
        rp->SetH1DrawOpt("HIST");
        rp->SetH2DrawOpt("HIST SAME");
        rp->SetGraphDrawOpt("lp");
        
        rp->Draw();
        rp->SetLowBottomMargin(0.3);
        rp->SetLeftMargin(0.15);
        rp->SetSplitFraction(0.35);
        rp->SetSeparationMargin(0.);
        rp->GetLowerRefYaxis()->SetTitle("#epsilon");
        rp->GetUpperRefYaxis()->SetTitle("# entries");
        rp->GetLowerRefYaxis()->SetTitleOffset(1.);
        rp->GetUpperRefYaxis()->SetTitleOffset(1.);
        rp->GetLowYaxis()->SetNdivisions(1005);
        rp->GetUpperRefYaxis()->SetRangeUser(0., fHistVPhaseSpace0[i]->GetMaximum()*1.1);
        rp->GetLowerRefGraph()->SetMarkerStyle(20);
        rp->GetLowerRefGraph()->SetMarkerSize(1);
        rp->GetLowerRefGraph()->SetLineColor(kBlack);
        if(i==0){
            legendRP->AddEntry(fHistVPhaseSpace0[i],"All","l");
            legendRP->AddEntry(fHistVPhaseSpace1[i],"After cut","l");
        }
       
        legendRP->Draw("same");
    }
    cPhaseSpace->WaitPrimitive();

    // --- Save canvas
    c2->Update();
    cPhaseSpace->Update();
    //c2->WaitPrimitive();
    std::string stageLabel = afterCut ? "1after" : "0before";
    c2->SaveAs(("OutputPlots/plot"+std::to_string(fPlotIndex)+stageLabel+fPlotDef.GetVarS()+".pdf").c_str());
    if(fDrawPhaseSpace) cPhaseSpace->SaveAs(("OutputPlots/PhaseSpace/zphaseSpacePlot"+std::to_string(fPlotIndex)+stageLabel+fPlotDef.GetVarS()+"PhaseSpace.pdf").c_str());
    TFile *fFile = new TFile("OutputPlots/OutputPlots.root","UPDATE");
    fFile->cd();
    c2->Write();
    fFile->Close();
    
    return;
    
}

#endif