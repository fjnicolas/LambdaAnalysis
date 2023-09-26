////////////////////////////////////////////////////////////////////////////
//
// \file TPCLineParameters.h
//
// \brief Definition of TPCLinesParameters
//
// \author fjnicolas@ugr.es
//
////////////////////////////////////////////////////////////////////////////
#ifndef TPC_EVENTHANDLE_PARAMETERS_H
#define TPC_EVENTHANDLE_PARAMETERS_H

#include <vector>

// ROOT includes
#include "TH2F.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TFile.h"
#include "TString.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TSystemDirectory.h"

#include "TPCSimpleEvents.h"

std::vector<TPad*> buildpadcanvas(int nx, int ny);
std::vector<TString> GetInputFileList(std::string file_name,  std::string ext);

class SEventId {
    public:
        SEventId(int r=0, int sr=0, int e=0):
            fRun(r),
            fSubRun(sr),
            fEvent(e),
            fEventLabel(std::to_string(r)+"_"+std::to_string(sr)+"_"+std::to_string(e))
            {};


        SEventId(const SEventId& other) : SEventId(other.fRun, other.fSubRun, other.fEvent) {}

        std::string Label(){return fEventLabel;};
        int Run(){return fRun;};
        int SubRun(){return fSubRun;};
        int Event(){return fEvent;};

        friend std::ostream& operator<<(std::ostream& os, const SEventId& statusPrinter) {
        os << " ***** R=" << statusPrinter.fRun << " SR=" << statusPrinter.fSubRun << " E=" << statusPrinter.fEvent << " *****" << std::endl;
        return os;
    }

    private:
        int fRun, fSubRun, fEvent;
        std::string fEventLabel;
};

class SEventSelection {
    public:
        SEventSelection(SEventId ev):
            fEvent(ev),
            fNSlices(0),
            fNSelected(0),
            fNNotSelected(0),
            fNSkipped(0)
            {};

        SEventSelection(){
            fEvent = SEventId();
            fNSlices=0;
            fNSelected=0;
            fNNotSelected=0;
            fNSkipped=0;
        }
        
        SEventSelection(const SEventSelection& other) : SEventSelection(other.fEvent) {}

        SEventId EventId(){return fEvent;};

        void AddSelected(){
            fNSlices++;
            fNSelected++;
        };

        void AddSkipped(){
            fNSlices++;
            fNSkipped++;
        }

        void AddNotSelected(){
            fNSlices++;
            fNNotSelected++;
        }

        int NSlices(){return fNSlices;};
        int NSelected(){return fNSelected;};
        int NNotSelected(){return fNNotSelected;};
        int NSkipped(){return fNSkipped;};

    private:
        SEventId fEvent;
        int fNSlices;
        int fNSelected;
        int fNNotSelected;
        int fNSkipped;
};


class EfficiencyCalculator {
public:
    EfficiencyCalculator():
        fEventList({}),
        nEvents(0),
        nEventsSkipped(0),
        nProcessedEvents(0),
        nEventsSelected(0)
        {

            hNOrigins = new TH1F("hNOrigins", "NOrigins;# origins; # entries", 10, 0, 10);
            hNOriginsMult1 = new TH1F("hNOriginsMult1", "NOrigins mult=1;# origins; # entries", 10, 0, 10);
            hNOriginsMult2 = new TH1F("hNOriginsMult2", "NOrigins mult=2;# origins; # entries", 10, 0, 10);
            hNOriginsMultGt3 = new TH1F("hNOriginsMultGt3", "NOrigins mult>=3;# origins; # entries", 10, 0, 10);
            hProf2DMatrix = new TH2F("hProf2DMatrix", "# origins; origin 1 multiplicity; origin 2 multiplicity", 10, 0, 10, 10, 0, 10);

            hHitDensity = new TH1F("hHitDensity", ";HitDensity [#hits/wire]; # entries", 30, 0, 10);
        };

    void AddEvent(SEventId ev){
        if(fEventList.find(ev.Label()) == fEventList.end()){
            SEventSelection sEv(ev);
            fEventList[ev.Label()] = SEventSelection(sEv);
        }
    }

    void UpdateSkipped(SEventId ev){
        AddEvent(ev);
        fEventList[ev.Label()].AddSkipped();
    }

    void UpdateSelected(SEventId ev){
        AddEvent(ev);
        fEventList[ev.Label()].AddSelected();
    }

    void UpdateNotSelected(SEventId ev){
        AddEvent(ev);
        fEventList[ev.Label()].AddNotSelected();

    }



    void UpdateHistograms(SEvent recoEv){
        hNOrigins->Fill(recoEv.GetNOrigins());
        hNOriginsMult1->Fill(recoEv.GetNOriginsMult(1));
        hNOriginsMult2->Fill(recoEv.GetNOriginsMult(2));
        hNOriginsMultGt3->Fill(recoEv.GetNOriginsMultGt(3));
        
        /*for(int m1=0; m1<=8; m1++){
            for(int m2=0; m2<=8; m2++){
                if(recoEv.GetNOriginsMult(m1)>0 && recoEv.GetNOriginsMult(m2)>0){
                    hProf2DMatrix->Fill(m1, m2, 1);
                }
                
            }
        }*/

        // use only the two main origins
        std::vector<SOrigin> origins = recoEv.GetOrigins();

        if(origins.size()==1){
            hProf2DMatrix->Fill(origins[0].Multiplicity(), 0 );
        }
        else{
            int maxMult = std::max(origins[0].Multiplicity(), origins[1].Multiplicity());
            int minMult = std::min(origins[0].Multiplicity(), origins[1].Multiplicity());
            hProf2DMatrix->Fill(maxMult, minMult);
        }
       
        hHitDensity->Fill(recoEv.HitDensity());
    }

    void UpdateValues(){
        nEvents = fEventList.size();
        nEventsSkipped=0;
        nProcessedEvents=0;
        nEventsSelected=0;
        for(auto &pair:fEventList){
            if(pair.second.NSkipped()==pair.second.NSlices()){
                nEventsSkipped++;
            }
            else{
                nProcessedEvents++;
                if(pair.second.NSelected()>0){
                    nEventsSelected++;
                }
            }
        }
    }

    friend std::ostream& operator<<(std::ostream& os, EfficiencyCalculator& statusPrinter) {
        statusPrinter.UpdateValues();
        os << " ********** Final status-" << std::endl;
        os << " ... NTotalEvents=" << statusPrinter.nEvents << " NSkipped=" << statusPrinter.nEventsSkipped << std::endl;
        os << " ... NProcessed=" << statusPrinter.nProcessedEvents << " NSelected=" << statusPrinter.nEventsSelected;
        os << " efficiency=" << static_cast<double>(statusPrinter.nEventsSelected) / statusPrinter.nProcessedEvents;
        os << " efficiency all=" << static_cast<double>(statusPrinter.nEventsSelected) / statusPrinter.nEvents << std::endl;
        return os;
    }

    void DrawHistograms(){
        gROOT->SetBatch(false);
        gStyle->SetOptStat(1);
        TCanvas c("cStats", "Final stats", 0, 0, 1400,900);
        std::vector<TPad*> padV = buildpadcanvas(3, 2);

        padV[1]->cd();
        hNOrigins->Draw();
        
        padV[2]->cd();
        hNOriginsMult1->Draw();

        padV[3]->cd();
        hNOriginsMult2->Draw();

        padV[4]->cd();
        hNOriginsMultGt3->Draw();

        padV[5]->cd();
        hHitDensity->Draw();

        padV[6]->cd();
        hProf2DMatrix->SetStats(0);
        padV[6]->SetRightMargin(0.12);
        // Normalize the histogram to the number of entries
        if (hProf2DMatrix->GetEntries() > 0) {
            hProf2DMatrix->Scale(1.0 / hProf2DMatrix->GetEntries());
        }

        hProf2DMatrix->Draw("colz");

        int maxMult=5;
        // Add percentage values to each bin as labels
        for (Int_t i = 2; i <= maxMult; ++i) {
            for (Int_t j = 1; j <maxMult; ++j) {
                Double_t binContent = hProf2DMatrix->GetBinContent(i, j);
                // two to replace the symmetric matrix, double counting
                Double_t percentage =  binContent * 100.0;

                // Create a TText object for the label
                TText* label = new TText(hProf2DMatrix->GetXaxis()->GetBinCenter(i), hProf2DMatrix->GetYaxis()->GetBinCenter(j), Form("%.2f%%", percentage));
                label->SetTextSize(0.05); // Adjust the text size as needed
                label->SetTextAlign(22);  // Centered alignment
                label->SetTextColor(kRed);
                label->Draw();
            }
        }

        hProf2DMatrix->GetXaxis()->SetRangeUser(1, 5);
        hProf2DMatrix->GetYaxis()->SetRangeUser(0, 4);
        hProf2DMatrix->GetXaxis()->SetNdivisions(4);
        hProf2DMatrix->GetYaxis()->SetNdivisions(4);

        

        c.Update();

        c.cd();
        c.WaitPrimitive();

    }


private:
    std::map<std::string, SEventSelection> fEventList;
    int nEvents;
    int nEventsSkipped;
    int nProcessedEvents;
    int nEventsSelected;

    TH1F *hNOrigins;
    TH1F *hNOriginsMult1;
    TH1F *hNOriginsMult2;
    TH1F *hNOriginsMultGt3;
    TH2F *hProf2DMatrix;
    TH1F *hHitDensity;
};


#endif // TPC_LINES_PARAMETERS_H