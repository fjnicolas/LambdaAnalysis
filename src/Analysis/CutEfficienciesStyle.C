
#ifndef CUTEFFICIENCIES_STYLE_H
#define CUTEFFICIENCIES_STYLE_H

#include "TString.h"

// CutStyler class
class CutStyler{
    public:
        CutStyler(bool batchMode):
            fBatchMode(batchMode),
            fColors({kRed+1, kBlue-4, kGreen+4, kOrange+5, kMagenta+2, kCyan+3}),
            fLineStyle({1, 2, 9, 3, 6, 10}),
            fMarkerStyle({20, 21, 22, 23, 24, 25, 26, 27, 28, 29})
        {
            SetGStyle();
        }
        void SetGStyle();
        size_t GetColor(size_t i){return fColors[i];}
        size_t GetLineStyle(size_t i){return fLineStyle[i];}
        size_t GetMarkerStyle(size_t i){return fMarkerStyle[i];}

    private:
        bool fBatchMode;
        std::vector<int> fColors;
        std::vector<int> fLineStyle;
        std::vector<int> fMarkerStyle;
};


void CutStyler::SetGStyle(){
    
  gStyle->SetTitleFont(62, "TXYZ");
  gStyle->SetTitleSize(0.05,"TXYZ");

  // LABELS SIZE AND FONT
  gStyle->SetLabelFont(62, "TXYZ");
  gStyle->SetLabelSize(0.05,"TXYZ");

  // AXIS OFFSETS AND SIZES
  gStyle->SetTitleXOffset(0.80);
  gStyle->SetTitleXSize(0.06);

  gStyle->SetTitleYOffset(1.2);
  gStyle->SetTitleYSize(0.06);

  gStyle->SetMarkerStyle(33);
  gStyle->SetMarkerSize(1.5);
  gStyle->SetLineColor(46);
  gStyle->SetLineWidth(1);
  
  // Draw horizontal and vertical grids
  //gStyle->SetPadGridX(kTRUE);                     
  //gStyle->SetPadGridY(kTRUE);

  // set left margin to 10%
  gStyle->SetPadLeftMargin(0.15);

  //gStyle->SetOptLogy(1);
}

#endif