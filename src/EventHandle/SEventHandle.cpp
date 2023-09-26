////////////////////////////////////////////////////////////////////////////
//
// \file TPCLineParameters.h
//
// \brief Definition of TPCLinesParameters
//
// \author fjnicolas@ugr.es
//
////////////////////////////////////////////////////////////////////////////
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

#include "SEventHandle.h"

std::vector<TPad*> buildpadcanvas(int nx, int ny){
    std::vector<TPad*> Tp;
    double x=0, y=1, dx=1./nx, dy=1./ny;
    TPad *pad = new TPad("","", 0, 0, 1, 1, -1, -1, -1);
    Tp.push_back(pad);
    for(int i=1; i<=nx; i++){
      y=1;
      for(int j=1; j<=ny; j++){
        TPad *pad = new TPad("","", x, y-dy, x+dx, y, -1, -1, -1);
        Tp.push_back(pad);
        y-=dy;
      }
      x+=dx;
    }
    for(int i=0; i<=nx*ny; i++){
      Tp.at(i)->Draw();
      Tp.at(i)->SetBottomMargin(0.15);
    }
    return Tp;
}


std::vector<TString> GetInputFileList(std::string file_name,  std::string ext){
  // Get the candidate Files
  std::vector<TString> fFilePaths;
  TSystemDirectory dir(".", ".");
  TList *files = dir.GetListOfFiles();
  TString targetFileName(file_name);
  TString targetExtension(ext);
  std::cout<<" Target file name "<<targetFileName<<std::endl;
  gSystem->Exec("ls");
  if (files){
      TSystemFile *file;
      TString fname;
      TIter next(files);
      while ((file=(TSystemFile*)next())) {
          fname = file->GetName();
          
          if (!file->IsDirectory() && fname.EndsWith(targetExtension)){
              if(fname.Contains(targetFileName)){
                  fFilePaths.push_back(fname);
              }
              std::cout << fname << " Target:" << targetFileName <<std::endl;
              std::cout<<" Contains: "<<fname.Contains(targetFileName)<<std::endl;
          }
      }
  }

  return fFilePaths;
}