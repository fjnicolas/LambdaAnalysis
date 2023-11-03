struct Binning {
    double x1;
    double x2;
    int nBins;

    Binning(const double _x1=0,
            const double _x2=1,
            const int _nBins=2)
        : x1(_x1),
          x2(_x2),
          nBins(_nBins) {
          }
};

struct PlotDef {

    std::string suffix;
    std::string var;
    std::string cut;
    std::string varLabel;
    std::string cutLabel;
    Binning bins;
    bool accumulateCut;
    bool log;

    PlotDef(const std::string& _suffix = "",
            const std::string& _var = "",
            const std::string& _cut = "",
            const std::string& _varLabel = "",
            const std::string& _cutLabel = "",
            const Binning& _bins = {0,1,2},
            bool _accumulateCut = true,
            bool _log = true)
        : suffix(_suffix),
          var(_var),
          cut(_cut),
          varLabel(_varLabel),
          cutLabel(_cutLabel),
          bins(_bins),
          accumulateCut(_accumulateCut),
          log(_log) {
          }
};


void generateAndCompileTeXTable(
    const std::vector<PlotDef>& cutDefs,
    const std::vector<std::string>& sampleDefNames,
    const std::vector<std::vector<int>>& histogramCounts,
    const std::string& fileName,
    const std::string& tableCaption
) {
    // Open a file stream for writing the .tex file
    std::ofstream texFile(fileName);

    // Check if the file is open
    if (!texFile.is_open()) {
        std::cerr << "Failed to open the output file." << std::endl;
        return;
    }

    // Write the LaTeX document preamble
    
    texFile << "\\documentclass{article}" << std::endl;
    texFile << "\\usepackage{graphicx}" << std::endl;
    texFile << "\\usepackage{pdflscape}"<<std::endl;
    texFile << "\\begin{document}" << std::endl;
    texFile << "\\begin{landscape}" << std::endl;


    // Write the LaTeX code for the table
    texFile << "\\begin{table}[h]" << std::endl;
    texFile << "\\centering" << std::endl;
    texFile << "\\begin{tabular}{|c|";
    for(size_t k = 0; k<sampleDefNames.size(); k++) texFile<<"c|";
    texFile << "}" << std::endl;
    texFile << "\\hline" << std::endl;

    // Write the header row
    texFile << " Cut & ";
    for(size_t k = 0; k<sampleDefNames.size(); k++){
        if(k==sampleDefNames.size()-1) texFile << sampleDefNames[k];
        else texFile << sampleDefNames[k] << " & ";
    }
    texFile << "\\\\ \\hline" << std::endl;

    // Write the data rows
    for (size_t i = 0; i < histogramCounts.size(); ++i) {
        texFile << "$ {\\rm " << cutDefs[i].cutLabel << "}$" << " & ";
        //texFile << cutDefs[i].varLabel << " & ";
        for (size_t j = 0; j < histogramCounts[i].size(); ++j) {
            std::ostringstream streamObjEff;
            streamObjEff <<  std::fixed << std::setprecision(2) << 100.*histogramCounts[i][j]/histogramCounts[1][j];
            std::cout<<"Efficiency: "<<streamObjEff.str()<<" %"<<std::endl;
            // Get string from out
            if(j==sampleDefNames.size()-1) texFile << histogramCounts[i][j] <<" ("<< streamObjEff.str() <<" \\%)" << std::endl;
            else texFile << histogramCounts[i][j] <<" ("<< streamObjEff.str() <<" \\%)"  << " & ";
        }
        texFile << "\\\\ \\hline" << std::endl;
    }
 
    texFile << "\\end{tabular}" << std::endl;
    texFile << "\\caption{" << tableCaption << "}" << std::endl;
    texFile << "\\end{table}" << std::endl;

    texFile << "\\end{landscape}" << std::endl;

    // Close the LaTeX document
    texFile << "\\end{document}" << std::endl;

    // Close the file
    texFile.close();

    // Compile the .tex file
    std::string compileCommand = "pdflatex " + fileName;
    int compileResult = std::system(compileCommand.c_str());

    if (compileResult != 0) {
        std::cerr << "Failed to compile the .tex file." << std::endl;
    }
}



void SetGStyle(){
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
  gStyle->SetLineWidth(2);
  // Draw horizontal and vertical grids
  gStyle->SetPadGridX(kTRUE);                     
  gStyle->SetPadGridY(kTRUE);

  // set left margin to 10%
  gStyle->SetPadLeftMargin(0.15);

  //gStyle->SetOptLogy(1);
}