// First create histogram ROOT file then run this script.
// First the validation script needs to be loaded:
//   .L Validation.C
//   .x comparePlots.C("output/Validation1.root", "output/Histos1.root", "label1",
//                     "output/Validation2.root", "output/Histos2.root", "label2", runValidation)

#include <cmath>
#include <iomanip>
#include <sstream>
#include <utility>
#include <vector>

// https://root.cern.ch/doc/master/RooCrystalBall_8cxx_source.html

double evaluateCrystalBallTail(double t, double alpha, double n)
{
    double a = std::pow(n / alpha, n) * std::exp(-0.5 * alpha * alpha);
    double b = (n / alpha) - alpha;

    return a / std::pow(b - t, n);
}

double RooCrystalBallFun(double *xArr, double *pars)
{
    const double x = xArr[0];
    const double norm = pars[0];
    const double x0 = pars[1];
    const double sigmaL = std::abs(pars[2]);
    const double sigmaR = std::abs(pars[3]);
    double alphaL = std::abs(pars[4]);
    double nL = pars[5];
    double alphaR = std::abs(pars[6]);
    double nR = pars[7];

    // If alphaL is negative, then the tail will be on the right side.
    // We follow the convention established by RooCBShape.
    //if(!alphaR && alphaL < 0.0) {
    //    std::swap(alphaL, alphaR);
    //    std::swap(nL, nR);
    //}

    const double t = (x - x0) / (x < x0 ? sigmaL : sigmaR);
    double value(0.0);

    if (t < -alphaL)
    {
        value = evaluateCrystalBallTail(t, alphaL, nL);
    }
    else if (t <= alphaR)
    {
        value = std::exp(-0.5 * t * t);
    }
    else
    {
        value = evaluateCrystalBallTail(-t, alphaR, nR);
    }

    return value * norm;
}

void FitCBFun(TH1 *hist, double theMean = -999.0, double theSigma1 = -999.0, double theSigma2 = -999.0, int color = 1)
{
    TString histName = hist->GetName();
    TString funName = histName + "Fun";
    TAxis *xAxis = hist->GetXaxis();
    const double xMin = xAxis->GetXmin();
    const double xMax = xAxis->GetXmax();
    const double mean = theMean > -999.0 ? theMean : hist->GetMean();
    const double sigmaL = theSigma1 > -999.0 ? theSigma1 : hist->GetStdDev();
    const double sigmaR = theSigma2 > -999.0 ? theSigma2 : hist->GetStdDev();
    const double norm = hist->GetMaximum();

    TF1 *fun = new TF1(funName.Data(), &RooCrystalBallFun, xMin, xMax, 8);
    fun->SetLineWidth(2);
    fun->SetLineColor(color);
    fun->SetParameters(norm, mean, sigmaL, sigmaR, 0.1, 1.0, 0.1, 1.0);
    fun->SetParNames("N", "#mu", "#sigma_{L}", "#sigma_{R}", "#alpha_{L}", "n_{L}", "#alpha_{R}", "n_{R}");
    hist->Fit(funName.Data());
}

TString PrintCBPars(TH1 *hist, const TString &preamble)
{
    TString line(preamble);

    if (hist)
    {

        const TString histName = hist->GetName();
        TString funName = histName + "Fun";
        TF1 *fun = hist->GetFunction(funName.Data());
        if (!fun)
        {
            return line;
        }

        const double mean = fun->GetParameter(1);
        const double sigmaL = std::abs(fun->GetParameter(2));
        const double sigmaR = std::abs(fun->GetParameter(3));

        // Convert cm to mm
        const double cm2mm = 10.0;

        std::ostringstream os;
        os << std::fixed << std::setprecision(1) << ": #mu = " << mean * cm2mm << ", #sigma_{L} = " << sigmaL * cm2mm
           << ", #sigma_{R} = " << sigmaR * cm2mm << " mm" << std::endl;
        line += os.str();
    }

    return line;
}

void SetHistColorStyle(TH1 *hist, int color, int style = -1)
{
    hist->SetLineColor(color);
    hist->SetMarkerColor(color);
    hist->SetLineWidth(2);
    if (style >= 0)
    {
        hist->SetMarkerStyle(style);
    }
}

void NormHist(TH1 *hist, TString yLabel)
{
    const double norm = hist->Integral();
    if (norm > 0.0)
    {
        hist->Scale(1.0 / norm);
    }
    hist->GetYaxis()->SetTitle(yLabel.Data());
}

TGraphErrors *getEfficiency(TTree *theTree, int interactionCode, int minProtons, int maxProtons, int color)
{

    std::cout << "getEfficiency for interactionCode = " << interactionCode << " with minProtons = " << minProtons
              << " and maxProtons = " << maxProtons << std::endl;

    int interactionType(0), isCorrectNu(0);
    theTree->SetBranchAddress("interactionType", &interactionType);
    theTree->SetBranchAddress("isCorrectNu", &isCorrectNu);

    // Number of passed & total events per number of protons for the given interaction
    // Number of graph points (at least 1)
    const int NPoints = maxProtons - minProtons + 1;
    // Vectors to keep track of the number of events for the interactionCode
    std::vector<double> passed(NPoints);
    std::vector<double> total(NPoints);
    // Vector of the number of protons for each interaction submode
    std::vector<double> nProtons(NPoints);
    for (int i = 0; i < NPoints; i++)
    {
        nProtons[i] = (minProtons + i) * 1.0;
    }

    const int minCode = interactionCode;
    const int maxCode = interactionCode + NPoints - 1;
    const int N = theTree->GetEntries();

    for (int i = 0; i < N; i++)
    {

        theTree->GetEntry(i);

        // Get required interaction mode
        if (interactionType >= minCode && interactionType <= maxCode)
        {

            const int pointInt = interactionType - interactionCode;
            total[pointInt] += 1.0;

            // Event is correct
            if (isCorrectNu == 1)
            {
                passed[pointInt] += 1.0;
            }
        }
    }

    // Create graph: eff vs N protons
    TGraphErrors *graph = new TGraphErrors(NPoints);
    graph->SetMarkerStyle(kFullCircle);
    graph->SetMarkerColor(color);
    graph->SetLineColor(color);
    graph->SetLineWidth(2);
    TString graphName("Graph");
    graphName += interactionCode;
    graph->SetName(graphName.Data());

    for (int iP = 0; iP < NPoints; iP++)
    {
        double NP = nProtons[iP];
        double eff = total[iP] > 0.0 ? passed[iP] / total[iP] : 0.0;
        double effErr = total[iP] > 0.0 ? sqrt(eff * (1.0 - eff) / total[iP]) : 0.0;
        graph->SetPoint(iP, NP, eff * 100.0);
        graph->SetPointError(iP, 0.0, effErr * 100.0);
    }

    return graph;
}

void plotCorrectFrac(const std::string &validFileName1, const std::string &label1, const std::string &validFileName2, const std::string &label2)
{

    std::cout << "Plotting correct event fraction vs number of protons" << std::endl;
    gStyle->SetGridWidth(1);
    gStyle->SetGridStyle(kDotted);
    TCanvas *theCanvasC = new TCanvas("theCanvasC", "", 1400, 700);
    theCanvasC->UseCurrentStyle();
    theCanvasC->Clear();
    theCanvasC->Divide(2, 2);

    TFile *validFile1 = TFile::Open(validFileName1.c_str(), "read");
    TTree *theTree1 = dynamic_cast<TTree *>(validFile1->Get("Validation"));
    TFile *validFile2 = TFile::Open(validFileName2.c_str(), "read");
    TTree *theTree2 = dynamic_cast<TTree *>(validFile2->Get("Validation"));

    // Interaction enum starting numbers; extra protons add +1 (0 to 5); enum = code line number - 20
    //https://github.com/PandoraPFA/LArContent/blob/master/larpandoracontent/LArHelpers/LArInteractionTypeHelper.h
    const int CCQEL_MU = 0;          // mu p QEL
    const int CCDIS_MU = 94;         // mu p DIS
    const int CCDIS_MU_PIPLUS = 100; // mu p pi
    const int CCDIS_MU_PIZERO = 112; // mu p pi0
    const int NCDIS_P = 118;         // p NC DIS
    const int NCDIS_PIPLUS = 123;    // pi+ NC DIS
    const int NCDIS_PIMINUS = 129;   // pi- ND DIS
    const int NCDIS_PIZERO = 141;    // pi0 ND DIS

    TLatex *latex = new TLatex();
    latex->SetTextSize(0.055);

    TH2F *nullHist = new TH2F("nullHist", "", 2, -1.0, 6.0, 2, 0.0, 100.0);
    nullHist->SetXTitle("Number of protons");
    nullHist->SetYTitle("Correct events (%)");
    nullHist->SetNdivisions(7, "X");

    TLegend nLeg(0.7, 0.625, 0.85, 0.825, "");
    nLeg.SetBorderSize(0);
    nLeg.SetTextSize(0.05);

    TGraphErrors *eff_CCQEL_MU1 = getEfficiency(theTree1, CCQEL_MU, 0, 5, kRed);
    TGraphErrors *eff_CCQEL_MU2 = getEfficiency(theTree2, CCQEL_MU, 0, 5, kBlue);
    theCanvasC->cd(1);
    gPad->SetGrid();
    nullHist->Draw();
    eff_CCQEL_MU1->Draw("samepl");
    eff_CCQEL_MU2->Draw("samepl");
    latex->DrawLatex(3.05, 80.0, "CCQE");
    latex->DrawLatex(3.05, 70.0, "#mu Np");
    nLeg.AddEntry(eff_CCQEL_MU1, label1.c_str());
    nLeg.AddEntry(eff_CCQEL_MU2, label2.c_str());
    nLeg.Draw();

    TGraphErrors *eff_CCDIS_MU1 = getEfficiency(theTree1, CCDIS_MU, 0, 5, kRed);
    TGraphErrors *eff_CCDIS_MU2 = getEfficiency(theTree2, CCDIS_MU, 0, 5, kBlue);
    theCanvasC->cd(2);
    gPad->SetGrid();
    nullHist->Draw();
    eff_CCDIS_MU1->Draw("samepl");
    eff_CCDIS_MU2->Draw("samepl");
    latex->DrawLatex(3.05, 80.0, "CCDIS");
    latex->DrawLatex(3.05, 70.0, "#mu Np");
    nLeg.Draw();

    TGraphErrors *eff_CCDIS_MU_PIPLUS1 = getEfficiency(theTree1, CCDIS_MU_PIPLUS, 0, 5, kRed);
    TGraphErrors *eff_CCDIS_MU_PIPLUS2 = getEfficiency(theTree2, CCDIS_MU_PIPLUS, 0, 5, kBlue);
    theCanvasC->cd(3);
    gPad->SetGrid();
    nullHist->Draw();
    eff_CCDIS_MU_PIPLUS1->Draw("samepl");
    eff_CCDIS_MU_PIPLUS2->Draw("samepl");
    latex->DrawLatex(3.05, 80.0, "CCDIS");
    latex->DrawLatex(3.05, 70.0, "#mu #pi^{#pm} Np");
    nLeg.Draw();

    TGraphErrors *eff_CCDIS_MU_PIZERO1 = getEfficiency(theTree1, CCDIS_MU_PIZERO, 0, 5, kRed);
    TGraphErrors *eff_CCDIS_MU_PIZERO2 = getEfficiency(theTree2, CCDIS_MU_PIZERO, 0, 5, kBlue);
    theCanvasC->cd(4);
    gPad->SetGrid();
    nullHist->Draw();
    eff_CCDIS_MU_PIZERO1->Draw("samepl");
    eff_CCDIS_MU_PIZERO2->Draw("samepl");
    latex->DrawLatex(3.05, 80.0, "CCDIS");
    latex->DrawLatex(3.05, 70.0, "#mu #pi^{0} Np");
    nLeg.Draw();

    theCanvasC->Print("CCEventFrac.png");
    theCanvasC->Print("CCEventFrac.pdf");

    // NC
    TGraphErrors *eff_NCDIS_P1 = getEfficiency(theTree1, NCDIS_P, 1, 5, kRed);
    TGraphErrors *eff_NCDIS_P2 = getEfficiency(theTree2, NCDIS_P, 1, 5, kBlue);
    theCanvasC->cd(1);
    gPad->SetGrid();
    nullHist->Draw();
    eff_NCDIS_P1->Draw("samepl");
    eff_NCDIS_P2->Draw("samepl");
    latex->DrawLatex(3.05, 80.0, "NCDIS");
    latex->DrawLatex(3.05, 70.0, "Np");
    nLeg.Draw();

    TGraphErrors *eff_NCDIS_PIPLUS1 = getEfficiency(theTree1, NCDIS_PIPLUS, 0, 5, kRed);
    TGraphErrors *eff_NCDIS_PIPLUS2 = getEfficiency(theTree2, NCDIS_PIPLUS, 0, 5, kBlue);
    theCanvasC->cd(2);
    gPad->SetGrid();
    nullHist->Draw();
    eff_NCDIS_PIPLUS1->Draw("samepl");
    eff_NCDIS_PIPLUS2->Draw("samepl");
    latex->DrawLatex(3.05, 80.0, "NCDIS");
    latex->DrawLatex(3.05, 70.0, "#pi^{+} Np");
    nLeg.Draw();

    TGraphErrors *eff_NCDIS_PIMINUS1 = getEfficiency(theTree1, NCDIS_PIMINUS, 0, 5, kRed);
    TGraphErrors *eff_NCDIS_PIMINUS2 = getEfficiency(theTree2, NCDIS_PIMINUS, 0, 5, kBlue);
    theCanvasC->cd(3);
    gPad->SetGrid();
    nullHist->Draw();
    eff_NCDIS_PIMINUS1->Draw("samepl");
    eff_NCDIS_PIMINUS2->Draw("samepl");
    latex->DrawLatex(3.05, 80.0, "NCDIS");
    latex->DrawLatex(3.05, 70.0, "#pi^{-} Np");
    nLeg.Draw();

    TGraphErrors *eff_NCDIS_PIZERO1 = getEfficiency(theTree1, NCDIS_PIZERO, 0, 5, kRed);
    TGraphErrors *eff_NCDIS_PIZERO2 = getEfficiency(theTree2, NCDIS_PIZERO, 0, 5, kBlue);
    theCanvasC->cd(4);
    gPad->SetGrid();
    nullHist->Draw();
    eff_NCDIS_PIZERO1->Draw("samepl");
    eff_NCDIS_PIZERO2->Draw("samepl");
    latex->DrawLatex(3.05, 80.0, "NCDIS");
    latex->DrawLatex(3.05, 70.0, "#pi^{0} Np");
    nLeg.Draw();

    theCanvasC->Print("NCEventFrac.png");
    theCanvasC->Print("NCEventFrac.pdf");

    validFile1->Close();
    validFile2->Close();
}

void comparePlots(const std::string &validFileName1 = "ValidationDLVtx.root", const std::string &histFileName1 = "HistosDLVtx.root",
    const std::string &label1 = "DLVtx", const std::string &validFileName2 = "Validation3D_DLVtx.root",
    const std::string &histFileName2 = "Histos3D_DLVtx.root", const std::string &label2 = "3D_DLVtx", bool runValidation = false)
{

    gROOT->SetStyle("Plain");
    gStyle->SetOptStat(0);
    gStyle->SetGridWidth(1);
    gStyle->SetGridStyle(kDotted);

    TCanvas *theCanvas = new TCanvas("theCanvas", "", 1800, 1400);
    theCanvas->UseCurrentStyle();

    // Run validation code to create histograms
    Parameters p;
    p.m_histogramOutput = true;
    p.m_vertexXCorrection = 0.0;
    // Reduce verbosity
    p.m_displayMatchedEvents = false;
    if (runValidation)
    {
        p.m_histFileName = histFileName1;
        Validation(validFileName1, p);
        p.m_histFileName = histFileName2;
        Validation(validFileName2, p);
    }

    // Plot histograms
    theCanvas->Divide(2, 2);

    TLatex *text = new TLatex();
    text->SetTextSize(0.075);
    text->SetNDC(true);

    TFile *histFile1 = TFile::Open(histFileName1.c_str(), "read");
    TFile *histFile2 = TFile::Open(histFileName2.c_str(), "read");

    // ALL HITS EFFICIENCY
    // maxMerged voxel limit
    double xMax = 10000;

    TLegend hitsEffLeg(0.715, 0.2, 0.825, 0.4, "");
    hitsEffLeg.SetBorderSize(0);
    hitsEffLeg.SetTextSize(0.05);

    theCanvas->cd(1);
    gPad->SetGrid();
    TH1 *muHitsEff1 = dynamic_cast<TH1 *>(histFile1->Get("ALL_INTERACTIONS_MUON_HitsEfficiency"));
    TH1 *muHitsEff2 = dynamic_cast<TH1 *>(histFile2->Get("ALL_INTERACTIONS_MUON_HitsEfficiency"));
    if (muHitsEff1 && muHitsEff2)
    {
        muHitsEff1->GetXaxis()->SetRangeUser(0.0, xMax);
        muHitsEff2->GetXaxis()->SetRangeUser(0.0, xMax);
        gPad->SetLogx();
        SetHistColorStyle(muHitsEff1, kRed);
        SetHistColorStyle(muHitsEff2, kBlue);
        muHitsEff1->Draw();
        muHitsEff2->Draw("same");
        text->DrawLatex(0.775, 0.45, "#mu");
        hitsEffLeg.AddEntry(muHitsEff1, label1.c_str());
        hitsEffLeg.AddEntry(muHitsEff2, label2.c_str());
        hitsEffLeg.Draw();
    }

    theCanvas->cd(2);
    gPad->SetGrid();
    TH1 *pHitsEff1 = dynamic_cast<TH1 *>(histFile1->Get("ALL_INTERACTIONS_PROTON1_HitsEfficiency"));
    TH1 *pHitsEff2 = dynamic_cast<TH1 *>(histFile2->Get("ALL_INTERACTIONS_PROTON1_HitsEfficiency"));
    if (pHitsEff1 && pHitsEff2)
    {
        pHitsEff1->GetXaxis()->SetRangeUser(0.0, xMax);
        pHitsEff2->GetXaxis()->SetRangeUser(0.0, xMax);
        gPad->SetLogx();
        SetHistColorStyle(pHitsEff1, kRed);
        SetHistColorStyle(pHitsEff2, kBlue);
        pHitsEff1->Draw();
        pHitsEff2->Draw("same");
        text->DrawLatex(0.775, 0.45, "p");
        hitsEffLeg.Draw();
    }

    theCanvas->cd(3);
    gPad->SetGrid();
    TH1 *pipHitsEff1 = dynamic_cast<TH1 *>(histFile1->Get("ALL_INTERACTIONS_PIPLUS_HitsEfficiency"));
    TH1 *pipHitsEff2 = dynamic_cast<TH1 *>(histFile2->Get("ALL_INTERACTIONS_PIPLUS_HitsEfficiency"));
    if (pipHitsEff1 && pipHitsEff2)
    {
        pipHitsEff1->GetXaxis()->SetRangeUser(0.0, xMax);
        pipHitsEff2->GetXaxis()->SetRangeUser(0.0, xMax);
        gPad->SetLogx();
        SetHistColorStyle(pipHitsEff1, kRed);
        SetHistColorStyle(pipHitsEff2, kBlue);
        pipHitsEff1->Draw();
        pipHitsEff2->Draw("same");
        text->DrawLatex(0.775, 0.45, "#pi^{+}");
        hitsEffLeg.Draw();
    }

    theCanvas->cd(4);
    gPad->SetGrid();
    TH1 *pimHitsEff1 = dynamic_cast<TH1 *>(histFile1->Get("ALL_INTERACTIONS_PIMINUS_HitsEfficiency"));
    TH1 *pimHitsEff2 = dynamic_cast<TH1 *>(histFile2->Get("ALL_INTERACTIONS_PIMINUS_HitsEfficiency"));
    if (pimHitsEff1 && pimHitsEff2)
    {
        pimHitsEff1->GetXaxis()->SetRangeUser(0.0, xMax);
        pimHitsEff2->GetXaxis()->SetRangeUser(0.0, xMax);
        gPad->SetLogx();
        SetHistColorStyle(pimHitsEff1, kRed);
        SetHistColorStyle(pimHitsEff2, kBlue);
        pimHitsEff1->Draw();
        pimHitsEff2->Draw("same");
        text->DrawLatex(0.775, 0.45, "#pi^{-}");
        hitsEffLeg.Draw();
    }
    theCanvas->Print("allHitsEff.png");
    theCanvas->Print("allHitsEff.pdf");

    // ALL MOMENTUM EFFICIENCY
    theCanvas->cd(1);
    gPad->SetLogx(0);
    gPad->SetGrid();
    TH1 *muMtmEff1 = dynamic_cast<TH1 *>(histFile1->Get("ALL_INTERACTIONS_MUON_MomentumEfficiency"));
    TH1 *muMtmEff2 = dynamic_cast<TH1 *>(histFile2->Get("ALL_INTERACTIONS_MUON_MomentumEfficiency"));
    if (muMtmEff1 && muMtmEff2)
    {
        SetHistColorStyle(muMtmEff1, kRed);
        SetHistColorStyle(muMtmEff2, kBlue);
        muMtmEff1->Draw();
        muMtmEff2->Draw("same");
        text->DrawLatex(0.775, 0.45, "#mu");
        hitsEffLeg.Draw();
    }

    theCanvas->cd(2);
    gPad->SetLogx(0);
    gPad->SetGrid();
    TH1 *pMtmEff1 = dynamic_cast<TH1 *>(histFile1->Get("ALL_INTERACTIONS_PROTON1_MomentumEfficiency"));
    TH1 *pMtmEff2 = dynamic_cast<TH1 *>(histFile2->Get("ALL_INTERACTIONS_PROTON1_MomentumEfficiency"));
    if (pMtmEff1 && pMtmEff2)
    {
        SetHistColorStyle(pMtmEff1, kRed);
        SetHistColorStyle(pMtmEff2, kBlue);
        pMtmEff1->Draw();
        pMtmEff2->Draw("same");
        text->DrawLatex(0.775, 0.45, "p");
        hitsEffLeg.Draw();
    }

    theCanvas->cd(3);
    gPad->SetLogx(0);
    gPad->SetGrid();
    TH1 *pipMtmEff1 = dynamic_cast<TH1 *>(histFile1->Get("ALL_INTERACTIONS_PIPLUS_MomentumEfficiency"));
    TH1 *pipMtmEff2 = dynamic_cast<TH1 *>(histFile2->Get("ALL_INTERACTIONS_PIPLUS_MomentumEfficiency"));
    if (pipMtmEff1 && pipMtmEff2)
    {
        SetHistColorStyle(pipMtmEff1, kRed);
        SetHistColorStyle(pipMtmEff2, kBlue);
        pipMtmEff1->Draw();
        pipMtmEff2->Draw("same");
        text->DrawLatex(0.775, 0.45, "#pi^{+}");
        hitsEffLeg.Draw();
    }

    theCanvas->cd(4);
    gPad->SetLogx(0);
    gPad->SetGrid();
    TH1 *pimMtmEff1 = dynamic_cast<TH1 *>(histFile1->Get("ALL_INTERACTIONS_PIMINUS_MomentumEfficiency"));
    TH1 *pimMtmEff2 = dynamic_cast<TH1 *>(histFile2->Get("ALL_INTERACTIONS_PIMINUS_MomentumEfficiency"));
    if (pimMtmEff1 && pimMtmEff2)
    {
        SetHistColorStyle(pimMtmEff1, kRed);
        SetHistColorStyle(pimMtmEff2, kBlue);
        pimMtmEff1->Draw();
        pimMtmEff2->Draw("same");
        text->DrawLatex(0.775, 0.45, "#pi^{-}");
        hitsEffLeg.Draw();
    }

    theCanvas->Print("allMtmEff.png");
    theCanvas->Print("allMtmEff.pdf");

    // ALL VERTICES
    gStyle->SetOptFit(0);
    gROOT->ForceStyle();
    theCanvas->UseCurrentStyle();
    theCanvas->Update();

    TLatex *latex = new TLatex();
    latex->SetTextSize(0.045);

    theCanvas->cd(1);
    gPad->SetGrid();
    TH1 *xVtxHist1 = dynamic_cast<TH1 *>(histFile1->Get("ALL_INTERACTIONS_VtxDeltaX"));
    TH1 *xVtxHist2 = dynamic_cast<TH1 *>(histFile2->Get("ALL_INTERACTIONS_VtxDeltaX"));
    NormHist(xVtxHist1, "Fraction of events");
    NormHist(xVtxHist2, "Fraction of events");
    //FitCBFun(xVtxHist1, xVtxHist1->GetMean(), xVtxHist1->GetStdDev(), xVtxHist1->GetStdDev(), kRed);
    //FitCBFun(xVtxHist2, xVtxHist2->GetMean(), xVtxHist2->GetStdDev(), xVtxHist2->GetStdDev(), kBlue);
    SetHistColorStyle(xVtxHist1, kRed);
    SetHistColorStyle(xVtxHist2, kBlue);
    xVtxHist1->GetYaxis()->SetTitleOffset(1.35);
    xVtxHist1->GetYaxis()->SetLabelSize(0.0325);

    if (xVtxHist1->GetMaximum() > xVtxHist2->GetMaximum())
    {
        xVtxHist1->Draw();
        xVtxHist2->Draw("same");
    }
    else
    {
        xVtxHist2->Draw();
        xVtxHist1->Draw("same");
    }

    // Print parameters
    latex->SetTextColor(kRed);
    // 0.535, 0.75 then 0.65
    //latex->SetTextSize(0.0325);
    //latex->DrawLatexNDC(0.1, 0.925, PrintCBPars(xVtxHist1, label1.c_str()));
    latex->DrawLatexNDC(0.7, 0.7, PrintCBPars(xVtxHist1, label1.c_str()));
    latex->SetTextColor(kBlue);
    //latex->DrawLatexNDC(0.525, 0.925, PrintCBPars(xVtxHist2, label2.c_str()));
    latex->DrawLatexNDC(0.7, 0.6, PrintCBPars(xVtxHist2, label2.c_str()));
    gPad->Update();

    theCanvas->cd(2);
    gPad->SetGrid();
    TH1 *yVtxHist1 = dynamic_cast<TH1 *>(histFile1->Get("ALL_INTERACTIONS_VtxDeltaY"));
    TH1 *yVtxHist2 = dynamic_cast<TH1 *>(histFile2->Get("ALL_INTERACTIONS_VtxDeltaY"));
    NormHist(yVtxHist1, "Fraction of events");
    NormHist(yVtxHist2, "Fraction of events");
    //FitCBFun(yVtxHist1, yVtxHist1->GetMean(), yVtxHist1->GetStdDev(), yVtxHist1->GetStdDev(), kRed);
    //FitCBFun(yVtxHist2, yVtxHist2->GetMean(), yVtxHist2->GetStdDev(), yVtxHist2->GetStdDev(), kBlue);
    SetHistColorStyle(yVtxHist1, kRed);
    SetHistColorStyle(yVtxHist2, kBlue);
    yVtxHist1->GetYaxis()->SetTitleOffset(1.35);
    yVtxHist1->GetYaxis()->SetLabelSize(0.0325);

    if (yVtxHist1->GetMaximum() > yVtxHist2->GetMaximum())
    {
        yVtxHist1->Draw();
        yVtxHist2->Draw("same");
    }
    else
    {
        yVtxHist2->Draw();
        yVtxHist1->Draw("same");
    }

    // Print parameters
    latex->SetTextColor(kRed);
    //latex->DrawLatexNDC(0.1, 0.925, PrintCBPars(yVtxHist1, label1.c_str()));
    latex->DrawLatexNDC(0.7, 0.7, PrintCBPars(yVtxHist1, label1.c_str()));
    latex->SetTextColor(kBlue);
    //latex->DrawLatexNDC(0.525, 0.925, PrintCBPars(yVtxHist2, label2.c_str()));
    latex->DrawLatexNDC(0.7, 0.6, PrintCBPars(yVtxHist2, label2.c_str()));
    gPad->Update();

    theCanvas->cd(3);
    gPad->SetGrid();
    TH1 *zVtxHist1 = dynamic_cast<TH1 *>(histFile1->Get("ALL_INTERACTIONS_VtxDeltaZ"));
    TH1 *zVtxHist2 = dynamic_cast<TH1 *>(histFile2->Get("ALL_INTERACTIONS_VtxDeltaZ"));
    NormHist(zVtxHist1, "Fraction of events");
    NormHist(zVtxHist2, "Fraction of events");
    //FitCBFun(zVtxHist1, zVtxHist1->GetMean(), zVtxHist1->GetStdDev(), zVtxHist1->GetStdDev(), kRed);
    //FitCBFun(zVtxHist2, zVtxHist2->GetMean(), zVtxHist2->GetStdDev(), zVtxHist2->GetStdDev(), kBlue);
    SetHistColorStyle(zVtxHist1, kRed);
    SetHistColorStyle(zVtxHist2, kBlue);
    zVtxHist1->GetYaxis()->SetTitleOffset(1.35);
    zVtxHist1->GetYaxis()->SetLabelSize(0.0325);
    if (zVtxHist1->GetMaximum() > zVtxHist2->GetMaximum())
    {
        zVtxHist1->Draw();
        zVtxHist2->Draw("same");
    }
    else
    {
        zVtxHist2->Draw();
        zVtxHist1->Draw("same");
    }

    // Print parameters
    latex->SetTextColor(kRed);
    //latex->DrawLatexNDC(0.1, 0.925, PrintCBPars(zVtxHist1, label1.c_str()));
    latex->DrawLatexNDC(0.7, 0.7, PrintCBPars(zVtxHist1, label1.c_str()));
    latex->SetTextColor(kBlue);
    //latex->DrawLatexNDC(0.525, 0.925, PrintCBPars(zVtxHist2, label2.c_str()));
    latex->DrawLatexNDC(0.7, 0.6, PrintCBPars(zVtxHist2, label2.c_str()));
    gPad->Update();

    theCanvas->cd(4);
    gPad->SetGrid();
    TH1 *rVtxHist1 = dynamic_cast<TH1 *>(histFile1->Get("ALL_INTERACTIONS_VtxDeltaR"));
    TH1 *rVtxHist2 = dynamic_cast<TH1 *>(histFile2->Get("ALL_INTERACTIONS_VtxDeltaR"));
    NormHist(rVtxHist1, "Fraction of events");
    NormHist(rVtxHist2, "Fraction of events");
    //FitCBFun(rVtxHist1, rVtxHist1->GetMean(), rVtxHist1->GetStdDev(), rVtxHist1->GetStdDev(), kRed);
    //FitCBFun(rVtxHist2, rVtxHist2->GetMean(), rVtxHist2->GetStdDev(), rVtxHist2->GetStdDev(), kBlue);
    SetHistColorStyle(rVtxHist1, kRed);
    SetHistColorStyle(rVtxHist2, kBlue);
    rVtxHist1->GetYaxis()->SetTitleOffset(1.35);
    rVtxHist1->GetYaxis()->SetLabelSize(0.0325);
    if (rVtxHist1->GetMaximum() > rVtxHist2->GetMaximum())
    {
        rVtxHist1->Draw();
        rVtxHist2->Draw("same");
    }
    else
    {
        rVtxHist2->Draw();
        rVtxHist1->Draw("same");
    }
    // Print parameters
    latex->SetTextColor(kRed);
    //latex->DrawLatexNDC(0.1, 0.925, PrintCBPars(rVtxHist1, label1.c_str()));
    latex->DrawLatexNDC(0.7, 0.7, PrintCBPars(rVtxHist1, label1.c_str()));
    latex->SetTextColor(kBlue);
    //latex->DrawLatexNDC(0.525, 0.925, PrintCBPars(rVtxHist2, label2.c_str()));
    latex->DrawLatexNDC(0.7, 0.6, PrintCBPars(rVtxHist2, label2.c_str()));
    gPad->Update();

    theCanvas->Print("allVtx.png");
    theCanvas->Print("allVtx.pdf");

    // Log scale
    gStyle->SetOptFit(0);
    gROOT->ForceStyle();
    theCanvas->UseCurrentStyle();
    theCanvas->Update();

    theCanvas->cd(1);
    gPad->SetGrid();
    gPad->SetLogy();
    SetHistColorStyle(xVtxHist1, kRed);
    SetHistColorStyle(xVtxHist2, kBlue);
    xVtxHist1->Draw();
    xVtxHist2->Draw("same");
    // Print parameters
    latex->SetTextColor(kRed);
    latex->DrawLatexNDC(0.7, 0.7, PrintCBPars(xVtxHist1, label1.c_str()));
    //latex->DrawLatexNDC(0.1, 0.925, PrintCBPars(xVtxHist1, label1.c_str()));
    latex->SetTextColor(kBlue);
    //latex->DrawLatexNDC(0.525, 0.925, PrintCBPars(xVtxHist2, label2.c_str()));
    latex->DrawLatexNDC(0.7, 0.6, PrintCBPars(xVtxHist2, label2.c_str()));
    gPad->Update();

    theCanvas->cd(2);
    gPad->SetGrid();
    gPad->SetLogy();
    SetHistColorStyle(yVtxHist1, kRed);
    SetHistColorStyle(yVtxHist2, kBlue);
    yVtxHist1->Draw();
    yVtxHist2->Draw("same");
    // Print parameters
    latex->SetTextColor(kRed);
    latex->DrawLatexNDC(0.7, 0.7, PrintCBPars(yVtxHist1, label1.c_str()));
    //latex->DrawLatexNDC(0.1, 0.925, PrintCBPars(yVtxHist1, label1.c_str()));
    latex->SetTextColor(kBlue);
    //latex->DrawLatexNDC(0.525, 0.925, PrintCBPars(yVtxHist2, label2.c_str()));
    latex->DrawLatexNDC(0.7, 0.6, PrintCBPars(yVtxHist2, label2.c_str()));
    gPad->Update();

    theCanvas->cd(3);
    gPad->SetGrid();
    gPad->SetLogy();
    SetHistColorStyle(zVtxHist1, kRed);
    SetHistColorStyle(zVtxHist2, kBlue);
    zVtxHist1->Draw();
    zVtxHist2->Draw("same");
    // Print parameters
    latex->SetTextColor(kRed);
    //latex->DrawLatexNDC(0.1, 0.925, PrintCBPars(zVtxHist1, label1.c_str()));
    latex->DrawLatexNDC(0.7, 0.7, PrintCBPars(zVtxHist1, label1.c_str()));
    latex->SetTextColor(kBlue);
    //latex->DrawLatexNDC(0.525, 0.925, PrintCBPars(zVtxHist2, label2.c_str()));
    latex->DrawLatexNDC(0.7, 0.6, PrintCBPars(zVtxHist2, label2.c_str()));
    gPad->Update();

    theCanvas->cd(4);
    gPad->SetGrid();
    gPad->SetLogy();
    SetHistColorStyle(rVtxHist1, kRed);
    SetHistColorStyle(rVtxHist2, kBlue);
    rVtxHist1->Draw();
    rVtxHist2->Draw("same");
    // Print parameters
    latex->SetTextColor(kRed);
    //latex->DrawLatexNDC(0.1, 0.925, PrintCBPars(rVtxHist1, label1.c_str()));
    latex->DrawLatexNDC(0.7, 0.7, PrintCBPars(rVtxHist1, label1.c_str()));
    latex->SetTextColor(kBlue);
    //latex->DrawLatexNDC(0.525, 0.925, PrintCBPars(rVtxHist2, label2.c_str()));
    latex->DrawLatexNDC(0.7, 0.6, PrintCBPars(rVtxHist2, label2.c_str()));
    gPad->Update();

    theCanvas->Print("logAllVtx.png");
    theCanvas->Print("logAllVtx.pdf");

    delete theCanvas;

    // CCQE: numu + Ar -> mu + p
    TCanvas *theCanvas2 = new TCanvas("theCanvas2", "", 1400, 700);
    theCanvas2->UseCurrentStyle();
    theCanvas2->Clear();
    theCanvas2->Divide(3, 2);

    theCanvas2->cd(1);
    gPad->SetGrid();

    TLegend CCQE_Leg(0.6, 0.2, 0.8, 0.4, "CCQE (#mu p)");
    CCQE_Leg.SetNColumns(2);
    CCQE_Leg.SetBorderSize(0);
    CCQE_Leg.SetTextSize(0.05);

    // Hits efficiency
    TH1 *CCQE_muHitsEff1 = dynamic_cast<TH1 *>(histFile1->Get("CCQEL_MU_P_MUON_HitsEfficiency"));
    TH1 *CCQE_pHitsEff1 = dynamic_cast<TH1 *>(histFile1->Get("CCQEL_MU_P_PROTON1_HitsEfficiency"));
    TH1 *CCQE_muHitsEff2 = dynamic_cast<TH1 *>(histFile2->Get("CCQEL_MU_P_MUON_HitsEfficiency"));
    TH1 *CCQE_pHitsEff2 = dynamic_cast<TH1 *>(histFile2->Get("CCQEL_MU_P_PROTON1_HitsEfficiency"));
    if (CCQE_muHitsEff1 && CCQE_pHitsEff1 && CCQE_muHitsEff2 && CCQE_pHitsEff2)
    {
        CCQE_muHitsEff1->GetXaxis()->SetRangeUser(0.0, xMax);
        gPad->SetLogx();
        SetHistColorStyle(CCQE_muHitsEff1, kRed, kFullCircle);
        CCQE_Leg.AddEntry(CCQE_muHitsEff1, "#mu");
        CCQE_muHitsEff1->Draw();

        CCQE_pHitsEff1->GetXaxis()->SetRangeUser(0.0, xMax);
        SetHistColorStyle(CCQE_pHitsEff1, kRed, kOpenCircle);
        CCQE_Leg.AddEntry(CCQE_pHitsEff1, "p");
        CCQE_pHitsEff1->Draw("same");

        CCQE_muHitsEff2->GetXaxis()->SetRangeUser(0.0, xMax);
        SetHistColorStyle(CCQE_muHitsEff2, kBlue, kFullCircle);
        CCQE_Leg.AddEntry(CCQE_muHitsEff2, "#mu");
        CCQE_muHitsEff2->Draw("same");

        CCQE_pHitsEff2->GetXaxis()->SetRangeUser(0.0, xMax);
        SetHistColorStyle(CCQE_pHitsEff2, kBlue, kOpenCircle);
        CCQE_Leg.AddEntry(CCQE_pHitsEff2, "p");
        CCQE_pHitsEff2->Draw("same");

        CCQE_Leg.Draw();
    }

    // Momentum efficiency
    theCanvas2->cd(2);
    gPad->SetGrid();
    TH1 *CCQE_muMtmEff1 = dynamic_cast<TH1 *>(histFile1->Get("CCQEL_MU_P_MUON_MomentumEfficiency"));
    TH1 *CCQE_pMtmEff1 = dynamic_cast<TH1 *>(histFile1->Get("CCQEL_MU_P_PROTON1_MomentumEfficiency"));
    TH1 *CCQE_muMtmEff2 = dynamic_cast<TH1 *>(histFile2->Get("CCQEL_MU_P_MUON_MomentumEfficiency"));
    TH1 *CCQE_pMtmEff2 = dynamic_cast<TH1 *>(histFile2->Get("CCQEL_MU_P_PROTON1_MomentumEfficiency"));
    if (CCQE_muMtmEff1 && CCQE_pMtmEff1 && CCQE_muMtmEff2 && CCQE_pMtmEff2)
    {
        SetHistColorStyle(CCQE_muMtmEff1, kRed, kFullCircle);
        CCQE_muMtmEff1->Draw();

        SetHistColorStyle(CCQE_pMtmEff1, kRed, kOpenCircle);
        CCQE_pMtmEff1->Draw("same");

        SetHistColorStyle(CCQE_muMtmEff2, kBlue, kFullCircle);
        CCQE_muMtmEff2->Draw("same");

        SetHistColorStyle(CCQE_pMtmEff2, kBlue, kOpenCircle);
        CCQE_pMtmEff2->Draw("same");

        CCQE_Leg.Draw();
    }

    // Completeness
    TLegend CCQE_Leg2(0.2, 0.6, 0.4, 0.8, "CCQE (#mu p)");
    CCQE_Leg2.SetNColumns(2);
    CCQE_Leg2.SetBorderSize(0);
    CCQE_Leg2.SetTextSize(0.05);

    theCanvas2->cd(4);
    gPad->SetGrid();
    TH1 *CCQE_muComp1 = dynamic_cast<TH1 *>(histFile1->Get("CCQEL_MU_P_MUON_Completeness"));
    TH1 *CCQE_pComp1 = dynamic_cast<TH1 *>(histFile1->Get("CCQEL_MU_P_PROTON1_Completeness"));
    TH1 *CCQE_muComp2 = dynamic_cast<TH1 *>(histFile2->Get("CCQEL_MU_P_MUON_Completeness"));
    TH1 *CCQE_pComp2 = dynamic_cast<TH1 *>(histFile2->Get("CCQEL_MU_P_PROTON1_Completeness"));
    if (CCQE_muComp1 && CCQE_pComp1 && CCQE_muComp2 && CCQE_pComp2)
    {
        SetHistColorStyle(CCQE_muComp1, kRed, kFullCircle);
        CCQE_muComp1->SetTitleOffset(1.3, "Y");
        CCQE_Leg2.AddEntry(CCQE_muComp1, "#mu");
        CCQE_muComp1->Draw();
        gPad->SetLogy();

        SetHistColorStyle(CCQE_pComp1, kRed, kOpenCircle);
        CCQE_Leg2.AddEntry(CCQE_pComp1, "p");
        CCQE_pComp1->Draw("same");

        SetHistColorStyle(CCQE_muComp2, kBlue, kFullCircle);
        CCQE_Leg2.AddEntry(CCQE_muComp2, "#mu");
        CCQE_muComp2->Draw("same");

        SetHistColorStyle(CCQE_pComp2, kBlue, kOpenCircle);
        CCQE_Leg2.AddEntry(CCQE_pComp2, "p");
        CCQE_pComp2->Draw("same");

        CCQE_Leg2.Draw();
    }

    // Purity
    theCanvas2->cd(5);
    gPad->SetGrid();
    TH1 *CCQE_muPure1 = dynamic_cast<TH1 *>(histFile1->Get("CCQEL_MU_P_MUON_Purity"));
    TH1 *CCQE_pPure1 = dynamic_cast<TH1 *>(histFile1->Get("CCQEL_MU_P_PROTON1_Purity"));
    TH1 *CCQE_muPure2 = dynamic_cast<TH1 *>(histFile2->Get("CCQEL_MU_P_MUON_Purity"));
    TH1 *CCQE_pPure2 = dynamic_cast<TH1 *>(histFile2->Get("CCQEL_MU_P_PROTON1_Purity"));
    if (CCQE_muPure1 && CCQE_pPure1 && CCQE_muPure2 && CCQE_pPure2)
    {
        SetHistColorStyle(CCQE_muPure1, kRed, kFullCircle);
        CCQE_muPure1->SetTitleOffset(1.3, "Y");
        CCQE_muPure1->Draw();
        gPad->SetLogy();

        SetHistColorStyle(CCQE_pPure1, kRed, kOpenCircle);
        CCQE_pPure1->Draw("same");

        SetHistColorStyle(CCQE_muPure2, kBlue, kFullCircle);
        CCQE_muPure2->Draw("same");

        SetHistColorStyle(CCQE_pPure2, kBlue, kOpenCircle);
        CCQE_pPure2->Draw("same");

        CCQE_Leg2.Draw();
    }

    // Vtx dR
    theCanvas2->cd(6);
    gPad->SetGrid();
    TH1 *CCQE_muVtxR1 = dynamic_cast<TH1 *>(histFile1->Get("CCQEL_MU_P_VtxDeltaR"));
    TH1 *CCQE_muVtxR2 = dynamic_cast<TH1 *>(histFile2->Get("CCQEL_MU_P_VtxDeltaR"));
    if (CCQE_muVtxR1 && CCQE_muVtxR2)
    {
        NormHist(CCQE_muVtxR1, "Fraction of events");
        //FitCBFun(CCQE_muVtxR1, 0.4, 0.15, 0.15, kRed);
        SetHistColorStyle(CCQE_muVtxR1, kRed);
        CCQE_muVtxR1->SetTitleOffset(1.3, "Y");

        NormHist(CCQE_muVtxR2, "Fraction of events");
        //FitCBFun(CCQE_muVtxR2, 0.3, 0.3, 0.3, kBlue);
        SetHistColorStyle(CCQE_muVtxR2, kBlue);

        if (CCQE_muVtxR1->GetMaximum() > CCQE_muVtxR2->GetMaximum())
        {
            CCQE_muVtxR1->Draw();
            CCQE_muVtxR2->Draw("same");
        }
        else
        {
            CCQE_muVtxR2->Draw();
            CCQE_muVtxR1->Draw("same");
        }

        // Print parameters
        latex->SetTextSize(0.045);
        latex->SetTextColor(kRed);
        latex->DrawLatexNDC(0.335, 0.75, PrintCBPars(CCQE_muVtxR1, label1.c_str()));
        latex->SetTextColor(kBlue);
        latex->DrawLatexNDC(0.335, 0.65, PrintCBPars(CCQE_muVtxR2, label2.c_str()));
    }

    theCanvas2->Print("CCQE.png");
    theCanvas2->Print("CCQE.pdf");

    // CCDIS: numu + Ar -> mu- p pi+
    theCanvas2->Clear();
    theCanvas2->Divide(3, 2);
    theCanvas2->cd(1);
    gPad->SetGrid();

    TLegend CCDIS_Leg(0.6, 0.2, 0.8, 0.45, "CCDIS (#mu p #pi)");
    CCDIS_Leg.SetNColumns(3);
    CCDIS_Leg.SetBorderSize(0);
    CCDIS_Leg.SetTextSize(0.05);

    // Hits efficiency
    TH1 *CCDIS_muHitsEff1 = dynamic_cast<TH1 *>(histFile1->Get("CCDIS_MU_P_PIPLUS_MUON_HitsEfficiency"));
    TH1 *CCDIS_pHitsEff1 = dynamic_cast<TH1 *>(histFile1->Get("CCDIS_MU_P_PIPLUS_PROTON1_HitsEfficiency"));
    TH1 *CCDIS_piHitsEff1 = dynamic_cast<TH1 *>(histFile1->Get("CCDIS_MU_P_PIPLUS_PIPLUS_HitsEfficiency"));
    TH1 *CCDIS_muHitsEff2 = dynamic_cast<TH1 *>(histFile2->Get("CCDIS_MU_P_PIPLUS_MUON_HitsEfficiency"));
    TH1 *CCDIS_pHitsEff2 = dynamic_cast<TH1 *>(histFile2->Get("CCDIS_MU_P_PIPLUS_PROTON1_HitsEfficiency"));
    TH1 *CCDIS_piHitsEff2 = dynamic_cast<TH1 *>(histFile2->Get("CCDIS_MU_P_PIPLUS_PIPLUS_HitsEfficiency"));
    if (CCDIS_muHitsEff1 && CCDIS_pHitsEff1 && CCDIS_piHitsEff1 && CCDIS_muHitsEff2 && CCDIS_pHitsEff2 && CCDIS_piHitsEff2)
    {

        CCDIS_muHitsEff1->GetXaxis()->SetRangeUser(0.0, xMax);
        gPad->SetLogx();
        SetHistColorStyle(CCDIS_muHitsEff1, kRed, kFullCircle);
        CCDIS_Leg.AddEntry(CCDIS_muHitsEff1, "#mu");
        CCDIS_muHitsEff1->Draw();

        CCDIS_pHitsEff1->GetXaxis()->SetRangeUser(0.0, xMax);
        SetHistColorStyle(CCDIS_pHitsEff1, kRed, kOpenCircle);
        CCDIS_Leg.AddEntry(CCDIS_pHitsEff1, "p");
        CCDIS_pHitsEff1->Draw("same");

        CCDIS_piHitsEff1->GetXaxis()->SetRangeUser(0.0, xMax);
        SetHistColorStyle(CCDIS_piHitsEff1, kOrange + 7, kFullDiamond);
        CCDIS_Leg.AddEntry(CCDIS_piHitsEff1, "#pi");
        CCDIS_piHitsEff1->Draw("same");

        CCDIS_muHitsEff2->GetXaxis()->SetRangeUser(0.0, xMax);
        SetHistColorStyle(CCDIS_muHitsEff2, kBlue, kFullCircle);
        CCDIS_Leg.AddEntry(CCDIS_muHitsEff2, "#mu");
        CCDIS_muHitsEff2->Draw("same");

        CCDIS_pHitsEff2->GetXaxis()->SetRangeUser(0.0, xMax);
        SetHistColorStyle(CCDIS_pHitsEff2, kBlue, kOpenCircle);
        CCDIS_Leg.AddEntry(CCDIS_pHitsEff2, "p");
        CCDIS_pHitsEff2->Draw("same");

        CCDIS_piHitsEff2->GetXaxis()->SetRangeUser(0.0, xMax);
        SetHistColorStyle(CCDIS_piHitsEff2, kGreen + 3, kFullDiamond);
        CCDIS_Leg.AddEntry(CCDIS_piHitsEff2, "#pi");
        CCDIS_piHitsEff2->Draw("same");

        CCDIS_Leg.Draw();
    }

    // Momentum efficiency
    theCanvas2->cd(2);
    gPad->SetGrid();
    TH1 *CCDIS_muMtmEff1 = dynamic_cast<TH1 *>(histFile1->Get("CCDIS_MU_P_PIPLUS_MUON_MomentumEfficiency"));
    TH1 *CCDIS_pMtmEff1 = dynamic_cast<TH1 *>(histFile1->Get("CCDIS_MU_P_PIPLUS_PROTON1_MomentumEfficiency"));
    TH1 *CCDIS_piMtmEff1 = dynamic_cast<TH1 *>(histFile1->Get("CCDIS_MU_P_PIPLUS_PIPLUS_MomentumEfficiency"));
    TH1 *CCDIS_muMtmEff2 = dynamic_cast<TH1 *>(histFile2->Get("CCDIS_MU_P_PIPLUS_MUON_MomentumEfficiency"));
    TH1 *CCDIS_pMtmEff2 = dynamic_cast<TH1 *>(histFile2->Get("CCDIS_MU_P_PIPLUS_PROTON1_MomentumEfficiency"));
    TH1 *CCDIS_piMtmEff2 = dynamic_cast<TH1 *>(histFile2->Get("CCDIS_MU_P_PIPLUS_PIPLUS_MomentumEfficiency"));
    if (CCDIS_muMtmEff1 && CCDIS_pMtmEff1 && CCDIS_piMtmEff1 && CCDIS_muMtmEff2 && CCDIS_pMtmEff2 && CCDIS_piMtmEff2)
    {

        SetHistColorStyle(CCDIS_muMtmEff1, kRed, kFullCircle);
        CCDIS_muMtmEff1->Draw();
        SetHistColorStyle(CCDIS_pMtmEff1, kRed, kOpenCircle);
        CCDIS_pMtmEff1->Draw("same");
        SetHistColorStyle(CCDIS_piMtmEff1, kOrange + 7, kFullDiamond);
        CCDIS_piMtmEff1->Draw("same");

        SetHistColorStyle(CCDIS_muMtmEff2, kBlue, kFullCircle);
        CCDIS_muMtmEff2->Draw("same");
        SetHistColorStyle(CCDIS_pMtmEff2, kBlue, kOpenCircle);
        CCDIS_pMtmEff2->Draw("same");
        SetHistColorStyle(CCDIS_piMtmEff2, kGreen + 3, kFullDiamond);
        CCDIS_piMtmEff2->Draw("same");

        CCDIS_Leg.Draw();
    }

    // Completeness
    TLegend CCDIS_Leg2(0.2, 0.625, 0.45, 0.875, "CCDIS (#mu p #pi)");
    CCDIS_Leg2.SetNColumns(3);
    CCDIS_Leg2.SetBorderSize(0);
    CCDIS_Leg2.SetTextSize(0.05);

    theCanvas2->cd(4);
    gPad->SetGrid();
    TH1 *CCDIS_muComp1 = dynamic_cast<TH1 *>(histFile1->Get("CCDIS_MU_P_PIPLUS_MUON_Completeness"));
    TH1 *CCDIS_pComp1 = dynamic_cast<TH1 *>(histFile1->Get("CCDIS_MU_P_PIPLUS_PROTON1_Completeness"));
    TH1 *CCDIS_piComp1 = dynamic_cast<TH1 *>(histFile1->Get("CCDIS_MU_P_PIPLUS_PIPLUS_Completeness"));
    TH1 *CCDIS_muComp2 = dynamic_cast<TH1 *>(histFile2->Get("CCDIS_MU_P_PIPLUS_MUON_Completeness"));
    TH1 *CCDIS_pComp2 = dynamic_cast<TH1 *>(histFile2->Get("CCDIS_MU_P_PIPLUS_PROTON1_Completeness"));
    TH1 *CCDIS_piComp2 = dynamic_cast<TH1 *>(histFile2->Get("CCDIS_MU_P_PIPLUS_PIPLUS_Completeness"));

    if (CCDIS_muComp1 && CCDIS_pComp1 && CCDIS_piComp1 && CCDIS_muComp2 && CCDIS_pComp2 && CCDIS_piComp2)
    {

        SetHistColorStyle(CCDIS_muComp1, kRed, kFullCircle);
        CCDIS_muComp1->SetTitleOffset(1.3, "Y");
        CCDIS_Leg2.AddEntry(CCDIS_muComp1, "#mu");
        CCDIS_muComp1->Draw();
        gPad->SetLogy();

        SetHistColorStyle(CCDIS_pComp1, kRed, kOpenCircle);
        CCDIS_Leg2.AddEntry(CCDIS_pComp1, "p");
        CCDIS_pComp1->Draw("same");

        SetHistColorStyle(CCDIS_piComp1, kOrange + 7, kFullDiamond);
        CCDIS_Leg2.AddEntry(CCDIS_piComp1, "#pi");
        CCDIS_piComp1->Draw("same");

        SetHistColorStyle(CCDIS_muComp2, kBlue, kFullCircle);
        CCDIS_muComp2->SetTitleOffset(1.3, "Y");
        CCDIS_Leg2.AddEntry(CCDIS_muComp2, "#mu");
        CCDIS_muComp2->Draw("same");

        SetHistColorStyle(CCDIS_pComp2, kBlue, kOpenCircle);
        CCDIS_Leg2.AddEntry(CCDIS_pComp2, "p");
        CCDIS_pComp2->Draw("same");

        SetHistColorStyle(CCDIS_piComp2, kGreen + 3, kFullDiamond);
        CCDIS_Leg2.AddEntry(CCDIS_piComp2, "#pi");
        CCDIS_piComp2->Draw("same");

        CCDIS_Leg2.Draw();
    }

    // Purity
    theCanvas2->cd(5);
    gPad->SetGrid();
    TH1 *CCDIS_muPure1 = dynamic_cast<TH1 *>(histFile1->Get("CCDIS_MU_P_PIPLUS_MUON_Purity"));
    TH1 *CCDIS_pPure1 = dynamic_cast<TH1 *>(histFile1->Get("CCDIS_MU_P_PIPLUS_PROTON1_Purity"));
    TH1 *CCDIS_piPure1 = dynamic_cast<TH1 *>(histFile1->Get("CCDIS_MU_P_PIPLUS_PIPLUS_Purity"));
    TH1 *CCDIS_muPure2 = dynamic_cast<TH1 *>(histFile2->Get("CCDIS_MU_P_PIPLUS_MUON_Purity"));
    TH1 *CCDIS_pPure2 = dynamic_cast<TH1 *>(histFile2->Get("CCDIS_MU_P_PIPLUS_PROTON1_Purity"));
    TH1 *CCDIS_piPure2 = dynamic_cast<TH1 *>(histFile2->Get("CCDIS_MU_P_PIPLUS_PIPLUS_Purity"));

    if (CCDIS_muPure1 && CCDIS_pPure1 && CCDIS_piPure1 && CCDIS_muPure2 && CCDIS_pPure2 && CCDIS_piPure2)
    {

        SetHistColorStyle(CCDIS_muPure1, kRed, kFullCircle);
        CCDIS_muPure1->SetTitleOffset(1.3, "Y");
        CCDIS_muPure1->Draw();
        gPad->SetLogy();

        SetHistColorStyle(CCDIS_pPure1, kRed, kOpenCircle);
        CCDIS_pPure1->Draw("same");

        SetHistColorStyle(CCDIS_piPure1, kOrange + 7, kFullDiamond);
        CCDIS_piPure1->Draw("same");

        SetHistColorStyle(CCDIS_muPure2, kBlue, kFullCircle);
        CCDIS_muPure2->Draw("same");

        SetHistColorStyle(CCDIS_pPure2, kBlue, kOpenCircle);
        CCDIS_pPure2->Draw("same");

        SetHistColorStyle(CCDIS_piPure2, kGreen + 3, kFullDiamond);
        CCDIS_piPure2->Draw("same");

        CCDIS_Leg2.Draw();
    }

    // Vtx dR
    theCanvas2->cd(6);
    gPad->SetGrid();
    TH1 *CCDIS_muVtxR1 = dynamic_cast<TH1 *>(histFile1->Get("CCDIS_MU_P_PIPLUS_VtxDeltaR"));
    TH1 *CCDIS_muVtxR2 = dynamic_cast<TH1 *>(histFile2->Get("CCDIS_MU_P_PIPLUS_VtxDeltaR"));
    if (CCDIS_muVtxR1 && CCDIS_muVtxR2)
    {
        NormHist(CCDIS_muVtxR1, "Fraction of events");
        NormHist(CCDIS_muVtxR2, "Fraction of events");
        //FitCBFun(CCDIS_muVtxR1, 0.4, 0.2, 0.2, kRed);
        //FitCBFun(CCDIS_muVtxR2, 0.3, 0.3, 0.3, kBlue);
        CCDIS_muVtxR1->SetTitleOffset(1.3, "Y");
        SetHistColorStyle(CCDIS_muVtxR1, kRed);
        SetHistColorStyle(CCDIS_muVtxR2, kBlue);

        if (CCDIS_muVtxR1->GetMaximum() > CCDIS_muVtxR2->GetMaximum())
        {
            CCDIS_muVtxR1->Draw();
            CCDIS_muVtxR2->Draw("same");
        }
        else
        {
            CCDIS_muVtxR2->Draw();
            CCDIS_muVtxR1->Draw("same");
        }

        // Print parameters
        latex->SetTextSize(0.045);
        latex->SetTextColor(kRed);
        latex->DrawLatexNDC(0.335, 0.75, PrintCBPars(CCDIS_muVtxR1, label1.c_str()));
        latex->SetTextColor(kBlue);
        latex->DrawLatexNDC(0.335, 0.65, PrintCBPars(CCDIS_muVtxR2, label2.c_str()));
    }

    theCanvas2->Print("CCDIS_pmupi.png");
    theCanvas2->Print("CCDIS_pmupi.pdf");

    // CCDIS: numu + Ar -> mu- p pi0
    theCanvas2->Clear();
    theCanvas2->Divide(3, 2);
    theCanvas2->cd(1);
    gPad->SetGrid();

    TLegend CCDISpi0_Leg(0.65, 0.4, 0.875, 0.7, "CCDIS (#mu p #pi^{0})");
    CCDISpi0_Leg.SetNColumns(2);
    CCDISpi0_Leg.SetBorderSize(0);
    CCDISpi0_Leg.SetTextSize(0.05);

    // Photon hits efficiency
    TH1 *CCDISpi0_g1HitsEff1 = dynamic_cast<TH1 *>(histFile1->Get("CCDIS_MU_P_PIZERO_PHOTON1_HitsEfficiency"));
    TH1 *CCDISpi0_g2HitsEff1 = dynamic_cast<TH1 *>(histFile1->Get("CCDIS_MU_P_PIZERO_PHOTON2_HitsEfficiency"));
    TH1 *CCDISpi0_g1HitsEff2 = dynamic_cast<TH1 *>(histFile2->Get("CCDIS_MU_P_PIZERO_PHOTON1_HitsEfficiency"));
    TH1 *CCDISpi0_g2HitsEff2 = dynamic_cast<TH1 *>(histFile2->Get("CCDIS_MU_P_PIZERO_PHOTON2_HitsEfficiency"));

    if (CCDISpi0_g1HitsEff1 && CCDISpi0_g2HitsEff1 && CCDISpi0_g1HitsEff2 && CCDISpi0_g2HitsEff2)
    {

        CCDISpi0_g1HitsEff1->GetXaxis()->SetRangeUser(0.0, xMax);
        gPad->SetLogx();
        SetHistColorStyle(CCDISpi0_g1HitsEff1, kRed, kFullCircle);
        CCDISpi0_Leg.AddEntry(CCDISpi0_g1HitsEff1, "#gamma_{1}");
        CCDISpi0_g1HitsEff1->Draw();

        CCDISpi0_g2HitsEff1->GetXaxis()->SetRangeUser(0.0, xMax);
        SetHistColorStyle(CCDISpi0_g2HitsEff1, kRed, kOpenCircle);
        CCDISpi0_Leg.AddEntry(CCDISpi0_g2HitsEff1, "#gamma_{2}");
        CCDISpi0_g2HitsEff1->Draw("same");

        CCDISpi0_g1HitsEff2->GetXaxis()->SetRangeUser(0.0, xMax);
        SetHistColorStyle(CCDISpi0_g1HitsEff2, kBlue, kFullCircle);
        CCDISpi0_Leg.AddEntry(CCDISpi0_g1HitsEff2, "#gamma_{1}");
        CCDISpi0_g1HitsEff2->Draw("same");

        CCDISpi0_g2HitsEff2->GetXaxis()->SetRangeUser(0.0, xMax);
        SetHistColorStyle(CCDISpi0_g2HitsEff2, kBlue, kOpenCircle);
        CCDISpi0_Leg.AddEntry(CCDISpi0_g2HitsEff2, "#gamma_{2}");
        CCDISpi0_g2HitsEff2->Draw("same");

        CCDISpi0_Leg.Draw();
    }

    // Momentum efficiency
    theCanvas2->cd(2);
    gPad->SetGrid();
    TH1 *CCDISpi0_g1MtmEff1 = dynamic_cast<TH1 *>(histFile1->Get("CCDIS_MU_P_PIZERO_PHOTON1_MomentumEfficiency"));
    TH1 *CCDISpi0_g2MtmEff1 = dynamic_cast<TH1 *>(histFile1->Get("CCDIS_MU_P_PIZERO_PHOTON2_MomentumEfficiency"));
    TH1 *CCDISpi0_g1MtmEff2 = dynamic_cast<TH1 *>(histFile2->Get("CCDIS_MU_P_PIZERO_PHOTON1_MomentumEfficiency"));
    TH1 *CCDISpi0_g2MtmEff2 = dynamic_cast<TH1 *>(histFile2->Get("CCDIS_MU_P_PIZERO_PHOTON2_MomentumEfficiency"));

    if (CCDISpi0_g1MtmEff1 && CCDISpi0_g2MtmEff1 && CCDISpi0_g1MtmEff2 && CCDISpi0_g2MtmEff2)
    {

        SetHistColorStyle(CCDISpi0_g1MtmEff1, kRed, kFullCircle);
        CCDISpi0_g1MtmEff1->Draw();
        SetHistColorStyle(CCDISpi0_g2MtmEff1, kRed, kOpenCircle);
        CCDISpi0_g2MtmEff1->Draw("same");

        SetHistColorStyle(CCDISpi0_g1MtmEff2, kBlue, kFullCircle);
        CCDISpi0_g1MtmEff2->Draw("same");
        SetHistColorStyle(CCDISpi0_g2MtmEff2, kBlue, kOpenCircle);
        CCDISpi0_g2MtmEff2->Draw("same");

        CCDISpi0_Leg.Draw();
    }

    // Completeness
    TLegend CCDISpi0_Leg2(0.2, 0.65, 0.4, 0.85, "CCDIS (#mu p #pi^{0})");
    CCDISpi0_Leg2.SetNColumns(2);
    CCDISpi0_Leg2.SetBorderSize(0);
    CCDISpi0_Leg2.SetTextSize(0.05);

    theCanvas2->cd(4);
    TH1 *CCDISpi0_g1Comp1 = dynamic_cast<TH1 *>(histFile1->Get("CCDIS_MU_P_PIZERO_PHOTON1_Completeness"));
    TH1 *CCDISpi0_g2Comp1 = dynamic_cast<TH1 *>(histFile1->Get("CCDIS_MU_P_PIZERO_PHOTON2_Completeness"));
    TH1 *CCDISpi0_g1Comp2 = dynamic_cast<TH1 *>(histFile2->Get("CCDIS_MU_P_PIZERO_PHOTON1_Completeness"));
    TH1 *CCDISpi0_g2Comp2 = dynamic_cast<TH1 *>(histFile2->Get("CCDIS_MU_P_PIZERO_PHOTON2_Completeness"));

    if (CCDISpi0_g1Comp1 && CCDISpi0_g2Comp1 && CCDISpi0_g1Comp2 && CCDISpi0_g2Comp2)
    {

        CCDISpi0_Leg2.AddEntry(CCDISpi0_g1Comp1, "#gamma_{1}");
        SetHistColorStyle(CCDISpi0_g1Comp1, kRed, kFullCircle);
        CCDISpi0_g1Comp1->SetTitleOffset(1.3, "Y");
        CCDISpi0_g1Comp1->Draw();
        gPad->SetLogy();

        SetHistColorStyle(CCDISpi0_g2Comp1, kRed, kOpenCircle);
        CCDISpi0_Leg2.AddEntry(CCDISpi0_g2Comp1, "#gamma_{2}");
        CCDISpi0_g2Comp1->Draw("same");

        SetHistColorStyle(CCDISpi0_g1Comp2, kBlue, kFullCircle);
        CCDISpi0_Leg2.AddEntry(CCDISpi0_g1Comp2, "#gamma_{1}");
        CCDISpi0_g1Comp2->Draw("same");

        SetHistColorStyle(CCDISpi0_g2Comp2, kBlue, kOpenCircle);
        CCDISpi0_Leg2.AddEntry(CCDISpi0_g2Comp2, "#gamma_{2}");
        CCDISpi0_g2Comp2->Draw("same");

        CCDISpi0_Leg2.Draw();
    }

    // Purity
    theCanvas2->cd(5);
    gPad->SetGrid();
    TH1 *CCDISpi0_g1Pure1 = dynamic_cast<TH1 *>(histFile1->Get("CCDIS_MU_P_PIZERO_PHOTON1_Purity"));
    TH1 *CCDISpi0_g2Pure1 = dynamic_cast<TH1 *>(histFile1->Get("CCDIS_MU_P_PIZERO_PHOTON2_Purity"));
    TH1 *CCDISpi0_g1Pure2 = dynamic_cast<TH1 *>(histFile2->Get("CCDIS_MU_P_PIZERO_PHOTON1_Purity"));
    TH1 *CCDISpi0_g2Pure2 = dynamic_cast<TH1 *>(histFile2->Get("CCDIS_MU_P_PIZERO_PHOTON2_Purity"));

    if (CCDISpi0_g1Pure1 && CCDISpi0_g2Pure1 && CCDISpi0_g1Pure1 && CCDISpi0_g2Pure1)
    {

        SetHistColorStyle(CCDISpi0_g1Pure1, kRed, kFullCircle);
        CCDISpi0_g1Pure1->SetTitleOffset(1.3, "Y");
        CCDISpi0_g1Pure1->Draw();
        gPad->SetLogy();

        SetHistColorStyle(CCDISpi0_g2Pure1, kRed, kOpenCircle);
        CCDISpi0_g2Pure1->Draw("same");

        SetHistColorStyle(CCDISpi0_g1Pure2, kBlue, kFullCircle);
        CCDISpi0_g1Pure2->Draw("same");

        SetHistColorStyle(CCDISpi0_g2Pure2, kBlue, kOpenCircle);
        CCDISpi0_g2Pure2->Draw("same");

        CCDISpi0_Leg2.Draw();
    }

    // Vtx dR
    theCanvas2->cd(6);
    gPad->SetGrid();
    TH1 *CCDISpi0_muVtxR1 = dynamic_cast<TH1 *>(histFile1->Get("CCDIS_MU_P_PIZERO_VtxDeltaR"));
    TH1 *CCDISpi0_muVtxR2 = dynamic_cast<TH1 *>(histFile2->Get("CCDIS_MU_P_PIZERO_VtxDeltaR"));

    if (CCDISpi0_muVtxR1 && CCDISpi0_muVtxR2)
    {
        NormHist(CCDISpi0_muVtxR1, "Fraction of events");
        NormHist(CCDISpi0_muVtxR2, "Fraction of events");
        //FitCBFun(CCDISpi0_muVtxR1, 0.3, 0.09, 0.2, kRed);
        //FitCBFun(CCDISpi0_muVtxR2, 0.3, 0.3, 0.3, kBlue);
        CCDISpi0_muVtxR1->SetTitleOffset(1.3, "Y");
        SetHistColorStyle(CCDISpi0_muVtxR1, kRed);
        SetHistColorStyle(CCDISpi0_muVtxR2, kBlue);

        if (CCDISpi0_muVtxR1->GetMaximum() > CCDISpi0_muVtxR2->GetMaximum())
        {
            CCDISpi0_muVtxR1->Draw();
            CCDISpi0_muVtxR2->Draw("same");
        }
        else
        {
            CCDISpi0_muVtxR2->Draw();
            CCDISpi0_muVtxR1->Draw("same");
        }

        // Print parameters
        latex->SetTextSize(0.045);
        latex->SetTextColor(kRed);
        latex->DrawLatexNDC(0.335, 0.75, label1.c_str());
        latex->DrawLatexNDC(0.335, 0.75, PrintCBPars(CCDISpi0_muVtxR1, label1.c_str()));
        latex->SetTextColor(kBlue);
        latex->DrawLatexNDC(0.335, 0.65, label2.c_str());
        latex->DrawLatexNDC(0.335, 0.65, PrintCBPars(CCDISpi0_muVtxR2, label2.c_str()));
    }

    theCanvas2->Print("CCDIS_pmupi0.png");
    theCanvas2->Print("CCDIS_pmupi0.pdf");

    plotCorrectFrac(validFileName1, label1, validFileName2, label2);
}
