// First create histogram ROOT file then run this script.
// First the validation script needs to be loaded:
//   .L Validation.C
//   .x createPlots.C("output/Validation.root", "output/Histos.root", true)

#include <cmath>
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
    if (!alphaR && alphaL < 0.0)
    {
        std::swap(alphaL, alphaR);
        std::swap(nL, nR);
    }

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

void FitCBFun(TH1 *hist, double theMean = -999.0, double theSigma = -999.0)
{
    TString histName = hist->GetName();
    TString funName = histName + "Fun";
    TAxis *xAxis = hist->GetXaxis();
    const double xMin = xAxis->GetXmin();
    const double xMax = xAxis->GetXmax();
    const double mean = theMean > -999.0 ? theMean : hist->GetMean();
    const double sigma = theSigma > -999.0 ? theSigma : hist->GetStdDev();
    const double norm = hist->GetMaximum();

    TF1 *fun = new TF1(funName.Data(), &RooCrystalBallFun, xMin, xMax, 8);
    fun->SetParameters(norm, mean, sigma, sigma, 1.0, 1.0, 1.0, 1.0);
    fun->SetParNames("N", "#mu", "#sigma_{L}", "#sigma_{R}", "#alpha_{L}", "n_{L}", "#alpha_{R}", "n_{R}");
    hist->Fit(funName.Data());
}

TGraphErrors *getEfficiency(TTree *theTree, int interactionCode, int minProtons, int maxProtons)
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

void plotCorrectFrac(const std::string &validFileName)
{

    std::cout << "Plotting correct event fraction vs number of protons" << std::endl;
    TCanvas *theCanvasC = new TCanvas("theCanvasC", "", 1400, 700);
    theCanvasC->UseCurrentStyle();
    theCanvasC->Clear();
    theCanvasC->Divide(2, 2);

    TFile *validFile = TFile::Open(validFileName.c_str(), "read");
    TTree *theTree = dynamic_cast<TTree *>(validFile->Get("Validation"));

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
    latex->SetTextSize(0.05);

    // CC
    TH2F *nullHist = new TH2F("nullHist", "", 2, -1.0, 6.0, 2, 0.0, 100.0);
    nullHist->SetXTitle("Number of protons");
    nullHist->SetYTitle("Correct events (%)");
    nullHist->SetNdivisions(7, "X");

    TGraphErrors *eff_CCQEL_MU = getEfficiency(theTree, CCQEL_MU, 0, 5);
    theCanvasC->cd(1);
    nullHist->Draw();
    eff_CCQEL_MU->Draw("same");
    latex->DrawLatex(3.5, 80.0, "CCQE #mu Np");

    TGraphErrors *eff_CCDIS_MU = getEfficiency(theTree, CCDIS_MU, 0, 5);
    theCanvasC->cd(2);
    nullHist->Draw();
    eff_CCDIS_MU->Draw("same");
    latex->DrawLatex(3.5, 80.0, "CCDIS #mu Np");

    TGraphErrors *eff_CCDIS_MU_PIPLUS = getEfficiency(theTree, CCDIS_MU_PIPLUS, 0, 5);
    theCanvasC->cd(3);
    nullHist->Draw();
    eff_CCDIS_MU_PIPLUS->Draw("same");
    latex->DrawLatex(3.5, 80.0, "CCDIS #mu #pi^{#pm} Np");

    TGraphErrors *eff_CCDIS_MU_PIZERO = getEfficiency(theTree, CCDIS_MU_PIZERO, 0, 5);
    theCanvasC->cd(4);
    nullHist->Draw();
    eff_CCDIS_MU_PIZERO->Draw("same");
    latex->DrawLatex(3.5, 80.0, "CCDIS #mu #pi^{0} Np");

    theCanvasC->Print("CCEventFrac.png");

    // NC
    TGraphErrors *eff_NCDIS_P = getEfficiency(theTree, NCDIS_P, 1, 5);
    theCanvasC->cd(1);
    nullHist->Draw();
    eff_NCDIS_P->Draw("same");
    latex->DrawLatex(3.5, 80.0, "NCDIS Np");

    TGraphErrors *eff_NCDIS_PIPLUS = getEfficiency(theTree, NCDIS_PIPLUS, 0, 5);
    theCanvasC->cd(2);
    nullHist->Draw();
    eff_NCDIS_PIPLUS->Draw("same");
    latex->DrawLatex(3.5, 80.0, "NCDIS #pi^{+} Np");

    TGraphErrors *eff_NCDIS_PIMINUS = getEfficiency(theTree, NCDIS_PIMINUS, 0, 5);
    theCanvasC->cd(3);
    nullHist->Draw();
    eff_NCDIS_PIMINUS->Draw("same");
    latex->DrawLatex(3.5, 80.0, "NCDIS #pi^{-} Np");

    TGraphErrors *eff_NCDIS_PIZERO = getEfficiency(theTree, NCDIS_PIZERO, 0, 5);
    theCanvasC->cd(4);
    nullHist->Draw();
    eff_NCDIS_PIZERO->Draw("same");
    latex->DrawLatex(3.5, 80.0, "NCDIS #pi^{0} Np");

    theCanvasC->Print("NCEventFrac.png");

    validFile->Close();
}

void createPlots(const std::string &validFileName = "Validation.root", const std::string &histFileName = "Histos.root", bool runValidation = false)
{

    gROOT->SetStyle("Plain");
    gStyle->SetOptStat(0);
    TCanvas *theCanvas = new TCanvas("theCanvas", "", 900, 700);
    theCanvas->UseCurrentStyle();

    // Run validation code to create histograms
    Parameters p;
    p.m_histogramOutput = true;
    p.m_histFileName = histFileName;
    // Set x vertex offset to 0
    p.m_vertexXCorrection = 0.0;
    // Specify map file
    p.m_mapFileName = "Map.txt";
    // Reduce verbosity
    p.m_displayMatchedEvents = false;
    if (runValidation)
    {
        Validation(validFileName, p);
    }

    // Plot histograms
    theCanvas->Divide(2, 2);

    TLatex *text = new TLatex();
    text->SetTextSize(0.075);
    text->SetNDC(true);

    TFile *histFile = TFile::Open(histFileName.c_str(), "read");

    // ALL HITS EFFICIENCY
    // maxMerged voxel limit
    double xMax = 10000;

    theCanvas->cd(1);
    TH1 *muHitsEff = dynamic_cast<TH1 *>(histFile->Get("ALL_INTERACTIONS_MUON_HitsEfficiency"));
    if (muHitsEff)
    {
        muHitsEff->GetXaxis()->SetRangeUser(0.0, xMax);
        gPad->SetLogx();
        muHitsEff->Draw();
        text->DrawLatex(0.775, 0.45, "#mu");
    }

    theCanvas->cd(2);
    TH1 *pHitsEff = dynamic_cast<TH1 *>(histFile->Get("ALL_INTERACTIONS_PROTON1_HitsEfficiency"));
    if (pHitsEff)
    {
        pHitsEff->GetXaxis()->SetRangeUser(0.0, xMax);
        gPad->SetLogx();
        pHitsEff->Draw();
        text->DrawLatex(0.775, 0.45, "p");
    }

    theCanvas->cd(3);
    TH1 *pipHitsEff = dynamic_cast<TH1 *>(histFile->Get("ALL_INTERACTIONS_PIPLUS_HitsEfficiency"));
    if (pipHitsEff)
    {
        pipHitsEff->GetXaxis()->SetRangeUser(0.0, xMax);
        gPad->SetLogx();
        pipHitsEff->Draw();
        text->DrawLatex(0.775, 0.45, "#pi^{+}");
    }

    theCanvas->cd(4);
    TH1 *pimHitsEff = dynamic_cast<TH1 *>(histFile->Get("ALL_INTERACTIONS_PIMINUS_HitsEfficiency"));
    if (pimHitsEff)
    {
        pimHitsEff->GetXaxis()->SetRangeUser(0.0, xMax);
        gPad->SetLogx();
        pimHitsEff->Draw();
        text->DrawLatex(0.775, 0.45, "#pi^{-}");
    }
    theCanvas->Print("allHitsEff.png");

    // ALL MOMENTUM EFFICIENCY
    theCanvas->cd(1);
    gPad->SetLogx(0);
    TH1 *muMtmEff = dynamic_cast<TH1 *>(histFile->Get("ALL_INTERACTIONS_MUON_MomentumEfficiency"));
    if (muMtmEff)
    {
        muMtmEff->Draw();
        text->DrawLatex(0.775, 0.45, "#mu");
    }

    theCanvas->cd(2);
    gPad->SetLogx(0);
    TH1 *pMtmEff = dynamic_cast<TH1 *>(histFile->Get("ALL_INTERACTIONS_PROTON1_MomentumEfficiency"));
    if (pMtmEff)
    {
        pMtmEff->Draw();
        text->DrawLatex(0.775, 0.45, "p");
    }

    theCanvas->cd(3);
    gPad->SetLogx(0);
    TH1 *pipMtmEff = dynamic_cast<TH1 *>(histFile->Get("ALL_INTERACTIONS_PIPLUS_MomentumEfficiency"));
    if (pipMtmEff)
    {
        pipMtmEff->Draw();
        text->DrawLatex(0.775, 0.45, "#pi^{+}");
    }

    theCanvas->cd(4);
    gPad->SetLogx(0);
    TH1 *pimMtmEff = dynamic_cast<TH1 *>(histFile->Get("ALL_INTERACTIONS_PIMINUS_MomentumEfficiency"));
    if (pimMtmEff)
    {
        pimMtmEff->Draw();
        text->DrawLatex(0.775, 0.45, "#pi^{-}");
    }

    theCanvas->Print("allMtmEff.png");

    // ALL VERTICES
    gStyle->SetOptFit(1111);
    gROOT->ForceStyle();
    theCanvas->UseCurrentStyle();
    theCanvas->Update();

    theCanvas->cd(1);
    TH1 *xVtxHist = dynamic_cast<TH1 *>(histFile->Get("ALL_INTERACTIONS_VtxDeltaX"));
    xVtxHist->GetYaxis()->SetTitleOffset(1.35);
    FitCBFun(xVtxHist);
    xVtxHist->Draw();
    gPad->Update();

    theCanvas->cd(2);
    TH1 *yVtxHist = dynamic_cast<TH1 *>(histFile->Get("ALL_INTERACTIONS_VtxDeltaY"));
    yVtxHist->GetYaxis()->SetTitleOffset(1.35);
    FitCBFun(yVtxHist);
    yVtxHist->Draw();
    gPad->Update();

    theCanvas->cd(3);
    TH1 *zVtxHist = dynamic_cast<TH1 *>(histFile->Get("ALL_INTERACTIONS_VtxDeltaZ"));
    zVtxHist->GetYaxis()->SetTitleOffset(1.35);
    FitCBFun(zVtxHist);
    zVtxHist->Draw();
    gPad->Update();

    theCanvas->cd(4);
    TH1 *rVtxHist = dynamic_cast<TH1 *>(histFile->Get("ALL_INTERACTIONS_VtxDeltaR"));
    rVtxHist->GetYaxis()->SetTitleOffset(1.35);
    FitCBFun(rVtxHist, 0.3, 0.3);
    rVtxHist->Draw();
    gPad->Update();

    theCanvas->Print("allVtx.png");

    // 2D VTX RESIDUALS
    theCanvas->Clear();
    theCanvas->Divide(2, 2);
    theCanvas->UseCurrentStyle();
    theCanvas->Update();

    theCanvas->cd(1);
    TH1 *xyVtxHist = dynamic_cast<TH1 *>(histFile->Get("ALL_INTERACTIONS_VtxDeltaXY"));
    xyVtxHist->GetYaxis()->SetTitleOffset(1.35);
    xyVtxHist->Draw("colz");
    gPad->SetLogz();
    gPad->Update();

    theCanvas->cd(2);
    TH1 *zyVtxHist = dynamic_cast<TH1 *>(histFile->Get("ALL_INTERACTIONS_VtxDeltaZY"));
    zyVtxHist->GetYaxis()->SetTitleOffset(1.35);
    zyVtxHist->Draw("colz");
    gPad->SetLogz();
    gPad->Update();

    theCanvas->cd(3);
    TH1 *zxVtxHist = dynamic_cast<TH1 *>(histFile->Get("ALL_INTERACTIONS_VtxDeltaZX"));
    zxVtxHist->GetYaxis()->SetTitleOffset(1.35);
    zxVtxHist->Draw("colz");
    gPad->SetLogz();
    gPad->Update();

    theCanvas->Print("vtx2D.png");

    delete theCanvas;

    // CCQE: numu + Ar -> mu + p
    TCanvas *theCanvas2 = new TCanvas("theCanvas2", "", 1400, 700);
    theCanvas2->UseCurrentStyle();
    theCanvas2->Clear();
    theCanvas2->Divide(3, 2);

    theCanvas2->cd(1);

    TLegend CCQE_Leg(0.6, 0.2, 0.8, 0.4, "CCQE (#mu p)");
    CCQE_Leg.SetBorderSize(0);
    CCQE_Leg.SetTextSize(0.05);

    // Hits efficiency
    TH1 *CCQE_muHitsEff = dynamic_cast<TH1 *>(histFile->Get("CCQEL_MU_P_MUON_HitsEfficiency"));
    TH1 *CCQE_pHitsEff = dynamic_cast<TH1 *>(histFile->Get("CCQEL_MU_P_PROTON1_HitsEfficiency"));
    if (CCQE_muHitsEff && CCQE_pHitsEff)
    {
        CCQE_muHitsEff->GetXaxis()->SetRangeUser(0.0, xMax);
        gPad->SetLogx();
        CCQE_muHitsEff->SetLineColor(kRed);
        CCQE_Leg.AddEntry(CCQE_muHitsEff, "#mu");
        CCQE_muHitsEff->Draw();

        CCQE_pHitsEff->GetXaxis()->SetRangeUser(0.0, xMax);
        CCQE_pHitsEff->SetLineColor(kBlue);
        CCQE_Leg.AddEntry(CCQE_pHitsEff, "p");
        CCQE_pHitsEff->Draw("same");
        CCQE_Leg.Draw();
    }

    // Momentum efficiency
    theCanvas2->cd(2);
    TH1 *CCQE_muMtmEff = dynamic_cast<TH1 *>(histFile->Get("CCQEL_MU_P_MUON_MomentumEfficiency"));
    TH1 *CCQE_pMtmEff = dynamic_cast<TH1 *>(histFile->Get("CCQEL_MU_P_PROTON1_MomentumEfficiency"));
    if (CCQE_muMtmEff && CCQE_pMtmEff)
    {
        CCQE_muMtmEff->SetLineColor(kRed);
        CCQE_muMtmEff->Draw();

        CCQE_pMtmEff->SetLineColor(kBlue);
        CCQE_pMtmEff->Draw("same");
        CCQE_Leg.Draw();
    }

    // Completeness
    TLegend CCQE_Leg2(0.2, 0.6, 0.4, 0.8, "CCQE (#mu p)");
    CCQE_Leg2.SetBorderSize(0);
    CCQE_Leg2.SetTextSize(0.05);

    theCanvas2->cd(4);
    TH1 *CCQE_muComp = dynamic_cast<TH1 *>(histFile->Get("CCQEL_MU_P_MUON_Completeness"));
    TH1 *CCQE_pComp = dynamic_cast<TH1 *>(histFile->Get("CCQEL_MU_P_PROTON1_Completeness"));
    if (CCQE_muComp && CCQE_pComp)
    {
        CCQE_muComp->SetLineColor(kRed);
        CCQE_muComp->SetTitleOffset(1.3, "Y");
        CCQE_Leg2.AddEntry(CCQE_muComp, "#mu");
        CCQE_muComp->Draw();
        gPad->SetLogy();

        CCQE_pComp->SetLineColor(kBlue);
        CCQE_Leg2.AddEntry(CCQE_pComp, "p");
        CCQE_pComp->Draw("same");
        CCQE_Leg2.Draw();
    }

    // Purity
    theCanvas2->cd(5);
    TH1 *CCQE_muPure = dynamic_cast<TH1 *>(histFile->Get("CCQEL_MU_P_MUON_Purity"));
    TH1 *CCQE_pPure = dynamic_cast<TH1 *>(histFile->Get("CCQEL_MU_P_PROTON1_Purity"));
    if (CCQE_muPure && CCQE_pPure)
    {
        CCQE_muPure->SetLineColor(kRed);
        CCQE_muPure->SetTitleOffset(1.3, "Y");
        CCQE_muPure->Draw();
        gPad->SetLogy();

        CCQE_pPure->SetLineColor(kBlue);
        CCQE_pPure->Draw("same");
        CCQE_Leg2.Draw();
    }

    // Vtx dR
    theCanvas2->cd(6);
    TH1 *CCQE_muVtxR = dynamic_cast<TH1 *>(histFile->Get("CCQEL_MU_P_VtxDeltaR"));
    if (CCQE_muVtxR)
    {
        FitCBFun(CCQE_muVtxR, 0.4, 0.15);
        CCQE_muVtxR->SetLineColor(kBlack);
        CCQE_muVtxR->Draw();
    }
    theCanvas2->Print("CCQE.png");

    // CCDIS: numu + Ar -> mu- p pi+
    theCanvas2->Clear();
    theCanvas2->Divide(3, 2);
    theCanvas2->cd(1);

    TLegend CCDIS_Leg(0.6, 0.2, 0.8, 0.5, "CCDIS (#mu p #pi)");
    CCDIS_Leg.SetBorderSize(0);
    CCDIS_Leg.SetTextSize(0.05);

    // Hits efficiency
    TH1 *CCDIS_muHitsEff = dynamic_cast<TH1 *>(histFile->Get("CCDIS_MU_P_PIPLUS_MUON_HitsEfficiency"));
    TH1 *CCDIS_pHitsEff = dynamic_cast<TH1 *>(histFile->Get("CCDIS_MU_P_PIPLUS_PROTON1_HitsEfficiency"));
    TH1 *CCDIS_piHitsEff = dynamic_cast<TH1 *>(histFile->Get("CCDIS_MU_P_PIPLUS_PIPLUS_HitsEfficiency"));
    if (CCDIS_muHitsEff && CCDIS_pHitsEff && CCDIS_piHitsEff)
    {
        CCDIS_muHitsEff->GetXaxis()->SetRangeUser(0.0, xMax);
        gPad->SetLogx();
        CCDIS_muHitsEff->SetLineColor(kRed);
        CCDIS_Leg.AddEntry(CCDIS_muHitsEff, "#mu");
        CCDIS_muHitsEff->Draw();

        CCDIS_pHitsEff->GetXaxis()->SetRangeUser(0.0, xMax);
        CCDIS_pHitsEff->SetLineColor(kBlue);
        CCDIS_Leg.AddEntry(CCDIS_pHitsEff, "p");
        CCDIS_pHitsEff->Draw("same");

        CCDIS_piHitsEff->GetXaxis()->SetRangeUser(0.0, xMax);
        CCDIS_piHitsEff->SetLineColor(kGreen + 3);
        CCDIS_Leg.AddEntry(CCDIS_piHitsEff, "#pi");
        CCDIS_piHitsEff->Draw("same");
        CCDIS_Leg.Draw();
    }

    // Momentum efficiency
    theCanvas2->cd(2);
    TH1 *CCDIS_muMtmEff = dynamic_cast<TH1 *>(histFile->Get("CCDIS_MU_P_PIPLUS_MUON_MomentumEfficiency"));
    TH1 *CCDIS_pMtmEff = dynamic_cast<TH1 *>(histFile->Get("CCDIS_MU_P_PIPLUS_PROTON1_MomentumEfficiency"));
    TH1 *CCDIS_piMtmEff = dynamic_cast<TH1 *>(histFile->Get("CCDIS_MU_P_PIPLUS_PIPLUS_MomentumEfficiency"));
    if (CCDIS_muMtmEff && CCDIS_pMtmEff && CCDIS_piMtmEff)
    {
        CCDIS_muMtmEff->SetLineColor(kRed);
        CCDIS_muMtmEff->Draw();
        CCDIS_pMtmEff->SetLineColor(kBlue);
        CCDIS_pMtmEff->Draw("same");
        CCDIS_piMtmEff->SetLineColor(kGreen + 3);
        CCDIS_piMtmEff->Draw("same");
        CCDIS_Leg.Draw();
    }

    // Completeness
    TLegend CCDIS_Leg2(0.2, 0.6, 0.4, 0.8, "CCDIS (#mu p #pi)");
    CCDIS_Leg2.SetBorderSize(0);
    CCDIS_Leg2.SetTextSize(0.05);

    theCanvas2->cd(4);
    TH1 *CCDIS_muComp = dynamic_cast<TH1 *>(histFile->Get("CCDIS_MU_P_PIPLUS_MUON_Completeness"));
    TH1 *CCDIS_pComp = dynamic_cast<TH1 *>(histFile->Get("CCDIS_MU_P_PIPLUS_PROTON1_Completeness"));
    TH1 *CCDIS_piComp = dynamic_cast<TH1 *>(histFile->Get("CCDIS_MU_P_PIPLUS_PIPLUS_Completeness"));

    if (CCDIS_muComp && CCDIS_pComp && CCDIS_piComp)
    {
        CCDIS_muComp->SetLineColor(kRed);
        CCDIS_muComp->SetTitleOffset(1.3, "Y");
        CCDIS_Leg2.AddEntry(CCDIS_muComp, "#mu");
        CCDIS_muComp->Draw();
        gPad->SetLogy();
        CCDIS_pComp->SetLineColor(kBlue);
        CCDIS_Leg2.AddEntry(CCDIS_pComp, "p");
        CCDIS_pComp->Draw("same");
        CCDIS_piComp->SetLineColor(kGreen + 3);
        CCDIS_Leg2.AddEntry(CCDIS_piComp, "#pi");
        CCDIS_piComp->Draw("same");
        CCDIS_Leg2.Draw();
    }

    // Purity
    theCanvas2->cd(5);
    TH1 *CCDIS_muPure = dynamic_cast<TH1 *>(histFile->Get("CCDIS_MU_P_PIPLUS_MUON_Purity"));
    TH1 *CCDIS_pPure = dynamic_cast<TH1 *>(histFile->Get("CCDIS_MU_P_PIPLUS_PROTON1_Purity"));
    TH1 *CCDIS_piPure = dynamic_cast<TH1 *>(histFile->Get("CCDIS_MU_P_PIPLUS_PIPLUS_Purity"));

    if (CCDIS_muPure && CCDIS_pPure && CCDIS_piPure)
    {
        CCDIS_muPure->SetLineColor(kRed);
        CCDIS_muPure->SetTitleOffset(1.3, "Y");
        CCDIS_muPure->Draw();
        gPad->SetLogy();
        CCDIS_pPure->SetLineColor(kBlue);
        CCDIS_pPure->Draw("same");
        CCDIS_piPure->SetLineColor(kGreen + 3);
        CCDIS_piPure->Draw("same");
        CCDIS_Leg2.Draw();
    }

    // Vtx dR
    theCanvas2->cd(6);
    TH1 *CCDIS_muVtxR = dynamic_cast<TH1 *>(histFile->Get("CCDIS_MU_P_PIPLUS_VtxDeltaR"));
    if (CCDIS_muVtxR)
    {
        //FitCBFun(CCDIS_muVtxR, 0.4, 0.15);
        CCDIS_muVtxR->SetLineColor(kBlack);
        CCDIS_muVtxR->Draw();
    }

    theCanvas2->Print("CCDIS_pmupi.png");

    // CCDIS: numu + Ar -> mu- p pi0
    theCanvas2->Clear();
    theCanvas2->Divide(3, 2);
    theCanvas2->cd(1);

    TLegend CCDISpi0_Leg(0.65, 0.4, 0.875, 0.7, "CCDIS (#mu p #pi^{0})");
    CCDISpi0_Leg.SetBorderSize(0);
    CCDISpi0_Leg.SetTextSize(0.05);

    // Hits efficiency
    /*TH1* CCDISpi0_muHitsEff = dynamic_cast<TH1*>(histFile->Get("CCDIS_MU_P_PIZERO_MUON_HitsEfficiency"));
    CCDISpi0_muHitsEff->GetXaxis()->SetRangeUser(0.0, xMax);
    gPad->SetLogx();
    CCDISpi0_muHitsEff->SetLineColor(kRe}d);
    CCDISpi0_Leg.AddEntry(CCDISpi0_muHitsEff, "#mu");
    CCDISpi0_muHitsEff->Draw();
    TH1* CCDISpi0_pHitsEff = dynamic_cast<TH1*>(histFile->Get("CCDIS_MU_P_PIZERO_PROTON1_HitsEfficiency"));
    CCDISpi0_pHitsEff->GetXaxis()->SetRangeUser(0.0, xMax);
    CCDISpi0_pHitsEff->SetLineColor(kBlue);
    CCDISpi0_Leg.AddEntry(CCDISpi0_pHitsEff, "p");
    CCDISpi0_pHitsEff->Draw("same");*/
    TH1 *CCDISpi0_g1HitsEff = dynamic_cast<TH1 *>(histFile->Get("CCDIS_MU_P_PIZERO_PHOTON1_HitsEfficiency"));
    TH1 *CCDISpi0_g2HitsEff = dynamic_cast<TH1 *>(histFile->Get("CCDIS_MU_P_PIZERO_PHOTON2_HitsEfficiency"));

    if (CCDISpi0_g1HitsEff && CCDISpi0_g2HitsEff)
    {
        CCDISpi0_g1HitsEff->GetXaxis()->SetRangeUser(0.0, xMax);
        gPad->SetLogx();
        CCDISpi0_g1HitsEff->SetLineColor(kBlue);
        CCDISpi0_Leg.AddEntry(CCDISpi0_g1HitsEff, "#gamma_{1}");
        CCDISpi0_g1HitsEff->Draw();
        CCDISpi0_g2HitsEff->GetXaxis()->SetRangeUser(0.0, xMax);
        CCDISpi0_g2HitsEff->SetLineColor(kOrange + 7);
        CCDISpi0_Leg.AddEntry(CCDISpi0_g2HitsEff, "#gamma_{2}");
        CCDISpi0_g2HitsEff->Draw("same");
        CCDISpi0_Leg.Draw();
    }

    // Momentum efficiency
    theCanvas2->cd(2);
    /*TH1* CCDISpi0_muMtmEff = dynamic_cast<TH1*>(histFile->Get("CCDIS_MU_P_PIZERO_MUON_MomentumEfficiency"));
    CCDISpi0_muMtmEff->SetLineColor(kRed);
    CCDISpi0_muMtmEff->Draw();
    TH1* CCDISpi0_pMtmEff = dynamic_cast<TH1*>(histFile->Get("CCDIS_MU_P_PIZERO_PROTON1_MomentumEfficiency"));
    CCDISpi0_pMtmEff->SetLineColor(kBlue);
    CCDISpi0_pMtmEff->Draw("same");*/
    TH1 *CCDISpi0_g1MtmEff = dynamic_cast<TH1 *>(histFile->Get("CCDIS_MU_P_PIZERO_PHOTON1_MomentumEfficiency"));
    TH1 *CCDISpi0_g2MtmEff = dynamic_cast<TH1 *>(histFile->Get("CCDIS_MU_P_PIZERO_PHOTON2_MomentumEfficiency"));

    if (CCDISpi0_g1MtmEff && CCDISpi0_g2MtmEff)
    {
        CCDISpi0_g1MtmEff->SetLineColor(kBlue);
        CCDISpi0_g1MtmEff->Draw();
        CCDISpi0_g2MtmEff->SetLineColor(kOrange + 7);
        CCDISpi0_g2MtmEff->Draw("same");
        CCDISpi0_Leg.Draw();
    }

    // Completeness
    TLegend CCDISpi0_Leg2(0.2, 0.6, 0.4, 0.8, "CCDIS (#mu p #pi^{0})");
    CCDISpi0_Leg2.SetBorderSize(0);
    CCDISpi0_Leg2.SetTextSize(0.05);

    theCanvas2->cd(4);
    /*TH1* CCDISpi0_muComp = dynamic_cast<TH1*>(histFile->Get("CCDIS_MU_P_PIZERO_MUON_Completeness"));
    CCDISpi0_muComp->SetLineColor(kRed);
    CCDISpi0_Leg2.AddEntry(CCDISpi0_muComp, "#mu");
    CCDISpi0_muComp->Draw();
    TH1* CCDISpi0_pComp = dynamic_cast<TH1*>(histFile->Get("CCDIS_MU_P_PIZERO_PROTON1_Completeness"));
    CCDISpi0_pComp->SetLineColor(kBlue);
    CCDISpi0_Leg2.AddEntry(CCDISpi0_pComp, "p");
    CCDISpi0_pComp->Draw("same");*/
    TH1 *CCDISpi0_g1Comp = dynamic_cast<TH1 *>(histFile->Get("CCDIS_MU_P_PIZERO_PHOTON1_Completeness"));
    TH1 *CCDISpi0_g2Comp = dynamic_cast<TH1 *>(histFile->Get("CCDIS_MU_P_PIZERO_PHOTON2_Completeness"));

    if (CCDISpi0_g1Comp && CCDISpi0_g2Comp)
    {
        CCDISpi0_g1Comp->SetLineColor(kBlue);
        CCDISpi0_Leg2.AddEntry(CCDISpi0_g1Comp, "#gamma_{1}");
        CCDISpi0_g1Comp->SetTitleOffset(1.3, "Y");
        CCDISpi0_g1Comp->Draw();
        gPad->SetLogy();
        CCDISpi0_g2Comp->SetLineColor(kOrange + 7);
        CCDISpi0_Leg2.AddEntry(CCDISpi0_g2Comp, "#gamma_{2}");
        CCDISpi0_g2Comp->Draw("same");
        CCDISpi0_Leg2.Draw();
    }

    // Purity
    theCanvas2->cd(5);
    /*TH1* CCDISpi0_muPure = dynamic_cast<TH1*>(histFile->Get("CCDIS_MU_P_PIZERO_MUON_Purity"));
    CCDISpi0_muPure->SetLineColor(kRed);
    CCDISpi0_muPure->Draw();
    TH1* CCDISpi0_pPure = dynamic_cast<TH1*>(histFile->Get("CCDIS_MU_P_PIZERO_PROTON1_Purity"));
    CCDISpi0_pPure->SetLineColor(kBlue);
    CCDISpi0_pPure->Draw("same");*/
    TH1 *CCDISpi0_g1Pure = dynamic_cast<TH1 *>(histFile->Get("CCDIS_MU_P_PIZERO_PHOTON1_Purity"));
    TH1 *CCDISpi0_g2Pure = dynamic_cast<TH1 *>(histFile->Get("CCDIS_MU_P_PIZERO_PHOTON2_Purity"));

    if (CCDISpi0_g1Pure && CCDISpi0_g2Pure)
    {
        CCDISpi0_g1Pure->SetLineColor(kBlue);
        CCDISpi0_g1Pure->SetTitleOffset(1.3, "Y");
        CCDISpi0_g1Pure->Draw();
        gPad->SetLogy();
        CCDISpi0_g2Pure->SetLineColor(kOrange + 7);
        CCDISpi0_g2Pure->Draw("same");
        CCDISpi0_Leg2.Draw();
    }

    // Vtx dR
    theCanvas2->cd(6);
    TH1 *CCDISpi0_muVtxR = dynamic_cast<TH1 *>(histFile->Get("CCDIS_MU_P_PIZERO_VtxDeltaR"));
    if (CCDISpi0_muVtxR)
    {
        //FitCBFun(CCDISpi0_muVtxR, 0.4, 0.15);
        CCDISpi0_muVtxR->SetLineColor(kBlack);
        CCDISpi0_muVtxR->Draw();
    }

    theCanvas2->Print("CCDIS_pmupi0.png");

    plotCorrectFrac(validFileName);
}
