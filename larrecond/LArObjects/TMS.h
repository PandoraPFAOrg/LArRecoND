/**
 *  @file   larrecond/LArObjects/TMS.h
 *
 *  @brief  Header file defining the TMS format
 *
 *  $Log: $
 */
#ifndef PANDORA_TMS_H
#define PANDORA_TMS_H 1

#include <iostream>
#include <vector>

#include "TChain.h"
#include "TFile.h"
#include "TROOT.h"

namespace lar_nd_reco
{

class TMS
{
public:
    /**
     *  @brief  Constructor requiring TTree pointer
     *
     *  @param  tree The TTree pointer
     */
    TMS(TTree *tree = nullptr);

    /**
     *  @brief  Destructor
     */
    virtual ~TMS();

    /**
     *  @brief  Get the corresponding data entry
     *
     *  @param  entry The entry integer index
     *
     *  @return Total number of bytes read
     */
    virtual Int_t GetEntry(Long64_t entry);

    /**
     *  @brief  Initialise using the input TTree
     *
     *  @param  tree The input TTree
     */
    virtual void Init(TTree *tree);

    TTree *m_fChain;  ///< pointer to the analyzed TTree or TChain
    Int_t m_fCurrent; ///< current Tree number in a TChain

    // MC truth info
    Int_t m_NTrueHits;
    float m_TrueHitX[100000];
    float m_TrueHitY[100000];
    float m_TrueHitZ[100000];
    float m_TrueHitT[100000];
    float m_TrueHitE[100000];
    float m_TrueHitPE[100000];
    float m_TrueHitPEAfterFibers[100000];
    float m_TrueHitPEAfterFibersLongPath[100000];
    float m_TrueHitPEAfterFibersShortPath[100000];
    TBranch *m_b_NTrueHits;
    TBranch *m_b_TrueHitX;
    TBranch *m_b_TrueHitY;
    TBranch *m_b_TrueHitZ;
    TBranch *m_b_TrueHitT;
    TBranch *m_b_TrueHitE;
    TBranch *m_b_TrueHitPE;
    TBranch *m_b_TrueHitPEAfterFibers;
    TBranch *m_b_TrueHitPEAfterFibersLongPath;
    TBranch *m_b_TrueHitPEAfterFibersShortPath;

    // Reco info
    Int_t m_NRecoHits;
    float m_RecoHitX[100000];
    float m_RecoHitY[100000];
    float m_RecoHitZ[100000];
    float m_RecoHitNotZ[100000];
    float m_RecoHitT[100000];
    float m_RecoHitE[100000];
    float m_RecoHitPE[100000];
    bool m_RecoHitIsPedSupped[100000];
    int m_RecoHitBar[100000];
    int m_RecoHitPlane[100000];
    float m_RecoHitDeadtimeStart[100000];
    float m_RecoHitDeadtimeStop[100000];
    TBranch *m_b_NRecoHits;
    TBranch *m_b_RecoHitX;
    TBranch *m_b_RecoHitY;
    TBranch *m_b_RecoHitZ;
    TBranch *m_b_RecoHitNotZ;
    TBranch *m_b_RecoHitT;
    TBranch *m_b_RecoHitE;
    TBranch *m_b_RecoHitPE;
    TBranch *m_b_RecoHitIsPedSupped;
    TBranch *m_b_RecoHitBar;
    TBranch *m_b_RecoHitPlane;
    TBranch *m_b_RecoHitDeadtimeStart;
    TBranch *m_b_RecoHitDeadtimeStop;
};

inline TMS::TMS(TTree *tree)
{
    if (tree == nullptr)
        std::cout << "Warning: null tree passed to TMS" << std::endl;

    Init(tree);
}

inline TMS::~TMS()
{
}

inline Int_t TMS::GetEntry(Long64_t entry)
{
    // Read contents of entry.
    if (!m_fChain)
        return 0;
    return m_fChain->GetEntry(entry);
}

inline void TMS::Init(TTree *tree)
{
    // The Init() function is called when the selector needs to initialize
    // a new tree or chain. Typically here the branch addresses and branch
    // pointers of the tree will be set.
    // It is normally not necessary to make changes to the generated
    // code, but the routine can be extended by the user if needed.
    // Init() will be called many times when running on PROOF
    // (once per file to be processed).

    // Set to default values
    m_fChain = nullptr;
    m_fCurrent = -1;

    // Set branch addresses and branch pointers
    if (!tree)
        return;

    std::cout << "Calling TMS::Init for TTree " << tree->GetName() << std::endl;

    m_fChain = tree;
    m_fChain->SetMakeClass(1);

    m_fChain->SetBranchAddress("NTrueHits", &m_NTrueHits, &m_b_NTrueHits);
    m_fChain->SetBranchAddress("TrueHitX", m_TrueHitX, &m_b_TrueHitX);
    m_fChain->SetBranchAddress("TrueHitY", m_TrueHitY, &m_b_TrueHitY);
    m_fChain->SetBranchAddress("TrueHitZ", m_TrueHitZ, &m_b_TrueHitZ);
    m_fChain->SetBranchAddress("TrueHitT", m_TrueHitT, &m_b_TrueHitT);
    m_fChain->SetBranchAddress("TrueHitE", m_TrueHitE, &m_b_TrueHitE);
    m_fChain->SetBranchAddress("TrueHitPE", m_TrueHitPE, &m_b_TrueHitPE);
    m_fChain->SetBranchAddress("TrueHitPEAfterFibers", m_TrueHitPEAfterFibers, &m_b_TrueHitPEAfterFibers);
    m_fChain->SetBranchAddress("TrueHitPEAfterFibersLongPath", m_TrueHitPEAfterFibersLongPath, &m_b_TrueHitPEAfterFibersLongPath);
    m_fChain->SetBranchAddress("TrueHitPEAfterFibersShortPath", m_TrueHitPEAfterFibersShortPath, &m_b_TrueHitPEAfterFibersShortPath);
    m_fChain->SetBranchAddress("NRecoHits", &m_NRecoHits, &m_b_NRecoHits);
    m_fChain->SetBranchAddress("RecoHitX", m_RecoHitX, &m_b_RecoHitX);
    m_fChain->SetBranchAddress("RecoHitY", m_RecoHitY, &m_b_RecoHitY);
    m_fChain->SetBranchAddress("RecoHitZ", m_RecoHitZ, &m_b_RecoHitZ);
    m_fChain->SetBranchAddress("RecoHitNotZ", m_RecoHitNotZ, &m_b_RecoHitNotZ);
    m_fChain->SetBranchAddress("RecoHitT", m_RecoHitT, &m_b_RecoHitT);
    m_fChain->SetBranchAddress("RecoHitE", m_RecoHitE, &m_b_RecoHitE);
    m_fChain->SetBranchAddress("RecoHitPE", m_RecoHitPE, &m_b_RecoHitPE);
    m_fChain->SetBranchAddress("RecoHitIsPedSupped", m_RecoHitIsPedSupped, &m_b_RecoHitIsPedSupped);
    m_fChain->SetBranchAddress("RecoHitBar", m_RecoHitBar, &m_b_RecoHitBar);
    m_fChain->SetBranchAddress("RecoHitPlane", m_RecoHitPlane, &m_b_RecoHitPlane);
    m_fChain->SetBranchAddress("RecoHitDeadtimeStart", m_RecoHitDeadtimeStart, &m_b_RecoHitDeadtimeStart);
    m_fChain->SetBranchAddress("RecoHitDeadtimeStop", m_RecoHitDeadtimeStop, &m_b_RecoHitDeadtimeStop);
}

} // namespace lar_nd_reco

#endif
