/**
 *  @file   larrecond/LArObjects/LArSlice3D.h
 *
 *  @brief  Header file for a slice containing 3D hits.
 *
 *  $Log: $
 */
#ifndef LAR_SLICE_3D_H
#define LAR_SLICE_3D_H 1

#include <vector>

#include "larpandoracontent/LArObjects/LArSlice.h"

namespace lar_content
{

/**
 *  @brief  Slice class
 */
class Slice3D : Slice
{
public:
    pandora::CaloHitList m_caloHitListU;  ///< The TPC_VIEW_U calo hit list
    pandora::CaloHitList m_caloHitListV;  ///< The TPC_VIEW_V calo hit list
    pandora::CaloHitList m_caloHitListW;  ///< The TPC_VIEW_W calo hit list
    pandora::CaloHitList m_caloHitList3D; ///< The TPC_3D calo hit list
};

typedef std::vector<Slice3D> Slice3DList;

} // namespace lar_content

#endif // #ifndef LAR_SLICE_3D_H
