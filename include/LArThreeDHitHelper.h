/**
 *  @file   include/LArThreeDHitHelper.h
 *
 *  @brief  Header file for the 3D hit helper.
 *
 *  $Log: $
 */
#ifndef LAR_THREE_D_HIT_HELPER_H
#define LAR_THREE_D_HIT_HELPER_H 1

#include "Objects/CaloHit.h"

#include <map>

namespace lar_content
{

/**
 *  @brief  LArThreeDHitHelper class
 */
class LArThreeDHitHelper
{
public:
    /**
   *  @brief  Default constructor
   */
    LArThreeDHitHelper();

    typedef std::map<const pandora::CaloHit *, std::map<const pandora::HitType, const pandora::CaloHit *>> CaloHitMapping;

    /**
   *  @brief  Get the 2D hits for a list of 3D hits
   *
   *  @param caloHitList3D the list of 3D hits
   *  @param caloHitListU the list of U-plane 2D hits
   *  @param caloHitListV the list of V-plane 2D hits
   *  @param caloHitListW the list of W-plane 2D hits
   *  @param caloHitMapping the mapping between 3D hits and their corresponding
   * 2D hits
   *
   */
    static void GetAllChild2DHits(const pandora::CaloHitList *const caloHitList3D, const pandora::CaloHitList *const caloHitListU,
        const pandora::CaloHitList *const caloHitListV, const pandora::CaloHitList *const caloHitListW, LArThreeDHitHelper::CaloHitMapping &caloHitMapping);
};

//------------------------------------------------------------------------------------------------------------------------------------------

} // namespace lar_content

#endif // #ifndef LAR_THREE_D_HIT_HELPER_H
