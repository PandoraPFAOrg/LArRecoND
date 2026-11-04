/**
 *  @file   LArRecoND/include/LArT0Helper.h
 *
 *  @brief  Header file for using t0 information from hits in ND-LAr
 *
 *  $Log: $
 */
#ifndef LAR_T0_HELPER_H
#define LAR_T0_HELPER_H 1

#include "Objects/CaloHit.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace lar_content
{

/**
 *  @brief  LArT0Helper class
 */
class LArT0Helper
{
public:
    /**
     *  @brief  Get the time for input calo hit list, here using median time
     *
     *  @param  inputHitList the list of input CaloHits
     *
     *  @return the median time
     */
    static float GetMedianTime(const pandora::CaloHitList &inputCaloHitList);

};

//------------------------------------------------------------------------------------------------------------------------------------------

float LArT0Helper::GetMedianTime(const pandora::CaloHitList &inputCaloHitList)
{
  std::vector<float> times;

  // Get calo hit times and sort them lowest to highest
  for ( auto const& calohit : inputCaloHitList)
    {
      times.push_back(calohit->GetTime());
    }
  std::sort(times.begin(), times.end());

  if ( size(times) == 0 )
    {
      std::cout << "Empty calo hit list being used in t0 helper." << std::endl;
      return -1.f;
    }
  else if ( size(times) == 1 )
    {
      return times[0];
    }
  // If an even number of hits, average the values on either side of the middle
  else if ( size(times)%2 == 0 )
    {
      float mid_f = float(size(times))/2.0;
      unsigned int mid = (unsigned int)mid_f;
      return (times[mid]+times[mid-1])/2.;
    }
  // If an odd number of hits, grab the middle entry of the vector
  else
    {
      unsigned int mid = (unsigned int)(size(times)/2.0);
      return times[mid];
    }

  return -1.f;
}

//------------------------------------------------------------------------------------------------------------------------------------------

} // namespace lar_content

#endif // #ifndef LAR_T0_HELPER_H
