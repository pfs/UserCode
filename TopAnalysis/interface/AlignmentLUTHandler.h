#ifndef RecoCTPPS_ProtonProducer_AlignmentLUTHandler_h
#define RecoCTPPS_ProtonProducer_AlignmentLUTHandler_h

#include "DataFormats/Provenance/interface/RunID.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "TopLJets2015/TopAnalysis/interface/ProtonReconstruction.h"

#include <regex>
#include <fstream>

namespace CTPPSAlCa
{
  class AlignmentLUTHandler
  {
    public:
      typedef unsigned short fill_t;
      typedef std::map<fill_t,RPAlignmentConstants> LUT;
    public:
      AlignmentLUTHandler( const char* file );

      LUT::const_iterator begin() const { return align_map_.begin(); }
      LUT::const_iterator end() const { return align_map_.end(); }

      const RPAlignmentConstants getAlignmentConstants( const fill_t& fill_num ) const;
      inline bool isValid() const { return valid_; }

    private:
      // FIXME to be replaced by a DB handler!
      void loadConstants( const char* file );
      RPAlignmentConstants& retrieveConstants( const fill_t& fill_num );

      bool valid_;
      LUT align_map_; // fill# -> parameters
      std::regex rgx_hdr_, rgx_algn_;
  };
}

#endif
