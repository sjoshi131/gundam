//
// Created by Adrien Blanchet on 29/11/2022.
//

/// @file DialInterface.cpp
/// @brief Implements the debug summary of DialInterface; its other methods
///        are inline in the header.
/// @ingroup DialDictionary
///
/// The public contract of each method is documented in DialInterface.h. The
/// comments here describe how each method does its job.

#include "DialInterface.h"
#include "GundamGlobals.h"

#include "Logger.h"


/// @par Implementation
/// Appends the parts in a fixed order: dial type, inputs, bin, supervisor,
/// current response, then (non-shallow only) the dial's own summary. The
/// response is recomputed here, not read from any cache.
std::string DialInterface::getSummary(bool shallow_) const {
  std::stringstream ss;
  ss << _dial_->getDialTypeName() << ":";

  if( _inputBufferPtr_ != nullptr ){
    ss << " input(" << _inputBufferPtr_->getSummary(shallow_) << ")";
  }

  // apply on?
  if( _dialBinRef_ != nullptr ){
    ss << " applyOn(" << _dialBinRef_->getSummary(shallow_) << ")";
  }

  if( _responseSupervisorPtr_ != nullptr ){
    ss << " supervisor(" << _responseSupervisorPtr_->getSummary(shallow_) << ")";
  }

  ss << " response=" << this->evalResponse();

  if( not shallow_ ){
    ss << std::endl;
    ss << _dial_->getSummary();
  }

  return ss.str();
}
