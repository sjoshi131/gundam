//
// Created by Adrien Blanchet on 10/04/2023.
//

/// @file RootFormula.cpp
/// @brief Implements RootFormula: formula parsing, evaluation and summary.
/// @ingroup DialDictionary
///
/// The public contract of each method is documented in RootFormula.h. The
/// comments here describe how each method does its job.

#include "RootFormula.h"


/// @par Implementation
/// Passes the address of the first buffer value as the variable array and
/// no parameter array, so TFormula uses the parameter values stored in it.
double RootFormula::evalResponse(const DialInputBuffer& input_) const {
  return _formula_.EvalPar(&input_.getInputBuffer()[0]);
}


/// @par Implementation
/// Builds a new TFormula whose name and expression are both `formulaStr_`,
/// assigns it to the member, then checks TFormula::IsValid().
void RootFormula::setFormulaStr(const std::string& formulaStr_){
  _formula_ = TFormula(formulaStr_.c_str(), formulaStr_.c_str());
  LogThrowIf(not _formula_.IsValid(), "\"" << formulaStr_ << "\": could not be parsed as formula expression.");
}

/// @par Implementation
/// Streams getDialTypeName() and TFormula::GetExpFormula() into a string.
std::string RootFormula::getSummary() const {
  std::stringstream ss;
  ss << this->getDialTypeName() << ": formula: " << _formula_.GetExpFormula();
  return ss.str();
};
