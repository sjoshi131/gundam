//
// Created by Adrien Blanchet on 10/04/2023.
//

/// @file RootFormula.h
/// @brief Declares RootFormula, a dial whose response is a ROOT TFormula of
///        the fit parameters.
/// @ingroup DialDictionary

#ifndef GUNDAM_ROOTFORMULA_H
#define GUNDAM_ROOTFORMULA_H

#include "DialBase.h"

#include "TFormula.h"

#include <string>


/// Dial whose response is a ROOT TFormula evaluated at the input parameter
/// values.
///
/// A dial is a response function of one or more fit parameters (see
/// DialBase). DialCollection::makeDial(const ConfigReader&) creates a
/// RootFormula for `dialType: Formula` (or `RootFormula`) collections, with
/// the expression from `dialConfig/formulaStr`. The dial inputs are the
/// formula variables `x[0]`, `x[1]`, ... (ROOT also accepts `x`, `y`, `z`,
/// `t` for the first four). They are the parameters listed in the
/// collection's `dialInputList`, in order; without that option, the single
/// parameter the dial is attached to.
///
/// When the collection has a `binning`, DialCollection creates one dial per
/// bin and sets, in each, the formula parameter named after each bin
/// variable to the bin centre (TFormula::SetParameter() through
/// getFormula()).
///
/// Not supported by Cache::Manager (no GPU implementation).
/// @ingroup DialDictionary
class RootFormula : public DialBase {

public:
  /// Creates a dial with an empty formula. Call setFormulaStr() before
  /// evaluating.
  RootFormula() = default;

  /// Copies this formula dial.
  /// @return A new RootFormula with a copy of the TFormula, including its
  ///         parameter values.
  [[nodiscard]] std::unique_ptr<DialBase> clone() const override { return std::make_unique<RootFormula>(*this); }

  /// Gives the type name of this dial.
  /// @return `"RootFormula"`.
  [[nodiscard]] std::string getDialTypeName() const override { return {"RootFormula"}; }

  /// Evaluates the formula with the dial inputs as its variables.
  /// @param input_ Buffer holding the formula variables `x[0]`, `x[1]`, ...
  ///        in `dialInputList` order.
  /// @return TFormula::EvalPar() of the input values, using the formula
  ///         parameters stored in this dial.
  [[nodiscard]] double evalResponse(const DialInputBuffer& input_) const override;

  /// Describes the formula.
  /// @return `"RootFormula: formula: "` followed by the expanded expression
  ///         (TFormula::GetExpFormula()).
  [[nodiscard]] std::string getSummary() const override;

  /// Parses a formula expression and stores it in the dial.
  /// @param formulaStr_ ROOT TFormula expression. It is also used as the
  ///        TFormula name.
  /// @post Any previous formula and its parameter values are replaced.
  ///       Throws (LogThrowIf()) if ROOT cannot parse the expression.
  void setFormulaStr(const std::string& formulaStr_);

  /// Gives mutable access to the formula, e.g. to set its parameters.
  /// DialCollection uses it to set the bin-centre parameters of binned
  /// formula dials.
  /// @return The stored TFormula.
  TFormula& getFormula(){ return _formula_; }


private:
  /// The formula evaluated by evalResponse().
  /// Set by setFormulaStr().
  TFormula _formula_{};

};


#endif //GUNDAM_ROOTFORMULA_H
