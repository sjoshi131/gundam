//
// Created by Nadrino on 07/03/2025.
//

/// @file DialManager.cpp
/// @brief Implements DialManager: builds the DialCollection objects from the
///        parameter definitions and keeps their dials up to date.
/// @ingroup DialDictionary
///
/// The public contract of each method is documented in DialManager.h. The
/// comments here describe how each method does its job.

#include "DialManager.h"


/// @par Implementation
/// Loops over the parameter sets in configuration order and appends one
/// DialCollection per dial-set definition. Every new collection gets a
/// pointer to the full list of parameter sets, its own position in
/// `_dialCollectionList_` as index, the index of the set it supervises and,
/// for per-parameter definitions, the index of that parameter. It is then
/// configured with its definition block, which creates the binned dials.
/// Event-by-event dials are only created later, while the data is loaded.
void DialManager::configureImpl(){
  LogThrowIf(_parametersManagerPtr_==nullptr, "Parameters manager not set.");

  // Rebuild from scratch: configure() may run more than once
  // (Propagator::clearContent() reconfigures before each reload).
  _dialCollectionList_.clear();
  for(size_t iParSet = 0 ; iParSet < _parametersManagerPtr_->getParameterSetsList().size() ; iParSet++ ){
    if( not _parametersManagerPtr_->getParameterSetsList()[iParSet].isEnabled() ){ continue; }
    // DEV / DialCollections
    // Set-wide form: the ParameterSet itself lists dial-set definitions.
    // Each one becomes a collection with one dial per parameter of the set,
    // so no supervised parameter index is given.
    if( not _parametersManagerPtr_->getParameterSetsList()[iParSet].getDialSetDefinitions().empty() ){
      for( auto& dialSetDef : _parametersManagerPtr_->getParameterSetsList()[iParSet].getDialSetDefinitions().loop() ){
        _dialCollectionList_.emplace_back(&_parametersManagerPtr_->getParameterSetsList());
        _dialCollectionList_.back().setIndex(int(_dialCollectionList_.size()) - 1);
        _dialCollectionList_.back().setSupervisedParameterSetIndex(int(iParSet) );
        _dialCollectionList_.back().configure( dialSetDef );
      }
    }
    else{
      // Per-parameter form: each enabled parameter carries its own
      // definitions, and each definition becomes a collection whose dials all
      // read that one parameter.

      for( auto& par : _parametersManagerPtr_->getParameterSetsList()[iParSet].getParameterList() ){ // mutable, since it may be disabled below
        if( not par.isEnabled() ){ continue; }

        // Check if no definition is present -> disable the parameter in that case
        if( par.getDialDefinitionsList().empty() ) {
          LogAlert << "Disabling \"" << par.getFullTitle() << "\": no dial definition." << std::endl;
          par.setIsEnabled(false);
          continue;
        }

        for( const auto& dialDefinitionConfig : par.getDialDefinitionsList() ){
          _dialCollectionList_.emplace_back(&_parametersManagerPtr_->getParameterSetsList());
          _dialCollectionList_.back().setIndex(int(_dialCollectionList_.size()) - 1);
          _dialCollectionList_.back().setSupervisedParameterSetIndex(int(iParSet) );
          _dialCollectionList_.back().setSupervisedParameterIndex(par.getParameterIndex() );
          _dialCollectionList_.back().configure( dialDefinitionConfig );
        }
      }
    }
  }
}
/// @par Implementation
/// Delegates to DialCollection::initialize() for each collection, in order.
void DialManager::initializeImpl(){
  for( auto& dialCollection : _dialCollectionList_ ){
    dialCollection.initialize();
  }
}

/// @par Implementation
/// Delegates to DialCollection::invalidateCachedInputBuffers(), which sets
/// every cached input value of the collection to NaN. NaN never compares
/// equal, so the next DialInputBuffer::update() always flags a change.
void DialManager::invalidateInputBuffers(){
  // be extra sure the dial input will request an update
  for( auto& dialCollection : _dialCollectionList_ ){
    dialCollection.invalidateCachedInputBuffers();
  }
}
/// @par Implementation
/// Only event-by-event collections are resized: binned collections have a
/// fixed number of dials set at configuration, so there is nothing to trim.
void DialManager::shrinkDialContainers(){
  LogInfo << "Resizing dial containers..." << std::endl;
  for( auto& dialCollection : _dialCollectionList_ ) {
    if( dialCollection.isEventByEvent() ){ dialCollection.resizeContainers(); }
  }
}
/// @par Implementation
/// Collections are selected by a non-empty dial leaf name rather than by
/// DialCollection::isEventByEvent(), so tabulated and kriged collections keep
/// their dials. The buffers are invalidated afterwards so the next update
/// re-evaluates every dial that remains.
void DialManager::clearEventByEventDials(){
  for( auto& dialCollection: _dialCollectionList_ ) {
    if( not dialCollection.getDialLeafName().empty() ) { dialCollection.clear(); }
  }
  invalidateInputBuffers();
}
/// @par Implementation
/// The scoped timer starts here and stops when the function returns, so the
/// whole loop is timed. Input buffers are refreshed before the callbacks run,
/// because the tabulated and kriged callbacks read the refreshed buffers.
void DialManager::updateDialState(){
  auto s{_updateTimer_.scopeTime()}; // scope guard; stops _updateTimer_ when it goes out of scope
  for( auto& dialCollection : _dialCollectionList_ ) {
    dialCollection.updateInputBuffers();
    dialCollection.update();
  }
}
/// @par Implementation
/// Builds a three-column GenericToolbox::TablePrinter. Collections are
/// already grouped by parameter set (configureImpl() creates them set by
/// set), so a separator is added whenever the supervised set index changes.
void DialManager::printSummaryTable() const{
  GenericToolbox::TablePrinter t; // the table being built, printed at the end
  t << "Dial Collection" << GenericToolbox::TablePrinter::NextColumn;
  t << "Type" << GenericToolbox::TablePrinter::NextColumn;
  t << "Options" << GenericToolbox::TablePrinter::NextLine;
  int parSetIdx{-1}; // parameter set of the previous row; -1 so the first row gets a separator
  for( auto& dialCollection : _dialCollectionList_ ) {
    // Blue background marks event-by-event collections (see the legend below the table).
    if( dialCollection.isEventByEvent() ){ t.setColorBuffer(LOGGER_STR_COLOR_BLUE_BG); }
    else{ t.setColorBuffer(""); }

    // New parameter set: draw a separator line before its first collection.
    if( parSetIdx != dialCollection.getSupervisedParameterSetIndex() ) {
      t.addSeparatorLine();
      parSetIdx = dialCollection.getSupervisedParameterSetIndex();
    }
    t << dialCollection.getTitle() << GenericToolbox::TablePrinter::NextColumn;
    t << dialCollection.getDialType() << GenericToolbox::TablePrinter::NextColumn;
    t << dialCollection.getDialOptions() << GenericToolbox::TablePrinter::NextColumn;
  }

  LogInfo << "Dial manager has " << _dialCollectionList_.size() << " dial collections:" << std::endl;
  t.printTable();
  LogInfo << LOGGER_STR_COLOR_BLUE_BG << "      " << LOGGER_STR_COLOR_RESET << ": Event-by-event dials." << std::endl;
}
