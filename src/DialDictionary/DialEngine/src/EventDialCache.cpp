//
// Created by Adrien Blanchet on 01/12/2022.
//

/// @file EventDialCache.cpp
/// @brief Implements EventDialCache: filling the indexed cache, turning it
///        into the reference cache, and reweighting one event.
/// @ingroup DialDictionary
///
/// The public contract of each method is documented in EventDialCache.h. The
/// comments here describe how each method does its job.

#include "EventDialCache.h"

#include "Logger.h"


/// @par Implementation
/// Works sample by sample so that each sample's events and entries can be
/// reordered together. Entries are matched to events by position only: the
/// event index written during loading is not used, and is overwritten with
/// the new position after the sort. The reference cache and each entry's
/// dial list are reserved to their final size before being filled, to avoid
/// reallocations. An out-of-range interface index is reported together with
/// the list of collections before the job stops.
/// @warning The sort comparator compares `a.getIndices().dataset` with
///          itself, so that test is always false and events are ordered by
///          TChain entry number only. GenericToolbox::getSortPermutation()
///          uses `std::sort`, which is not stable, so events of one sample
///          with the same entry number (for example from different data
///          sets) end up in an unspecified order.
/// @note `_indexedCache_.clear()` keeps the vector's capacity, so the memory
///       of the indexed cache is not released here.
void EventDialCache::buildReferenceCache( SampleSet& sampleSet_, std::vector<DialCollection>& dialCollectionList_){
  LogInfo << "Building event dial cache..." << std::endl;

  LogInfo << "Indexed cache size: " << _indexedCache_.size() << std::endl;
  LogInfo << "Sorting events in sync with indexed cache..." << std::endl;

  size_t nCacheSlots{0}; // total number of reference-cache entries to create (one per event)
  std::vector<std::vector<IndexedCacheEntry>> sampleIndexCacheList{sampleSet_.getSampleList().size()}; // valid entries, per sample

  {
    LogScopeIndent;
    LogInfo << "Breaking down indexed cache per sample..." << std::endl;
    for( auto& entry : _indexedCache_ ){
      if( entry.event.sampleIndex == size_t(-1) ){ continue; }
      if( entry.event.eventIndex == size_t(-1) ){ continue; }

      // Copy the entry, then keep only its dial slots that point at a dial.
      sampleIndexCacheList[entry.event.sampleIndex].emplace_back( entry );
      sampleIndexCacheList[entry.event.sampleIndex].back().dials.clear();
      for( auto& dial : entry.dials ){
        if( dial.collectionIndex == size_t(-1) ){ continue; }
        if( dial.interfaceIndex == size_t(-1)  ){ continue; }
        sampleIndexCacheList[entry.event.sampleIndex].back().dials.emplace_back(dial);
      }
    }

    LogInfo << "Cleaning up the index cache..." << std::endl;
    _indexedCache_.clear();

    LogInfo << "Performing per sample sorting..." << std::endl;
    int iSample{-1}; // index of the current sample; incremented at the top of each iteration
    for( auto& sample : sampleSet_.getSampleList() ){
      iSample++;

      // Note: the first test compares a with itself (always false), so only
      // the entry number takes part in the ordering.
      auto p = GenericToolbox::getSortPermutation( // sort permutation of the sample's events
          sample.getEventList(), []( const Event& a, const Event& b) {
            if( a.getIndices().dataset < a.getIndices().dataset ){ return true; }
            if( a.getIndices().entry < b.getIndices().entry ){ return true; }
            return false;
          });

      LogThrowIf(
          sampleIndexCacheList[iSample].size() != sample.getEventList().size(),
          std::endl << "MISMATCH cache and event list for sample: #" << sample.getIndex() << " " << sample.getName()
              << std::endl << GET_VAR_NAME_VALUE(sampleIndexCacheList[iSample].size())
              << " <-> " << GET_VAR_NAME_VALUE(sample.getEventList().size())
      );
      nCacheSlots += sampleIndexCacheList[iSample].size();

      // Same permutation on both vectors: entry k keeps matching event k.
      GenericToolbox::applyPermutation( sample.getEventList(), p );
      GenericToolbox::applyPermutation( sampleIndexCacheList[iSample],     p );

      // now update the event indices
      for( size_t iEvent = 0 ; iEvent < sample.getEventList().size() ; iEvent++ ){
        sampleIndexCacheList[iSample][iEvent].event.eventIndex = iEvent;
      }
    }

  }

  // Counts the dial slots that point at a dial (both indices set).
  auto countValidDials = [](std::vector<DialIndexCacheEntry>& dialIndices_){
    return std::count_if(dialIndices_.begin(), dialIndices_.end(),
      []( DialIndexCacheEntry& dialIndex_){
        if( dialIndex_.collectionIndex == size_t(-1) ){ return false; }
        if( dialIndex_.interfaceIndex == size_t(-1) ){ return false; }
        return true;
      });
  };

  LogInfo << "Filling up the " << nCacheSlots << " cache dial with references..." << std::endl;
  _cache_.reserve( nCacheSlots );

  for( auto& sampleIndexCache : sampleIndexCacheList ){
    for( auto& indexCache : sampleIndexCache ){

      _cache_.emplace_back();
      auto& cacheEntry = _cache_.back();

      cacheEntry.event =
          &sampleSet_.getSampleList().at(
              indexCache.event.sampleIndex
          ).getEventList().at(
              indexCache.event.eventIndex
          );

      // make sure we don't need extra allocation while emplace_back
      cacheEntry.dialResponseCacheList.reserve( countValidDials(indexCache.dials) );

      // filling up the dial references
      for( auto& dialIndex : indexCache.dials ){
        if( dialIndex.collectionIndex == size_t(-1) or dialIndex.interfaceIndex == size_t(-1) ){ continue; }

        if( dialIndex.interfaceIndex >= dialCollectionList_.at(dialIndex.collectionIndex).getDialInterfaceList().size() ){
          LogError << "BAD: dialIndex.interfaceIndex >= dialCollectionList_.at(dialIndex.collectionIndex).getDialInterfaceList().size()" << std::endl;
          LogError << GET_VAR_NAME_VALUE(dialIndex.interfaceIndex) << " (+1 for the size)" << std::endl;
          LogError << "Selected collection: " << dialCollectionList_.at(dialIndex.collectionIndex).getSummary() << std::endl;
          LogError << "Nb of defined interfaces: " << dialCollectionList_.at(dialIndex.collectionIndex).getDialInterfaceList().size() << std::endl;

          LogError << "Listing available collections:" << std::endl;
          for( auto& dialCol : dialCollectionList_ ){
            LogDebug << dialCol.getSummary() << std::endl;
          }

          LogThrow("DEV ERROR: Please report this issue to github!! This should not happen");
        }

        cacheEntry.dialResponseCacheList.emplace_back(
            dialCollectionList_.at(dialIndex.collectionIndex)
            .getDialInterfaceList().at(dialIndex.interfaceIndex)
        );
      }
    }
  }

  LogInfo << "Reference cache has been setup." << std::endl;
}
/// @par Implementation
/// A single resize() appends the new entries, each a copy of one prototype
/// entry with unset indices and `nDialsMaxPerEvent_` unset dial slots.
void EventDialCache::allocateCacheEntries( size_t nEvent_, size_t nDialsMaxPerEvent_) {
    _indexedCache_.resize(
        _indexedCache_.size() + nEvent_,
        {{std::size_t(-1),std::size_t(-1)},
         std::vector<DialIndexCacheEntry>(nDialsMaxPerEvent_,
                                          {std::size_t(-1),std::size_t(-1)})} );
}
/// @par Implementation
/// resize() drops the trailing entries; shrink_to_fit() then asks the vector
/// to free the capacity they used.
void EventDialCache::shrinkIndexedCache(){
  _indexedCache_.resize(_fillIndex_+1);
  _indexedCache_.shrink_to_fit();
}
/// @par Implementation
/// Allocates one entry per event, with no dial slot, then claims them with
/// fetchNextCacheEntry() sample by sample and event by event. The entries of
/// each sample are therefore in the same order as its events, as
/// buildReferenceCache() requires.
void EventDialCache::fillCacheEntries(const SampleSet& sampleSet_){
  _fillIndex_ = 0;
  allocateCacheEntries( sampleSet_.getNbOfEvents(), 0 );
  for( size_t iSample = 0 ; iSample < sampleSet_.getSampleList().size() ; iSample++ ){
    for( size_t iEvent = 0 ; iEvent < sampleSet_.getSampleList()[iSample].getEventList().size() ; iEvent++ ){
      auto* entry = fetchNextCacheEntry();
      entry->event.sampleIndex = iSample;
      entry->event.eventIndex = iEvent;
    }
  }
}

/// @par Implementation
/// Returns the address of the element at the fill index, then
/// post-increments the index.
EventDialCache::IndexedCacheEntry* EventDialCache::fetchNextCacheEntry(){
  // Warning warning Will Robinson!
  // This only works IFF the indexed cache is not resized.
  LogThrowIf(_fillIndex_ >= _indexedCache_.size());
  return &_indexedCache_[_fillIndex_++];
}


/// @par Implementation
/// The product is accumulated in a local variable so that the cap applies to
/// the whole product before the event is touched. The weight is reset to the
/// base weight and then multiplied, so repeated calls do not compound.
void EventDialCache::reweightEntry( EventDialCache::CacheEntry& entry_){
  // storing the reweight factor in a temporary buffer
  // this allows to perform capping of the value
  double tempReweight{1};

  // calculate the dial responses
  for( auto& dialResponseCache : entry_.dialResponseCacheList ){
    tempReweight *= dialResponseCache.getResponse();
  }

  // applying event weight cap if defined
  _globalEventReweightCap_.process( tempReweight );

  entry_.event->getWeights().resetCurrentWeight(); // reset to the base weight
  entry_.event->getWeights().current *= tempReweight; // apply the reweight factor
}
