This is a port of the DDPlanarDigi processor. The original header and source
files can be found in
https://github.com/iLCSoft/MarlinTrkProcessors/blob/master/source/Digitisers/include/DDPlanarDigiProcessor.h
and
https://github.com/iLCSoft/MarlinTrkProcessors/blob/master/source/Digitisers/src/DDPlanarDigiProcessor.cc

It's very similar to the original one. The (possibly incomplete) list of differences would be:
- Different random number generator and seed (makes it impossible to compare on
  an event-by-event level). GSL was used in the processor, TRandom3 is being
  used in the algorithm use the UniqueIDGenSvc to seed each event.
- The detector instance from DD4hep is obtained with GeoSvc in the algorithm
- The algorithm has an extra option to only use a certain number of bits for the
  cell IDs that are obtained from the hits (see
  https://github.com/key4hep/k4Reco/pull/25)
