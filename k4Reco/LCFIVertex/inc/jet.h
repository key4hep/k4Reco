/*
 * Copyright (c) 2020-2024 Key4hep-Project.
 *
 * This file is part of Key4hep.
 * See https://key4hep.github.io/key4hep-doc/ for further info.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
#ifndef LCFIJET_H
#define LCFIJET_H

#include "../util/inc/recalculate.h"
#include "../util/inc/vector3.h"
#include <map>
#include <vector>

namespace vertex_lcfi {
using namespace util;

// Forward Declarations
class Track;
class Event;
class DecayChain;

//! Simple Jet Class
/*!
Description
\author Ben Jeffery (b.jeffery1@physics.ox.ac.uk)
*/

class Jet {
public:
  //! Default Constructor
  Jet();
  Jet(const Jet&) = default;
  Jet& operator=(const Jet&) = default;
  //! Full Constructor
  /*!
  \param Event Pointer to the event containing this jet
  \param Tracks Vector of pointers to the tracks in this jet
  \return Return
  */
  Jet(Event* Event, const std::vector<Track*>& Tracks, double _Energy, Vector3 Momentum, void* TrackingNum);

  //! Event
  /*!
  \return A pointer to this jets event
  */
  Event* event() const;

  //! Tracks
  /*!
  \return A Vector of pointers to the tracks in the jet
  */
  const std::vector<Track*>& tracks() const;

  //! Tracking Number
  /*!
  \return Integer number
  */
  void* trackingNum() const;

  //! Add Track
  /*!
  Add a track to the jet
  \param Track Pointer to the track to be added
  */
  void addTrack(Track* Track);

  //! Remove Track
  /*!
  Remove a track from the jet
  \param TrackR Pointer to track to be removed
  \return 1 if track was found and removed, 0 if not found
  */
  bool removeTrack(Track* TrackR);

  //! Does the Jet have this Track?
  /*!
  Does this jet have a perticular track?
  \param Track Pointer to track to be checked
  \return 1 if track was found, 0 if not found
  */
  bool hasTrack(Track* Track) const;

  //! Momentum
  /*!
  The average perigee momentum of the tracks in the jet
  \return Vector3 of the momentum
  */
  inline Vector3 momentum() const { return _Momentum; }

  //! Energy
  /*!
  Sum energy of the jet
  \return double of the energy
  */
  inline double energy() const { return _Energy; }

private:
  Event* _Event = nullptr;
  std::vector<vertex_lcfi::Track*> _Tracks{};
  void* _TrackingNum = nullptr;
  double _Energy = 0.0;
  Vector3 _Momentum{};
};

} // namespace vertex_lcfi
#endif // LCFIJET_H
