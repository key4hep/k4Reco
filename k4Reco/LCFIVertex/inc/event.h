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
#ifndef LCFIEVENT_H
#define LCFIEVENT_H

#include "../util/inc/matrix.h"
#include "../util/inc/vector3.h"
#include <vector>

namespace vertex_lcfi {
using namespace util;

// Forward Declarations
class Track;
class Jet;
class Vertex;

//! Event
/*!
Description
\author Ben Jeffery (b.jeffery1@physics.ox.ac.uk)
*/

class Event {
public:
  //! Default Constructor
  /*!
  Sets ip to origin and ip error to 10micron spherical
  */
  Event();
  Event(const Event&) = default;
  Event& operator=(const Event&) = default;
  //! Construct with interaction point
  /*!
  \param Position ip position
  \param Error ip error
  */
  Event(const Vector3& Position, const SymMatrix3x3& Error);

  //! Construct with interaction point as vertex
  /*!
  \param Pointer to IP Vertex
  */
  Event(Vertex* ipVertex);

  //! Add Jet
  /*!
  \param Jet Pointer to jet
  */
  void addJet(Jet* Jet);

  //! Get Jets
  /*!
  \return Vector of pointers to jets in the event
  */
  const std::vector<Jet*>& jets() const;

  //! Add Track
  /*!
  \param Track Pointer to track
  */
  void addTrack(Track* Track);

  //! Get Tracks
  /*!
  \return Vector of pointers to tracks in the event
  */
  const std::vector<Track*>& tracks() const;

  //! Replace Primary Vertex
  /*!
  Replace the current Primary Vertex
  */
  inline void replacePrimaryVertex(Vertex* NewPrimary) { _IPVertex = NewPrimary; }

  //! Interaction Point position
  /*!
  \return Vector3 of interaction point position
  */
  const Vector3& interactionPoint() const;

  //! Interaction Point position error
  /*!
  \return SymMatrix3x3 of interaction point position error
  */
  const SymMatrix3x3& interactionPointError() const;

  //! Event Vertex
  /*!
  \return Pointer to this events primary vertex (IP)
  */
  Vertex* ipVertex() const;

private:
  Vertex* _IPVertex = nullptr;
  // Depreceated in favour of holding a complete vertex object
  // Vector3 _IPPosition;
  // SymMatrix3x3 _IPError;
  std::vector<Jet*> _Jets{};
  std::vector<Track*> _Tracks{};
};

} // namespace vertex_lcfi
#endif // LCFIEVENT_H
