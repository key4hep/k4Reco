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
#ifndef LCFIALGO_H
#define LCFIALGO_H

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
// For exception
#include "lcio.h"
#include "util/inc/memorymanager.h"

using std::string;

namespace vertex_lcfi {
//! Algorithm interface for decay chain construction or vertexing etc
/*!
Description
*/
template <class INTYPE, class OUTTYPE>
class Algo {

public:
  typedef OUTTYPE OutType;
  typedef INTYPE InType;

  virtual ~Algo() {}

  //! Name
  /*!
  String name of the algorithm
  \return String name
  */
  virtual string name() const = 0;

  //! Parameter Names
  /*!
  A vector of the names of the algorithms parameters
  \return vector of string names
  */
  virtual std::vector<string> parameterNames() const = 0;

  //! Parameter Values
  /*!
  A vector of the values of the algorithms parameters, in the same order as parameter names
  \return vector of string values
  */
  virtual std::vector<string> parameterValues() const = 0;

  //! Set String Parameter
  /*!
  Set a string parameter
  \param Parameter String of parameter name
  \param Value String of parameter value
  */
  virtual void setStringParameter(const string& Parameter, const string& Value) = 0;

  //! Set Double Parameter
  /*!
  Set a double parameter
  \param Parameter String of parameter name
  \param Value double of parameter value
  */
  virtual void setDoubleParameter(const string& Parameter, const double Value) = 0;

  //! Set Pointer Parameter
  /*!
  Set a pointer parameter
  \param Parameter String of parameter name
  \param Value pointer to void
  */
  virtual void setPointerParameter(const string& Parameter, void* Value) = 0;

  //! Run the algorithm on a jet
  /*!
  Calculate the Output of the Algo
  \param Jet Pointer to object to be analysed
  \return Output of the algorithm
  */
  virtual OUTTYPE calculateFor(INTYPE Input) const = 0;

protected:
  void badParameter(std::string Parameter) {
    std::stringstream Msg;
    Msg << this->name() << " does not have parameter " << Parameter << "." << std::endl;
    Msg << "(or parameter was of wrong type)" << std::endl;
    Msg << "Avaliable parameters are:" << std::endl;
    std::vector<std::string> Names = this->parameterNames();
    for (std::vector<std::string>::const_iterator iP = Names.begin(); iP != Names.end(); ++iP)
      Msg << (*iP) << std::endl;
    vertex_lcfi::MetaMemoryManager::Run()->delAllObjects();
    vertex_lcfi::MetaMemoryManager::Event()->delAllObjects();
    // Replace with your systems exception if not LCIO
    throw lcio::Exception(Msg.str());
  }
};
} // namespace vertex_lcfi
#endif // LCFIALGO_H
