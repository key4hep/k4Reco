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
#include <sstream>
#include <string>
using std::ostringstream;
using std::string;

namespace vertex_lcfi {
namespace util {

  string makeString(const double param) {
    ostringstream os;
    // We don't need this but here it is for rememberance
    // os.setf(ios::fixed);
    // os.precision(2);
    os << param;
    return os.str();
  }

  string makeString(const void* param) {
    ostringstream os;
    os << param;
    return os.str();
  }

  string makeString(const bool param) {
    if (param == 1)
      return "TRUE";
    else
      return "FALSE";
  }

  string makeString(const string param) { return param; }
} // namespace util
} // namespace vertex_lcfi
