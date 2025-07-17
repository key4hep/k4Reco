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
#ifndef SortingFunctions_hh
#define SortingFunctions_hh 1

/* --------------------------------------------------------------------------
   sorting of hits with respect to their energies
   -------------------------------------------------------------------------- */
// template <class Lhs, class Rhs> inline bool compareByValue(const Lhs& lhs, const Rhs& rhs) {
//   return lhs.second < rhs.second;
// }
template <class T>
inline bool compareByValue(const T& lhs, const T& rhs) {
  return lhs.second < rhs.second;
}

#endif // SortingFunctions_hh
