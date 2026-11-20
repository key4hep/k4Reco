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
#ifndef LCFIMEMMANAGE_H
#define LCFIMEMMANAGE_H

#include <vector>

namespace vertex_lcfi {
//! Base class for all MemoryManagers
class MemoryManagerType {
public:
  //! Empty Destructor
  virtual ~MemoryManagerType() {}
  //! Delete all objects that this MemoryManger has pointers to
  virtual void delAll() = 0;
};

//! MemoryManager Controller - see MemoryManager
/*!
Keeps track of MemoryManagers for each type, and tells them to delete
their objects when delAllObjects is called.
*/
class MetaMemoryManager {
public:
  //! Returns the Event duration singleton instance of the controller
  static MetaMemoryManager* Event();
  //! Returns the Run duration singleton instance of the controller
  static MetaMemoryManager* Run();
  //! Delete all objects of all types held by this instance
  void delAllObjects();
  //! Used by the MemoryManager of each type to alert the controller of its existance
  void registerType(MemoryManagerType* Type);

protected:
  //! Do not use
  MetaMemoryManager();
  //! Do not use
  MetaMemoryManager(const MetaMemoryManager&);
  //! Do not use
  MetaMemoryManager& operator=(const MetaMemoryManager&);

private:
  std::vector<MemoryManagerType*> _Types{};
};

//! Memory management
/*!
Memory in the vertex package is handled on a run and event basis
If you wish to have an object that lasts the length of the event then call:
<br><pre>myType* myObject = new myType(construction parameters);</pre>
<br><pre>MemoryManager<myType>::Event()->registerObject(myObject);</pre>
<br>At the end of the event to free all objects of all types made using the above call:
<br><pre>MetaMemoryManager::Event()->delAllObjects();</pre>
<br>Similarly for run lifetime objects, replacing %Event with Run.
*/
template <class T>
class MemoryManager : public MemoryManagerType {
public:
  //! Destructor - will delete all held objects
  virtual ~MemoryManager();
  //! Returns the Event duration singleton instance of the MemoryManager for type T
  static MemoryManager<T>* Event();
  //! Returns the Run duration singleton instance of the MemoryManager for type T
  static MemoryManager<T>* Run();
  //! Register an object for memory management
  void registerObject(T* pointer);
  //! Delete all objects held by this MemoryManager
  void delAll();
  // Protect the constructor, copy and assignment to prevent usage.
protected:
  //! Do not use
  MemoryManager() {}
  //! Do not use
  MemoryManager(const MemoryManager<T>&) {}
  //! Do not use
  MemoryManager<T>& operator=(const MemoryManager<T>&) { return MemoryManager<T>(); }

private:
  std::vector<T*> _Objects{};
};

template <class T>
MemoryManager<T>::~MemoryManager() {
  // Delete all in case the user hasn't done so
  this->delAll();
}

template <class T>
MemoryManager<T>* MemoryManager<T>::Event() {
  static MemoryManager<T> eventInstance;
  MetaMemoryManager::Event()->registerType(&eventInstance);
  return &eventInstance;
}

template <class T>
MemoryManager<T>* MemoryManager<T>::Run() {
  static MemoryManager<T> runInstance;
  MetaMemoryManager::Run()->registerType(&runInstance);
  return &runInstance;
}

template <class T>
void MemoryManager<T>::registerObject(T* pointer) {
  _Objects.push_back(pointer);
}

template <class T>
void MemoryManager<T>::delAll() {
  for (typename std::vector<T*>::iterator iP = _Objects.begin(); iP != _Objects.end(); ++iP) {
    // std::cout << "Delete object @" << (*iP) << std::endl;
    delete (*iP);
  }
  _Objects.clear();
}

} // namespace vertex_lcfi
#endif // LCFIMEMMANAGE_H
