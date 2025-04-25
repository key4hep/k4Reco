#
# Copyright (c) 2020-2024 Key4hep-Project.
#
# This file is part of Key4hep.
# See https://key4hep.github.io/key4hep-doc/ for further info.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
from Configurables import ConformalTracking
def configure_conformal_tracking_steps(tracking: ConformalTracking, parameters: dict):
    """
    Configure the ConformalTracking steps based on the provided parameters.
    """
    tracking.stepCollections = [elem["collections"] for elem in parameters.values()]
    tracking.stepParametersNames = [list(elem["params"].keys()) for elem in parameters.values()]
    tracking.stepParametersValues = [list(elem["params"].values()) for elem in parameters.values()]
    tracking.stepParametersFlags = [elem["flags"] for elem in parameters.values()]
    tracking.stepParametersFunctions = [elem["functions"] for elem in parameters.values()]
