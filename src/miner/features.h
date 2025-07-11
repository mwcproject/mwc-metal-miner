// Copyright 2025 The MWC Developers
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.


#ifndef FEATURES_H
#define FEATURES_H


// metrics make sense for debugging only to detect if data spill over
//#define USE_METRICS

// Normally collapsed needs to be tracked. BUT this feature preventing as from the consistency tracking
#define TRACK_COLLAPSED


// Testing features, METAL_DO_WAITS needs to be defined
//#define STAGE1_TESTS
//#define STAGE2_TESTS
//#define STAGE_TRIM_TESTS
#define STAGE_TRIM_TESTS_START_STEP 0
//#define PHASE2_TESTS

// Related to tests. Validates all the data at every stage. It is very slow, can be used for data corruptions source
//#define  ALL_BUCKETS_DATA_VALIDATION

// Printing of some debug/progress messages
//#define SHOW_TRACING
//#define SHOW_NETWORK

// If defined, don't wait on commands to complete. That feature will disable the tests.
//#define METAL_DO_WAITS

// If defined, will add complete callback so performance events will be emitted
//#define METALL_CALLBACKS_PERF

#ifndef METAL_DO_WAITS
  // no debugging with no wait
  #undef STAGE1_TESTS
  #undef STAGE2_TESTS
  #undef STAGE_TRIM_TESTS
  #undef PHASE2_TESTS
#endif

#endif //FEATURES_H
