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

#ifndef INPLACE_BUKET_SORING_TEST_H
#define INPLACE_BUKET_SORING_TEST_H


// We have some issue here, related to the backup buffer. Need special test to find the bug
void test_bucket_read_write();

void test_inplace_sorting1();
void test_inplace_sorting2();

#endif //INPLACE_BUKET_SORING_TEST_H
