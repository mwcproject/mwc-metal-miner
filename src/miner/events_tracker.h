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

#ifndef EVENTS_TRACKER_H
#define EVENTS_TRACKER_H

#include <string>
#include <vector>

// Performance metrics support. Normally should be deactivated for production build.

#ifdef TRACK_EVENTS

void addEvent(const char * event, const char * base_event, const std::string & thread_name, const std::string & comment);
void addEvent(const std::string & thread_name, const std::string & comment);
void resetReport();
void generateReport();

#define REPORT_EVENT(event, base_event, thread_name, comment)   addEvent(event, base_event, thread_name, comment)
#define REPORT_COMMENT(thread_name, comment)   addEvent(thread_name, comment)

#define REPORT_RESET     resetReport();
#define REPORT_GENERATE  generateReport();

#else

#define REPORT_EVENT(event, base_event, thread_name, comment)
#define REPORT_COMMENT(thread_name, comment)

#define REPORT_RESET
#define REPORT_GENERATE

#endif


#endif //EVENTS_TRACKER_H
