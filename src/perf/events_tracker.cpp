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

#include "../miner/events_tracker.h"
#include <chrono>
#include <iostream>

// Event to track performance
struct Event {
    std::chrono::time_point<std::chrono::steady_clock> timestamp; // time
    std::string event;   // this event name
    std::string base_event; // parent event name
    std::string thread_name; // thread name (not used)
    std::string comment; // comment (not used much)

    Event() = default;
    Event(const char * _event, const char * _base_event, const std::string & _thread_name, const std::string & _comment) :
            event(_event), base_event(_base_event), thread_name(_thread_name), comment(_comment) {
        timestamp = std::chrono::steady_clock::now();
    }
    Event(const std::string & _thread_name, const std::string & _comment) :
            thread_name(_thread_name), comment(_comment) {
        timestamp = std::chrono::steady_clock::now();
    }
};

// collection of the events
struct EventsTracker {
private:
    std::vector<Event> events;

public:
    EventsTracker() {}

    void addEvent(const char * event, const char * base_event, const std::string & thread_name, const std::string & comment) {
        events.push_back(Event(event, base_event, thread_name, comment));
    }
    void addEvent(const std::string & thread_name, const std::string & comment) {
        events.push_back(Event(thread_name, comment));
    }

    // Report generation. Currently every even with available parant generates a report line
    // We might think about aggreagtion with some exclusions (first event better to be excluded because oif metal warm up)
    void generateReport() {
        std::cout << "Events:" << std::endl;

        for (int i = 0; i < events.size(); i++) {
            const Event & evt = events[i];
            const Event * base = nullptr;
            if (!evt.base_event.empty()) {
                for (int j= i-1; j >= 0; j--) {
                    if (evt.base_event == events[j].event && evt.thread_name == events[j].thread_name) {
                        base = &events[j];
                        break;
                    }
                }
            }
            if (!evt.comment.empty() || base!=nullptr) {
                std::cout << evt.event;
                if (base!=nullptr) {
                    std::cout << " => " << evt.base_event << "  " << std::chrono::duration_cast<std::chrono::microseconds>(evt.timestamp - base->timestamp).count() << "  " << evt.thread_name << "  " << evt.comment << std::flush;
                }
                std::cout << std::endl;
            }
        }
    }

    void reset() {
        events.clear();
    }
};

static EventsTracker eventTracker;

void addEvent(const char * event, const char * base_event, const std::string & thread_name, const std::string & comment) {
    eventTracker.addEvent(event, base_event, thread_name, comment);
}

void addEvent(const std::string & thread_name, const std::string & comment) {
    eventTracker.addEvent(thread_name, comment);
}

void resetReport() {
    eventTracker.reset();
}

void generateReport() {
    eventTracker.generateReport();
}

