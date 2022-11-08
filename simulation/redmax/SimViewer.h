#pragma once
#include "Common.h"
#include "Utils.h"
#include "opengl_viewer.h"

namespace redmax {

class Simulation;
class SimViewer;

class SimViewerTimer : public opengl_viewer::Timer {
public:
    SimViewerTimer(Simulation* sim, int fps = 25, dtype speed = 1.) : opengl_viewer::Timer() {
        _sim = sim;
        _fps = fps;
        _dt = 1. / fps;
        _speed = speed;
        _paused = false;
    }   

    const float CurrentTime();

    void reset();

    void reverse_run_or_pause();
    void slow_down();
    void speed_up();
         
private:
    Simulation* _sim;
    int _fps;
    dtype _dt;
    dtype _speed;
    bool _paused;
    float _current_time;
};

class SimViewerKeyboardHandler : public opengl_viewer::KeyboardHandler {
public:
    SimViewerKeyboardHandler(SimViewer* sim_viewer) {
        _sim_viewer = sim_viewer;
    }

    void KeyCallback(const int key, const int action);

private:
    SimViewer* _sim_viewer;
};

class SimViewer {
public:
    SimViewer(Simulation* sim);
    ~SimViewer();

    void initialize();

    void run();

private:
    Simulation* _sim;

    SimViewerTimer* _timer;
    opengl_viewer::KeyboardHandler* _keyboard_handler;
    std::vector<opengl_viewer::Animator*> _animator_list;

    void reverse_run_or_pause();
    void slow_down();
    void speed_up();

    friend class SimViewerKeyboardHandler;
};

}