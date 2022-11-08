#include "SimViewer.h"
#include "Simulation.h"
#include <thread>

namespace redmax {

// SimViewerTimer
const float SimViewerTimer::CurrentTime() {
    if (!_paused) {
        auto t_start = std::chrono::high_resolution_clock::now();

        int num_steps = max(1, int(_dt * _speed / _sim->_options->_h + 0.5));
        bool done = _sim->advance_viewer_step(num_steps);

        auto t_end = std::chrono::high_resolution_clock::now();
        auto t_span = std::chrono::duration_cast<std::chrono::duration<float>>(t_end - t_start);
        std::this_thread::sleep_for(std::chrono::milliseconds((int)((_dt - t_span.count()) * 1e3f)));

        _current_time += (float)_dt;

        if (done && !_sim->_viewer_options->_infinite) {
            _current_time = -1000000.;
        }
    }

    return _current_time;
}

void SimViewerTimer::reset() {
    _current_time = 0.;
}

void SimViewerTimer::reverse_run_or_pause() {
    _paused = !_paused;
}

void SimViewerTimer::slow_down() {
    _speed /= 1.25;
    std::cerr << "speed = " << _speed << "x" << std::endl;
}

void SimViewerTimer::speed_up() {
    _speed *= 1.25;
    std::cerr << "speed = " << _speed << "x" << std::endl;
}


// SimViewerKeyboardHandler
void SimViewerKeyboardHandler::KeyCallback(const int key, const int action) {
    if (key == GLFW_KEY_SPACE && action == GLFW_PRESS) {
        _sim_viewer->reverse_run_or_pause();
    } else if (key == GLFW_KEY_S && action == GLFW_PRESS) {
        _sim_viewer->slow_down();
    } else if (key == GLFW_KEY_F && action == GLFW_PRESS) {
        _sim_viewer->speed_up();
    }
}

// SimViewer
SimViewer::SimViewer(Simulation* sim) : _sim(sim) {
}

SimViewer::~SimViewer() {
    opengl_viewer::Viewer& viewer = opengl_viewer::Viewer::GetViewer();
    viewer.Cleanup();

    delete _timer;
}

void SimViewer::initialize() {
    opengl_viewer::Viewer& viewer = opengl_viewer::Viewer::GetViewer();
    
    opengl_viewer::Option option;

    option.SetVectorOption("background color", _sim->_viewer_options->_bg_color.cast<float>());
    
    option.SetFloatOption("camera field of view", 60.0f);
    option.SetFloatOption("camera aspect ratio", 4.0f / 3.0f);
    // option.SetFloatOption("camera aspect ratio", 16.0f / 9.0f);
    option.SetVectorOption("camera pos", _sim->_viewer_options->_camera_pos.cast<float>());
    option.SetVectorOption("camera up", _sim->_viewer_options->_camera_up.cast<float>());
    option.SetVectorOption("camera look at", _sim->_viewer_options->_camera_lookat.cast<float>());
    option.SetFloatOption("camera zoom speed", 0.5f);
    option.SetFloatOption("camera pan speed", 0.005f);
    
    option.SetIntOption("height", 768);
    option.SetIntOption("width", 1024);
    // option.SetIntOption("height", 1080);
    // option.SetIntOption("width", 1920);
    option.SetFloatOption("shadow acne bias", 0.005f);
    option.SetFloatOption("shadow sampling angle", 0.57f);
    option.SetIntOption("shadow sampling number", 2);
    option.SetBoolOption("record", _sim->_viewer_options->_record);
    if (_sim->_viewer_options->_record) {
        option.SetStringOption("record folder", _sim->_viewer_options->_record_folder);
    }

    // Set up timer and keyboard callback.
    _timer = new SimViewerTimer(_sim, _sim->_viewer_options->_fps, _sim->_viewer_options->_speed);
    option.SetPointerOption("timer", static_cast<void*>(_timer));

    // Set up keyboard handler
    _keyboard_handler = new SimViewerKeyboardHandler(this);
    option.SetPointerOption("keyboard handler", static_cast<void*>(_keyboard_handler));

    viewer.Initialize(option);

    // A checkerboard.
    const int checker_size = 512, square_size = 32;
    std::vector<std::vector<Eigen::Vector3f>> checker_image(checker_size);
    for (int i = 0; i < checker_size; ++i) {
        checker_image[i] = std::vector<Eigen::Vector3f>(checker_size);
        for (int j = 0; j < checker_size; ++j) {
            // Determine the color of the checker.
            if ((i / square_size - j / square_size) % 2) {
                checker_image[i][j] = Eigen::Vector3f(157.0f, 150.0f, 143.0f) / 255.0f;
            } else {
                checker_image[i][j] = Eigen::Vector3f(216.0f, 208.0f, 197.0f) / 255.0f;
            }
        }
    }
    opengl_viewer::Image checker_texture;
    checker_texture.Initialize(checker_image);

    if (_sim->_viewer_options->_ground) {
        // A simple ground.
        Eigen::Matrix3Xf vertex = (Eigen::Matrix<float, 3, 4>()
            << -4.0f, 4.0f, -4.0f, 4.0f,
            4.0f, 4.0f, -4.0f, -4.0f,
            0.0f, 0.0f, 0.0f, 0.0f).finished();
            // -2.0f, -2.0f, -2.0f, -2.0f).finished();
        Eigen::Matrix4f E_g = _sim->_viewer_options->_E_g.cast<float>();
        for (int i = 0;i < vertex.cols();i++) {
            vertex.col(i) = E_g.topLeftCorner(3, 3) * vertex.col(i) + E_g.topRightCorner(3, 1);
        }

        const Eigen::Matrix3Xi face = (Eigen::Matrix<int, 3, 2>()
            << 0, 0,
            2, 3,
            3, 1).finished();
        const Eigen::Matrix2Xf uv = (Eigen::Matrix<float, 2, 4>()
            << 0.0f, 0.0f, 1.0f, 1.0f,
            0.0f, 1.0f, 0.0f, 1.0f).finished();
        opengl_viewer::Option object_option;
        object_option.SetVectorOption("ambient", 0.7f, 0.7f, 0.7f);
        object_option.SetVectorOption("diffuse", 1.0f, 1.0f, 1.0f);
        object_option.SetVectorOption("specular", 1.0f, 1.0f, 1.0f);
        object_option.SetFloatOption("shininess", 1.5f);
        object_option.SetMatrixOption("uv", uv);
        object_option.SetMatrixOption("texture", checker_texture.rgb_data());
        object_option.SetIntOption("texture row num", checker_size);
        object_option.SetIntOption("texture col num", checker_size);
        object_option.SetStringOption("texture mag filter", "nearest");
        viewer.AddStaticObject(vertex, face, object_option);
    }
    
    // Add point lights.
    opengl_viewer::Option light_option;
    light_option.SetVectorOption("ambient", 0.13f, 0.12f, 0.12f);
    light_option.SetVectorOption("diffuse", 0.37f, 0.35f, 0.34f);
    light_option.SetVectorOption("specular", 0.3f, 0.3f, 0.3f);
    float base_z = (float)_sim->_viewer_options->_E_g(2, 3);
    for (float x: {-0.8f, 0.8f}) {
        for (float y: {-0.8f, 0.8f}) {
            for (float z: {0., 2.}) {
                viewer.AddStaticPointLight(Eigen::Vector3f(x, y, base_z + z), light_option);
            }
        }
    }

    opengl_viewer::Option light_option2;
    light_option2.SetVectorOption("ambient", 0.1f, 0.1f, 0.1f);
    light_option2.SetVectorOption("diffuse", 0.4f, 0.4f, 0.4f);
    light_option2.SetVectorOption("specular", 0.02f, 0.01f, 0.01f);
    const Eigen::Vector3f light_position(0.f, 0.f, 1.5);
    viewer.AddStaticPointLight(light_position, light_option2);

    // Add objects
    std::vector<Matrix3Xf> vertex_list;
    std::vector<Matrix3Xi> face_list;
    std::vector<opengl_viewer::Option> option_list;
    _sim->get_rendering_objects(vertex_list, face_list, option_list, _animator_list);
    for (int i = 0;i < vertex_list.size();i++) {
        viewer.AddDynamicObject(vertex_list[i], face_list[i], _animator_list[i], 
            option_list[i]);
    }
}

void SimViewer::reverse_run_or_pause() {
    if (_timer) {
        _timer->reverse_run_or_pause();
    }
}

void SimViewer::slow_down() {
    if (_timer) {
        _timer->slow_down();
    }
}

void SimViewer::speed_up() {
    if (_timer) {
        _timer->speed_up();
    }
}

void SimViewer::run() {
    _timer->reset();

    opengl_viewer::Viewer& viewer = opengl_viewer::Viewer::GetViewer();

    viewer.Run();

    viewer.Cleanup();
}

}