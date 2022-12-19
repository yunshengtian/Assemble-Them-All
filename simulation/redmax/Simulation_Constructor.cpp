/*************************************
 * Constructors for Simulation class
 *************************************/

#include "Simulation.h"
#include "Robot.h"
#include "Body/Body.h"
#include "Joint/Joint.h"
#include "Body/BodyCuboid.h"
#include "Body/BodyMeshObj.h"
#include "Body/BodySDFObj.h"
#include "Body/BodyBVHObj.h"
#include "Body/BodySphere.h"
#include "Body/BodyCylinder.h"
#include "Body/BodyAbstract.h"
#include "Joint/JointFixed.h"
#include "Joint/JointRevolute.h"
#include "Joint/JointPrismatic.h"
#include "Joint/JointFree2D.h"
#include "Joint/JointFree3DEuler.h"
#include "Joint/JointFree3DExp.h"
#include "Joint/JointPlanar.h"
#include "Sensor/Tactile.h"
#include "Sensor/TactileRectArray.h"
#include "Force/ForceGroundContact.h"
#include "Force/ForceCuboidCuboidContact.h"
#include "Force/ForceGeneralPrimitiveContact.h"
#include "Force/ForceGeneralSDFContact.h"
#include "Force/ForceGeneralBVHContact.h"
#include "Actuator/Actuator.h"
#include "Actuator/ActuatorMotor.h"
#include "EndEffector/EndEffector.h"
#include "VirtualObject/VirtualObjectSphere.h"
#include "SimViewer.h"
#include "Utils.h"
#include "Common.h"

namespace redmax {

VectorX Simulation::str_to_eigen(std::string str) {
    std::stringstream iss(str);

    dtype value;
    std::vector<dtype> values;
    while (iss >> value) {
        values.push_back(value);
    }

    VectorX vec(values.size());
    for (int i = 0;i < values.size();i++)
        vec(i) = values[i];
    
    return vec;
}

VectorXi Simulation::str_to_eigen_int(std::string str) {
    std::stringstream iss(str);

    int value;
    std::vector<int> values;
    while (iss >> value) {
        values.push_back(value);
    }

    VectorXi vec(values.size());
    for (int i = 0;i < values.size();i++)
        vec(i) = values[i];
    
    return vec;
}

std::vector<Vector3> Simulation::parse_contact_points(std::string str) {
    std::vector<Vector3> contacts; 
    contacts.clear();
    std::string filename = str;
    if (filename.find(_asset_folder) == std::string::npos) // relative path
        filename = _asset_folder + "//" + filename;
    FILE* fp = fopen(filename.c_str(), "r");
    int n;
    int res = fscanf(fp, "%d", &n);
    for (int i = 0;i < n;i++) {
        float x, y, z;
        res = fscanf(fp, "%f %f %f", &x, &y, &z);
        contacts.push_back(Vector3(x, y, z));
    }
    fclose(fp);
    return contacts;
}

Simulation::Simulation(Simulation::Options *options, std::string name): 
    _options(options), _viewer(nullptr) {
    _name = name;
    _robot = NULL;
    _viewer_options = new ViewerOptions();
}

Simulation::Simulation(Simulation::Options *options, Simulation::ViewerOptions *viewer_options, std::string name): 
    _options(options), _viewer_options(viewer_options), _viewer(nullptr) {
    _name = name;
    _robot = NULL;
}

Simulation::Simulation(std::string xml_file_path, bool verbose) {
    _asset_folder = directory_of(xml_file_path);

    pugi::xml_document doc;
    pugi::xml_parse_result result = doc.load_file(xml_file_path.c_str());
    
    if (!result) {
        throw_error("Input Model File (.xml) is incorrect.");
    }

    _joint_map.clear();

    // parse from file
    int joint_cnt = 0;
    parse_from_xml_file(doc.child("redmax"), doc.child("redmax"), nullptr, joint_cnt, verbose);
    
    // init simulation
    this->init(verbose);
}

Simulation::Simulation(std::string xml_string, std::string asset_folder, bool verbose) {
    _asset_folder = asset_folder;

    pugi::xml_document doc;
    pugi::xml_parse_result result = doc.load_string(xml_string.c_str());

    if (!result) {
        throw_error("Input xml string is incorrect.");
    }

    _joint_map.clear();

    // parse from file
    int joint_cnt = 0;
    parse_from_xml_file(doc.child("redmax"), doc.child("redmax"), nullptr, joint_cnt, verbose);
    
    // init simulation
    this->init(verbose);
}

Joint* Simulation::parse_from_xml_file(pugi::xml_node root, pugi::xml_node node,
                                        Joint* parent_joint, int &joint_cnt, bool verbose) {
    if (verbose) {
        std::cerr << "enter: " << node.name() << std::endl;
    }

    if ((std::string)(node.name()) == "redmax") {
        if (node.attribute("model")) {
            this->_name = node.attribute("model").value();
        }
        // parse options
        this->_options = new Options();
        if (node.child("option")) {
            pugi::xml_node option_node = node.child("option");
            if (option_node.attribute("gravity")) {
                this->_options->_gravity = str_to_eigen(option_node.attribute("gravity").value());
            }
            if (option_node.attribute("timestep")) {
                this->_options->_h = (dtype)option_node.attribute("timestep").as_float();
            }
            if (option_node.attribute("integrator")) {
                this->_options->_integrator = option_node.attribute("integrator").value();
            }
            if (option_node.attribute("unit")) {
                this->_options->_unit = option_node.attribute("unit").value();
            }
        }
        // viewer options
        this->_viewer_options = new ViewerOptions();
        // parse ground
        if (node.child("ground")) {
            pugi::xml_node ground_node = node.child("ground");
            _ground = true;
            this->_viewer_options->_ground = true;
            Vector3 pos = str_to_eigen(ground_node.attribute("pos").value());
            Vector3 normal = str_to_eigen(ground_node.attribute("normal").value());
            Vector3 nz = normal.normalized();
            Eigen::Quaternion<dtype> quat;
            quat.setFromTwoVectors(Vector3::UnitZ(), nz);
            Vector3 nx = quat * Vector3::UnitX();
            Vector3 ny = quat * Vector3::UnitY();
            this->_E_g = Matrix4::Identity();
            this->_E_g.topRightCorner(3, 1) = pos;
            this->_E_g.block(0, 0, 3, 1) = nx;
            this->_E_g.block(0, 1, 3, 1) = ny;
            this->_E_g.block(0, 2, 3, 1) = nz;
            this->_viewer_options->_E_g = this->_E_g;
            if (this->_options->_unit == "cm-g")
                this->_viewer_options->_E_g.topRightCorner(3, 1) /= 10.; // scale for better visualization
            else
                this->_viewer_options->_E_g.topRightCorner(3, 1) *= 10.; // scale for better visualization
        }
        this->_E_g = this->_viewer_options->_E_g;
        if (this->_options->_unit == "cm-g")
            this->_E_g.topRightCorner(3, 1) *= 10.;
        else
            this->_E_g.topRightCorner(3, 1) /= 10.;
        
        // parse rendering
        if (node.child("render")) {
            pugi::xml_node render_node = node.child("render");
            if (render_node.attribute("bg_color")) {
                this->_viewer_options->_bg_color = str_to_eigen(render_node.attribute("bg_color").value());
            }
            if (render_node.attribute("speed")) {
                this->_viewer_options->_speed = (dtype)render_node.attribute("speed").as_float();
            }
        }

        // robot
        Robot* robot = new Robot();
        this->addRobot(robot);

        for (auto child : node.children())
            if ((std::string)(child.name()) == "robot") {
                robot->_root_joints.push_back(parse_from_xml_file(root, child, parent_joint, joint_cnt, verbose));
            }

        // parse actuator
        for (auto actuator : node.children())
            if ((std::string)(actuator.name()) == "actuator") {
                for (auto child : actuator.children()) {
                    if ((std::string)(child.name()) == "motor") {
                        auto it = _joint_map.find(child.attribute("joint").value());
                        if (it == _joint_map.end()) {
                            std::string error_msg = "Actuator joint name error: " + (std::string)(child.attribute("joint").value());
                            throw_error(error_msg);
                        }
                        std::string type = child.attribute("ctrl").value() ;
                        if (type == "force") {
                            Vector2 ctrl_range = str_to_eigen(child.attribute("ctrl_range").value());
                            Actuator* actuator = new ActuatorMotor(it->first + "-motor-force", it->second, ActuatorMotor::ControlMode::FORCE, 
                                                                    ctrl_range(0), ctrl_range(1));
                            robot->add_actuator(actuator);
                        } else if (type == "position") {
                            dtype ctrl_P = (dtype)child.attribute("P").as_float();
                            dtype ctrl_D = (dtype)child.attribute("D").as_float();
                            Vector2 ctrl_range(INT_MIN, INT_MAX);
                            if (child.attribute("ctrl_range"))
                                ctrl_range = str_to_eigen(child.attribute("ctrl_range").value());
                            Actuator* actuator = new ActuatorMotor(it->first + "-motor-position", it->second, ActuatorMotor::ControlMode::POS,
                                                                    ctrl_range(0), ctrl_range(1), ctrl_P, ctrl_D);
                            robot->add_actuator(actuator);
                        }
                    }
                }
            }

        // parse contact
        for (auto contact : node.children()) 
            if ((std::string)(contact.name()) == "contact") {
                for (auto child : contact.children()) {
                     dtype kn = 0.;
                    if (child.attribute("kn")) {
                        kn = (dtype)(child.attribute("kn").as_float());
                    } else if (root.child("default").child(child.name()).attribute("kn")) {
                        kn = (dtype)(root.child("default").child(child.name()).attribute("kn").as_float());
                    }
                    dtype kt = 0.;
                    if (child.attribute("kt")) {
                        kt = (dtype)(child.attribute("kt").as_float());
                    } else if (root.child("default").child(child.name()).attribute("kt")) {
                        kt = (dtype)(root.child("default").child(child.name()).attribute("kt").as_float());
                    }
                    dtype mu = 0.;
                    if (child.attribute("mu")) {
                        mu = (dtype)(child.attribute("mu").as_float());
                    } else if (root.child("default").child(child.name()).attribute("mu")) {
                        mu = (dtype)(root.child("default").child(child.name()).attribute("mu").as_float());
                    }
                    dtype damping = 0.;
                    if (child.attribute("damping")) {
                        damping = (dtype)(child.attribute("damping").as_float());
                    } else if (root.child("default").child(child.name()).attribute("damping")) {
                        damping = (dtype)(root.child("default").child(child.name()).attribute("damping").as_float());
                    }
                    if ((std::string)(child.name()) == "ground_contact") {
                        auto it = _body_map.find(child.attribute("body").value());
                        if (it == _body_map.end()) {
                            std::string error_msg = "Ground contact body name error: " + (std::string)(child.attribute("body").value());
                            throw_error(error_msg);
                        }
                        ForceGroundContact* force = new ForceGroundContact(this, it->second, this->_E_g, kn, kt, mu, damping);
                        robot->add_force(force);
                    } else if ((std::string)(child.name()) == "general_primitive_contact") {
                        auto it1 = _body_map.find(child.attribute("general_body").value());
                        if (it1 == _body_map.end()) {
                            std::string error_msg = "General contact body name error: " + (std::string)(child.attribute("general_body").value());
                            throw_error(error_msg);
                        }
                        auto it2 = _body_map.find(child.attribute("primitive_body").value());
                        if (it2 == _body_map.end()) {
                            std::string error_msg = "Primitive contact body name error: " + (std::string)(child.attribute("primitive_body").value());
                            throw_error(error_msg);
                        }
                        ForceGeneralPrimitiveContact* force = new ForceGeneralPrimitiveContact(this, it1->second, it2->second, kn, kt, mu, damping);
                        robot->add_force(force);
                    } else if ((std::string)(child.name()) == "general_SDF_contact") {
                        auto it1 = _body_map.find(child.attribute("general_body").value());
                        if (it1 == _body_map.end()) {
                            std::string error_msg = "General contact body name error: " + (std::string)(child.attribute("general_body").value());
                            throw_error(error_msg);
                        }
                        auto it2 = _body_map.find(child.attribute("SDF_body").value());
                        if (it2 == _body_map.end()) {
                            std::string error_msg = "SDF contact body name error: " + (std::string)(child.attribute("SDF_body").value());
                            throw_error(error_msg);
                        }
                        ForceGeneralSDFContact* force = new ForceGeneralSDFContact(this, it1->second, it2->second, kn, kt, mu, damping);
                        robot->add_force(force);
                    } else if ((std::string)(child.name()) == "general_BVH_contact") {
                        auto it1 = _body_map.find(child.attribute("general_body").value());
                        if (it1 == _body_map.end()) {
                            std::string error_msg = "General contact body name error: " + (std::string)(child.attribute("general_body").value());
                            throw_error(error_msg);
                        }
                        auto it2 = _body_map.find(child.attribute("BVH_body").value());
                        if (it2 == _body_map.end()) {
                            std::string error_msg = "BVH contact body name error: " + (std::string)(child.attribute("BVH_body").value());
                            throw_error(error_msg);
                        }
                        ForceGeneralBVHContact* force = new ForceGeneralBVHContact(this, it1->second, it2->second, kn, kt, mu, damping);
                        robot->add_force(force);
                    }
                }
            }

        // parser variables
        for (auto variable : node.children())
            if ((std::string)(variable.name()) == "variable") {
                for (auto child : variable.children()) {
                    if ((std::string)(child.name()) == "endeffector") {
                        auto it = _joint_map.find(child.attribute("joint").value());
                        if (it == _joint_map.end()) {
                            std::string error_msg = "Endeffector joint name error: " + (std::string)(child.attribute("joint").value());
                            throw_error(error_msg);
                        }
                        Vector3 pos = str_to_eigen(child.attribute("pos").value());
                        dtype radius = 0.1;
                        if (child.attribute("radius")) {
                            radius = (dtype)(child.attribute("radius").as_float());
                        } else if (root.child("default").child("endeffector").attribute("radius")) {
                            radius = (dtype)(root.child("default").child("endeffector").attribute("radius").as_float());
                        }

                        EndEffector* end_effector = new EndEffector(it->second, pos, radius);

                        if (child.attribute("rgba")) {
                            Vector4 rgba = str_to_eigen(child.attribute("rgba").value());
                            end_effector->set_color(rgba.head(3));
                        } else if (root.child("default").child("endeffector").attribute("rgba")) {
                            Vector4 rgba = str_to_eigen(root.child("default").child("endeffector").attribute("rgba").value());
                            end_effector->set_color(rgba.head(3));
                        }

                        robot->add_end_effector(end_effector);
                    }
                }
            }
        
        // parse virtual objects
        for (auto virtuals : node.children())
            if ((std::string)(virtuals.name()) == "virtual") {
                for (auto child : virtuals.children()) {
                    if ((std::string)(child.name()) == "sphere") {
                        std::string name = child.attribute("name").value();
                        Vector3 pos = str_to_eigen(child.attribute("pos").value());
                        dtype radius = 0.1;
                        if (child.attribute("radius")) {
                            radius = (dtype)(child.attribute("radius").as_float());
                        }
                        Vector4 rgba = (Vector4() << 0., 1., 0., 1.).finished();
                        if (child.attribute("rgba")) {
                            rgba = str_to_eigen(child.attribute("rgba").value());
                        }

                        VirtualObject* virtual_object = new VirtualObjectSphere(name, pos, radius, rgba.head(3));
                        
                        robot->add_virtual_object(virtual_object);
                    }
                }
            }

        return nullptr;
    } else if ((std::string)(node.name()) == "link") {
        bool design_params_1 = false, design_params_2 = false, design_params_3 = false, design_params_4 = false, design_params_5 = false, design_params_6 = false; 
        if (node.attribute("design_params")) {
            int design_params_mask = node.attribute("design_params").as_int();
            if (design_params_mask & (1 << 0)) {
                design_params_1 = true;
            }
            if (design_params_mask & (1 << 1)) {
                design_params_2 = true;
            }
            if (design_params_mask & (1 << 2)) {
                design_params_3 = true;
            }
            if (design_params_mask & (1 << 3)) {
                design_params_4 = true;
            }
            if (design_params_mask & (1 << 4)) {
                design_params_5 = true;
            }
            if (design_params_mask & (1 << 5)) {
                design_params_6 = true;
            }
        }

        // parse joint
        Joint *joint;
        if (node.child("joint")) {
            pugi::xml_node joint_node = node.child("joint");
            std::string type = joint_node.attribute("type").value();

            // extract common parameters
            Vector3 pos = str_to_eigen(joint_node.attribute("pos").value());
            Vector4 quat = str_to_eigen(joint_node.attribute("quat").value());
            Matrix3 R = math::quat2mat(quat);

            Joint::Frame frame = Joint::Frame::LOCAL;
            if (joint_node.attribute("frame")) {
                std::string frame_str = joint_node.attribute("frame").value();
                if (frame_str == "LOCAL")
                    frame = Joint::Frame::LOCAL;
                else if (frame_str == "WORLD")
                    frame = Joint::Frame::WORLD;
                else {
                    std::string error_msg = "Frame type error: " + frame_str;
                    throw_error(error_msg);
                }
            }
            if (type == "revolute") {
                Vector3 axis = str_to_eigen(joint_node.attribute("axis").value());
                joint = new JointRevolute(this, joint_cnt ++, axis, parent_joint, R, pos, frame);
            } else if (type == "fixed") {
                joint = new JointFixed(this, joint_cnt ++, parent_joint, R, pos, frame);
            } else if (type == "prismatic") {
                Vector3 axis = str_to_eigen(joint_node.attribute("axis").value());
                joint = new JointPrismatic(this, joint_cnt ++, axis, parent_joint, R, pos, frame);
            } else if (type == "free2d") {
                joint = new JointFree2D(this, joint_cnt ++, parent_joint, R, pos, frame);
            } else if (type == "free3d" || type == "free3d-euler") {
                joint = new JointFree3DEuler(this, joint_cnt ++, parent_joint, R, pos, JointSphericalEuler::Chart::XYZ, frame);
            } else if (type == "free3d-exp") {
                joint = new JointFree3DExp(this, joint_cnt ++, parent_joint, R, pos, frame);
            } else if (type == "spherical" || type == "spherical-euler") {
                joint = new JointSphericalEuler(this, joint_cnt ++, parent_joint, R, pos, JointSphericalEuler::Chart::XYZ, frame);
            } else if (type == "spherical-exp") {
                joint = new JointSphericalExp(this, joint_cnt ++, parent_joint, R, pos, frame);
            } else if (type == "translational") {
                joint = new JointTranslational(this, joint_cnt ++, parent_joint, R, pos, frame);
            } else if (type == "planar") {
                Vector3 axis0 = str_to_eigen(joint_node.attribute("axis0").value());
                Vector3 axis1 = str_to_eigen(joint_node.attribute("axis1").value());
                joint = new JointPlanar(this, joint_cnt++, axis0, axis1, parent_joint, R, pos, frame);
            } else {
                std::string error_msg = "Joint type error: " + type;
                throw_error(error_msg);
            }

            // parse name
            if (joint_node.attribute("name")) {
                joint->_name = joint_node.attribute("name").value();
                _joint_map[joint->_name] = joint;
            }

            // parse damping
            if (joint_node.attribute("damping")) {
                dtype damping = (dtype)joint_node.attribute("damping").as_float();
                joint->set_damping(damping);
            } else if (root.child("default").child("joint").attribute("damping")) {
                dtype damping = (dtype)root.child("default").child("joint").attribute("damping").as_float();
                joint->set_damping(damping);
            }
            
            // // parse stiffness
            // if (joint_node.attribute("stiffness")) {
            //     dtype stiffness = (dtype)joint_node.attribute("stiffness").as_float();
            //     joint->set_stiffness(stiffness);
            // } else if (root.child("default").child("joint").attribute("stiffness")) {
            //     dtype stiffness = (dtype)root.child("default").child("joint").attribute("stiffness").as_float();
            //     joint->set_stiffness(stiffness);
            // }

            // parse joint limit
            if (joint_node.attribute("lim")) {
                Vector2 joint_limit = str_to_eigen(joint_node.attribute("lim").value());
                dtype joint_limit_stiffness = 0.0;
                if (joint_node.attribute("lim_stiffness")) {
                    joint_limit_stiffness = (dtype)joint_node.attribute("lim_stiffness").as_float();
                } else if (root.child("default").child("joint").attribute("lim_stiffness")) {
                    joint_limit_stiffness = (dtype)root.child("default").child("joint").attribute("lim_stiffness").as_float();
                }
                joint->set_joint_limit(joint_limit(0), joint_limit(1), joint_limit_stiffness);
            }
            // activate design parameters
            if (design_params_1) {
                joint->activate_design_parameters_type_1(true);
            }

            if (design_params_5) {
                joint->activate_design_parameters_type_5(true);
            }
        }

        // parse body
        Body *body;
        if (node.child("body")) {
            pugi::xml_node body_node = node.child("body");
            std::string type = body_node.attribute("type").value();

            // extract common parameters
            Vector3 pos = str_to_eigen(body_node.attribute("pos").value());
            Vector4 quat = str_to_eigen(body_node.attribute("quat").value());
            Matrix3 R = math::quat2mat(quat);
            dtype density = 1.0;
            if (body_node.attribute("density")) {
                density = (dtype)body_node.attribute("density").as_float();
            } else if (root.child("default").child("body").attribute("density")) {
                density = (dtype)root.child("default").child("body").attribute("density").as_float();
            }

            if (type == "cuboid") {
                Vector3 length = str_to_eigen(body_node.attribute("size").value());
                body = new BodyCuboid(this, joint, length, R, pos, density);
            } else if (type == "sphere") {
                dtype radius = (dtype)body_node.attribute("radius").as_float();
                body = new BodySphere(this, joint, radius, R, pos, density);
            } else if (type == "cylinder") {
                dtype radius = (dtype)body_node.attribute("radius").as_float();
                dtype height = (dtype)body_node.attribute("height").as_float();
                body = new BodyCylinder(this, joint, height, radius, R, pos, density);
            } else if (type == "mesh" || type == "SDF" || type == "BVH") {
                BodyMeshObj::TransformType transform_type = BodyMeshObj::TransformType::BODY_TO_JOINT;
                if (body_node.attribute("transform_type")) {
                    std::string transform_type_str = body_node.attribute("transform_type").value();
                    if (transform_type_str == "OBJ_TO_JOINT")
                        transform_type = BodyMeshObj::TransformType::OBJ_TO_JOINT;
                    else if (transform_type_str == "OBJ_TO_WORLD")
                        transform_type = BodyMeshObj::TransformType::OBJ_TO_WOLRD;
                    else if (transform_type_str == "BODY_TO_JOINT")
                        transform_type = BodyMeshObj::TransformType::BODY_TO_JOINT;
                    else {
                        std::string error_msg = "Transform type error: " + transform_type_str;
                        throw_error(error_msg);
                    }
                }
                Vector3 scale = Vector3::Ones();
                if (body_node.attribute("scale")) {
                    scale = str_to_eigen(body_node.attribute("scale").value());
                }
                bool adaptive_sample = false;
                if (body_node.attribute("adaptive_sample")) {
                    adaptive_sample = body_node.attribute("adaptive_sample").as_bool();
                }
                std::string filename = body_node.attribute("filename").value();
                if (filename.find(_asset_folder) == std::string::npos) // relative path
                    filename = _asset_folder + "//" + filename;
                if (type == "mesh") {
                    body = new BodyMeshObj(this, joint, filename, R, pos, transform_type, density, scale, adaptive_sample);
                } else if (type == "SDF" || type == "BVH") {
                    dtype dx = 0.05;
		    if (body_node.attribute("dx")) {
			dx = (dtype)body_node.attribute("dx").as_float();
		    }
		    int res = 20;
		    if (body_node.attribute("res")) {
			res = body_node.attribute("res").as_int();
		    }
                    dtype col_th = 0.01;
                    if (body_node.attribute("col_th")) {
                        col_th = (dtype)body_node.attribute("col_th").as_float();
                    }
                    bool load_sdf = false;
                    if (body_node.attribute("load_sdf")) {
                        load_sdf = body_node.attribute("load_sdf").as_bool();
                    }
                    bool save_sdf = false;
                    if (body_node.attribute("save_sdf")) {
                        save_sdf = body_node.attribute("save_sdf").as_bool();
                    }
                    if (type == "SDF") {
                        body = new BodySDFObj(this, joint, filename, R, pos, dx, res, col_th, transform_type, density, scale, adaptive_sample, load_sdf, save_sdf);
                    } else if (type == "BVH") {
                        std::string filename_BVH_mesh = filename;
                        if (body_node.attribute("BVH_mesh_filename")) {
                            filename_BVH_mesh = body_node.attribute("BVH_mesh_filename").value();
                        }
                        body = new BodyBVHObj(this, joint, filename, filename_BVH_mesh, R, pos, dx, res, col_th, transform_type, density, scale, adaptive_sample, load_sdf, save_sdf);
                    }
                }
                
            } else if (type == "abstract") {
                dtype mass = (dtype)body_node.attribute("mass").as_float();
                VectorX inertia = str_to_eigen(body_node.attribute("inertia").value());

                // parse visual mesh
                std::string visual_mesh_filename = "";
                Vector3 visual_frame_pos = Vector3::Zero();
                Vector4 visual_frame_quat = Vector4(1, 0, 0, 0);
                if (body_node.attribute("mesh")) {
                    visual_mesh_filename = body_node.attribute("mesh").value();
                    if (visual_mesh_filename.find(_asset_folder) == std::string::npos) // relative path
                        visual_mesh_filename = _asset_folder + "//" + visual_mesh_filename;
                } else if (body_node.child("visual") && body_node.child("visual").attribute("mesh")) {
                    visual_mesh_filename = body_node.child("visual").attribute("mesh").value();
                    if (visual_mesh_filename.find(_asset_folder) == std::string::npos) // relative path
                        visual_mesh_filename = _asset_folder + "//" + visual_mesh_filename;
                    if (body_node.child("visual").attribute("pos")) {
                        visual_frame_pos = str_to_eigen(body_node.child("visual").attribute("pos").value());
                    }
                    if (body_node.child("visual").attribute("quat")) {
                        visual_frame_quat = str_to_eigen(body_node.child("visual").attribute("quat").value());
                    }
                }
                
                // parse collision contact points
                Vector3 collision_frame_pos = Vector3::Zero();
                Vector4 collision_frame_quat = Vector4(1, 0, 0, 0);
                std::vector<Vector3> collision_points;
                if (body_node.attribute("contacts")) {
                    collision_points = parse_contact_points(body_node.attribute("contacts").value());
                } else if (body_node.child("collision") && body_node.child("collision").attribute("contacts")) {
                    collision_points = parse_contact_points(body_node.child("collision").attribute("contacts").value());
                    if (body_node.child("collision").attribute("pos")) {
                        collision_frame_pos = str_to_eigen(body_node.child("collision").attribute("pos").value());
                    }
                    if (body_node.child("collision").attribute("quat")) {
                        collision_frame_quat = str_to_eigen(body_node.child("collision").attribute("quat").value());
                    }
                }

                // parse scale of mesh and contact points
                Vector3 scale = Vector3::Ones();
                if (body_node.attribute("scale")) {
                    scale = str_to_eigen(body_node.attribute("scale").value());
                }
                body = new BodyAbstract(this, joint, R, pos, mass, inertia, 
                                        visual_mesh_filename, math::quat2mat(visual_frame_quat), visual_frame_pos, 
                                        collision_points, math::quat2mat(collision_frame_quat), collision_frame_pos, scale);
            } else {
                std::string error_msg = "Body type error: " + type;
                throw_error(error_msg);
            }

            // parse name
            if (body_node.attribute("name")) {
                body->_name = body_node.attribute("name").value();
                _body_map[body->_name] = body;
            }

            // parse rgba
            if (body_node.attribute("texture")) {
                body->set_texture(body_node.attribute("texture").value());
            } else if (body_node.attribute("rgba")) {
                Vector4 rgba = str_to_eigen(body_node.attribute("rgba").value());
                body->set_color(rgba.head(3));
            } else if (root.child("default").child("body").attribute("rgba")) {
                Vector4 rgba = str_to_eigen(root.child("default").child("body").attribute("rgba").value());
                body->set_color(rgba.head(3));
            }

            // activate design parameters
            if (design_params_2) {
                body->activate_design_parameters_type_2(true);
            }
            if (design_params_3) {
                body->activate_design_parameters_type_3(true);
            }
            if (design_params_4) {
                body->activate_design_parameters_type_4(true);
            }
            if (design_params_6) {
                body->activate_design_parameters_type_6(true);
            }
        }

        // parse tactile
        if (node.child("tactile")) {
            pugi::xml_node tactile_node = node.child("tactile");
            std::string type = tactile_node.attribute("type").value();
            std::string name = tactile_node.attribute("name").value();
            if (type == "rect_array") {
                Vector3 rect_pos0 = str_to_eigen(tactile_node.attribute("rect_pos0").value());
                Vector3 rect_pos1 = str_to_eigen(tactile_node.attribute("rect_pos1").value());
                Vector3 axis0 = str_to_eigen(tactile_node.attribute("axis0").value());
                Vector3 axis1 = str_to_eigen(tactile_node.attribute("axis1").value());
                Vector2i resolution = str_to_eigen_int(tactile_node.attribute("resolution").value());
                Tactile* tactile = new TactileRectArray(_robot, body, name, rect_pos0, rect_pos1, axis0, axis1, resolution);
                _robot->add_tactile(tactile);
            } else {
                throw_error("Tactile type " + type + " has not been implemented yet");
            }
        }

        for (auto child : node.children())
            if ((std::string)(child.name()) == "link") {
                parse_from_xml_file(root, child, joint, joint_cnt);
            }

        return joint;
    } else if ((std::string)(node.name()) == "robot") {
        Joint* root_joint;
        for (auto child : node.children())
            if ((std::string)(child.name()) == "link") {
                root_joint = parse_from_xml_file(root, child, parent_joint, joint_cnt);
            }
        return root_joint;
    }

    return nullptr;
}

}
