#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include "qpOASES.hpp"

using namespace std;
using namespace Eigen;
using namespace qpOASES;

// 车辆参数
struct VehicleParams {
    double L = 2.5;          // 轴距(m)
    double max_tan_phi = 0.5; // 最大转向角正切值
    double max_a = 3.0;      // 最大加速度(m/s²)
    double max_decel = -4.0; // 最大减速度(m/s²)
    double max_delta_u = 1.0;// 控制变化率限制
    double dt = 0.1;         // 时间步长(s)
    int N = 40;              // 预测步数(4秒)
};

// 车辆状态
struct State {
    double x = 0;
    double y = 0;
    double theta = 0;
    double v = 0;

    Vector4d toVector() const {
        return Vector4d(x, y, theta, v);
    }
};

// 控制输入
struct Control {
    double tan_phi = 0;
    double a = 0;

    Vector2d toVector() const {
        return Vector2d(tan_phi, a);
    }
};

// 障碍物信息
struct Obstacle {
    double x;    // 位置x
    double y;    // 位置y
    double vx;   // 速度x分量
    double vy;   // 速度y分量
    double r;    // 半径
};

class VehicleMPC {
private:
    VehicleParams params;
    MatrixXd reference_trajectory;
    std::vector<Obstacle> obstacles;
    
    // 目标函数权重
    double w_x = 50.0;     
    double w_v = 30.0;     
    double w_u = 10.0;     
    double w_du = 20.0;    
    double w_lambda = 500.0;

public:
    VehicleMPC(VehicleParams p) : params(p) {}

    void setReferenceTrajectory(const MatrixXd& ref) {
        if (ref.rows() != params.N || ref.cols() != 4) {
            throw invalid_argument("参考轨迹维度错误: 应为 N x 4");
        }
        reference_trajectory = ref;
    }

    void setObstacles(const std::vector<Obstacle>& obs) {
        obstacles = obs;
    }

    // 车辆动力学模型
    State dynamics(const State& x, const Control& u) {
        State x_next;
        x_next.x = x.x + x.v * cos(x.theta) * params.dt;
        x_next.y = x.y + x.v * sin(x.theta) * params.dt;
        x_next.theta = x.theta + (x.v / params.L) * u.tan_phi * params.dt;
        x_next.v = x.v + u.a * params.dt;
        x_next.v = max(x_next.v, 0.01);
        return x_next;
    }

    // 计算SCBF值
    double computeSCBFValue(const State& x, const Obstacle& obs) {
        double dx = obs.x - x.x;
        double dy = obs.y - x.y;
        double d = sqrt(dx*dx + dy*dy);
        double delta = 0.5 + obs.r;
        
        double vA_x = x.v * cos(x.theta);
        double vA_y = x.v * sin(x.theta);
        double v_rel_x = obs.vx - vA_x;
        double v_rel_y = obs.vy - vA_y;
        double v_ij = (v_rel_x * dx + v_rel_y * dy) / max(d, 1e-6);
        
        return d - delta - (v_ij * v_ij) / (2 * params.max_a);
    }

    // 求解MPC问题
    std::vector<Control> solve(const State& current_state) {
        int M = obstacles.size();
        int n_states = 4;
        int n_controls = 2;
        
        // 优化变量和约束数量
        int n_var = n_controls * params.N + M * params.N;
        int n_con = n_states * (params.N - 1) + 2 * params.N + 2 * (params.N - 1) + M * params.N;
        
        // 初始化QP问题，使用兼容旧版本的设置
        QProblem qp(n_var, n_con);
        Options options;
        // 替换setToRobust()，兼容旧版本qpOASES
        options.setToDefault();
        options.enableRegularisation = BT_TRUE; // 增加数值稳定性
        // options.epsIterRefinement = 1e-6;
        options.printLevel = PL_NONE;
        qp.setOptions(options);
        
        // 预测未来状态
        std::vector<State> predicted_states(params.N);
        predicted_states[0] = current_state;
        for (int k = 1; k < params.N; ++k) {
            Control u_guess;
            u_guess.tan_phi = 0;
            u_guess.a = min(1.0, params.max_a);
            predicted_states[k] = dynamics(predicted_states[k-1], u_guess);
        }

        // 1. 目标函数 H 和 g
        MatrixXd H = MatrixXd::Zero(n_var, n_var);
        VectorXd g = VectorXd::Zero(n_var);
        
        for (int k = 0; k < params.N; ++k) {
            // 状态跟踪成本
            Vector4d x_ref = reference_trajectory.row(k);
            Vector4d x_err = predicted_states[k].toVector() - x_ref;
            int u_idx = k * n_controls;
            g[u_idx] += -2 * w_x * x_err[2] * (predicted_states[k].v / params.L) * params.dt;
            
            // 控制输入成本
            H(u_idx, u_idx) += w_u;
            H(u_idx + 1, u_idx + 1) += w_u;

            // 控制变化率成本
            if (k > 0) {
                int u_prev_idx = (k-1) * n_controls;
                H(u_idx, u_idx) += w_du;
                H(u_idx, u_prev_idx) -= w_du;
                H(u_prev_idx, u_idx) -= w_du;
                H(u_prev_idx, u_prev_idx) += w_du;

                H(u_idx + 1, u_idx + 1) += w_du;
                H(u_idx + 1, u_prev_idx + 1) -= w_du;
                H(u_prev_idx + 1, u_idx + 1) -= w_du;
                H(u_prev_idx + 1, u_prev_idx + 1) += w_du;
            }

            // 松弛变量成本
            for (int i = 0; i < M; ++i) {
                int lambda_idx = n_controls * params.N + i * params.N + k;
                H(lambda_idx, lambda_idx) = w_lambda;
            }
        }

        // 2. 约束矩阵 A, lb, ub
        MatrixXd A = MatrixXd::Zero(n_con, n_var);
        VectorXd lb(n_con), ub(n_con);
        VectorXd var_lb(n_var), var_ub(n_var);
        int con_idx = 0;

        // 2.1 动力学约束
        for (int k = 0; k < params.N - 1; ++k) {
            State x_k = predicted_states[k];
            Control u0; u0.tan_phi = 0; u0.a = 0;
            State x_k1_guess = dynamics(x_k, u0);
            Vector4d x_ref_k1 = reference_trajectory.row(k+1);

            MatrixXd J_u(n_states, n_controls);
            J_u << 0, cos(x_k.theta) * params.dt,
                   0, sin(x_k.theta) * params.dt,
                   (x_k.v / params.L) * params.dt, 0,
                   0, params.dt;

            for (int i = 0; i < n_states; ++i) {
                for (int j = 0; j < n_controls; ++j) {
                    A(con_idx + i, k * n_controls + j) = J_u(i, j);
                }
                double rhs = x_ref_k1[i] - x_k1_guess.toVector()[i];
                lb(con_idx + i) = rhs - 0.1;
                ub(con_idx + i) = rhs + 0.1;
            }
            con_idx += n_states;
        }

        // 2.2 控制输入边界约束
        for (int k = 0; k < params.N; ++k) {
            int u_idx = k * n_controls;
            var_lb[u_idx] = -params.max_tan_phi;
            var_ub[u_idx] = params.max_tan_phi;
            var_lb[u_idx + 1] = params.max_decel;
            var_ub[u_idx + 1] = params.max_a;

            if (k > 0) {
                int u_prev_idx = (k-1) * n_controls;
                A(con_idx, u_idx) = 1;
                A(con_idx, u_prev_idx) = -1;
                lb(con_idx) = -params.max_delta_u;
                ub(con_idx) = params.max_delta_u;
                con_idx++;

                A(con_idx, u_idx + 1) = 1;
                A(con_idx, u_prev_idx + 1) = -1;
                lb(con_idx) = -params.max_delta_u;
                ub(con_idx) = params.max_delta_u;
                con_idx++;
            }
        }

        // 2.3 障碍物约束（软约束）
        for (int k = 0; k < params.N; ++k) {
            for (int i = 0; i < M; ++i) {
                double h_scbf = computeSCBFValue(predicted_states[k], obstacles[i]);
                int lambda_idx = n_controls * params.N + i * params.N + k;

                A(con_idx, lambda_idx) = 1.0;
                lb(con_idx) = max(-h_scbf, 0.0);
                ub(con_idx) = 1e10;
                con_idx++;
            }
        }

        // 2.4 松弛变量非负约束
        for (int i = 0; i < M; ++i) {
            for (int k = 0; k < params.N; ++k) {
                int lambda_idx = n_controls * params.N + i * params.N + k;
                var_lb[lambda_idx] = 0;
                var_ub[lambda_idx] = 1e10;
            }
        }

        // 3. 求解QP问题（修复参数传递错误）
        int nWSR = 2000;
        returnValue status = qp.init(H.data(), g.data(), A.data(), 
                                    var_lb.data(), var_ub.data(), 
                                    lb.data(), ub.data(), nWSR);

        // 若初始求解失败，尝试热启动（修复参数顺序）
        if (status != SUCCESSFUL_RETURN) {
            cout << "尝试热启动求解..." << endl;
            nWSR = 2000; // 重新设置迭代次数
            status = qp.hotstart(H.data(), g.data(), A.data(), 
                                var_lb.data(), var_ub.data(), 
                                nWSR, nullptr, nullptr, nullptr);
        }

        if (status != SUCCESSFUL_RETURN) {
            cerr << "QP求解失败，状态码: " << status << endl;
            // 返回默认控制
            std::vector<Control> default_controls(params.N);
            for (auto& u : default_controls) {
                u.tan_phi = 0;
                u.a = 0.5;
            }
            return default_controls;
        }

        // 提取控制序列
        VectorXd result(n_var);
        qp.getPrimalSolution(result.data());
        
        std::vector<Control> controls;
        for (int k = 0; k < params.N; ++k) {
            controls.push_back({
                result[k * n_controls],
                result[k * n_controls + 1]
            });
        }
        
        return controls;
    }
};

int main() {
    VehicleParams params;
    VehicleMPC mpc(params);
    
    // 参考轨迹（缓慢加速）
    MatrixXd ref(params.N, 4);
    for (int i = 0; i < params.N; ++i) {
        double t = i * params.dt;
        double v_ref = min(3.0, 1.5 * t); // 1.5m/s²加速
        ref(i, 0) = 0.5 * 1.5 * t * t;    // 位移
        ref(i, 1) = 0.0;
        ref(i, 2) = 0.0;
        ref(i, 3) = v_ref;
    }
    mpc.setReferenceTrajectory(ref);
    
    // 障碍物设置
    std::vector<Obstacle> obstacles;
    obstacles.push_back({30.0, 5.0, 0, 0, 0.5});
    mpc.setObstacles(obstacles);
    
    // 当前状态
    State current_state;
    current_state.x = 0;
    current_state.y = 0;
    current_state.theta = 0;
    current_state.v = 5.0;
    
    // 求解MPC
    std::vector<Control> controls = mpc.solve(current_state);
    
    // 输出结果
    cout << "前10个时刻的控制指令:" << endl;
    for (int i = 0; i < 10; ++i) {
        cout << "时刻 " << i*params.dt << "s: "
             << "tan_phi = " << controls[i].tan_phi << ", "
             << "a = " << controls[i].a << endl;
    }
    
    State x = current_state;
    cout << "\n前10步车辆位置:" << endl;
    for (int i = 0; i < 10; ++i) {
        x = mpc.dynamics(x, controls[i]);
        cout << "时刻 " << (i+1)*params.dt << "s: "
             << "x = " << x.x << "m, y = " << x.y << "m, v = " << x.v << "m/s" << endl;
    }
    
    return 0;
}
