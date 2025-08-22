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
    double max_delta_tan_phi = 0.5; // 转向角变化率限制 (rad/s), 假设 dt=0.1s
    double max_delta_a = 0.5;     // 加速度变化率限制 (m/s³), 假设 dt=0.1s
    double dt = 0.1;         // 时间步长(s)
    int N = 40;              // 预测步数(4秒)
};

// 车辆状态
struct State {
    double x = 0;
    double y = 0;
    double theta = 0;
    double v = 5.0;

    Vector4d toVector() const {
        return Vector4d(x, y, theta, v);
    }

    static State fromVector(const Vector4d& vec) {
        return {vec(0), vec(1), vec(2), vec(3)};
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
    double w_x = 500.0;      // 轨迹跟踪权重 (x, y)
    double w_theta = 5.0;   // 轨迹跟踪权重 (theta)
    double w_v = 0.0;     // 速度跟踪权重
    double w_u_tan_phi = 10.0; // 控制量大小权重
    double w_u_a = 50.0;        // 控制量大小权重
    double w_du_tan_phi = 0.0;// 控制变化率权重
    double w_du_a = 0.0;      // 控制变化率权重
    double w_lambda = 5.0;  // 松弛变量权重（需要较大值以优先满足硬约束）

public:
    VehicleMPC(VehicleParams p) : params(p) {}

    void setReferenceTrajectory(const MatrixXd& ref) {
        if (ref.rows() != params.N + 1 || ref.cols() != 4) {
            throw invalid_argument("参考轨迹维度错误: 应为 (N+1) x 4");
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
        x_next.v = max(x_next.v, 0.01); // 避免速度为零
        return x_next;
    }

    // 计算动力学模型的雅可比矩阵 A_k = df/dx, B_k = df/du
    // FIX: Replaced ambiguous 'Matrix' with 'Eigen::Matrix'
    void compute_dynamics_jacobians(const State& x, const Control& u, Matrix4d& A, Eigen::Matrix<double, 4, 2>& B) {
        double c = cos(x.theta);
        double s = sin(x.theta);
        
        A << 1.0, 0.0, -x.v * s * params.dt, c * params.dt,
             0.0, 1.0,  x.v * c * params.dt, s * params.dt,
             0.0, 0.0, 1.0, (u.tan_phi / params.L) * params.dt,
             0.0, 0.0, 0.0, 1.0;

        B << 0.0, 0.0,
             0.0, 0.0,
             (x.v / params.L) * params.dt, 0.0,
             0.0, params.dt;
    }

    // 计算SCBF函数值 h(x)
    double computeSCBFValue(const State& x, const Obstacle& obs) {
        double dx = obs.x - x.x;
        double dy = obs.y - x.y;
        double d = sqrt(dx*dx + dy*dy);
        double safe_dist = 1.0 + obs.r; // 包含车辆自身半径和障碍物半径的安全距离
        
        double vA_x = x.v * cos(x.theta);
        double vA_y = x.v * sin(x.theta);
        double v_rel_x = obs.vx - vA_x;
        double v_rel_y = obs.vy - vA_y;
        double v_ij = (v_rel_x * dx + v_rel_y * dy) / max(d, 1e-6);
        
        // h(x) = d - d_safe - v_ij^2 / (2*a_max)
        // 仅当 v_ij > 0 (正在靠近) 时，制动项才有意义
        if (v_ij > 0) {
            return d - safe_dist - (v_ij * v_ij) / (2 * params.max_a);
        } else {
            return d - safe_dist;
        }
    }

    // 计算SCBF函数关于状态x的梯度 ∇h(x)
    Vector4d compute_scbf_gradient(const State& x, const Obstacle& obs) {
        double dx = obs.x - x.x;
        double dy = obs.y - x.y;
        double d_sq = dx*dx + dy*dy;
        double d = sqrt(d_sq);

        double vA_x = x.v * cos(x.theta);
        double vA_y = x.v * sin(x.theta);
        double v_rel_x = obs.vx - vA_x;
        double v_rel_y = obs.vy - vA_y;
        
        double d_inv = 1.0 / max(d, 1e-6);
        double d_cub_inv = d_inv * d_inv * d_inv;
        
        double v_ij = (v_rel_x * dx + v_rel_y * dy) * d_inv;

        Vector4d grad_h = Vector4d::Zero();
        double common_term = 0;

        if (v_ij > 0) {
            common_term = -v_ij / params.max_a;
        }

        //  ∂h/∂x
        double dvij_dx = (-v_rel_x * d_inv) - (v_rel_x * dx + v_rel_y * dy) * d_cub_inv * (-dx);
        grad_h(0) = -dx * d_inv + common_term * dvij_dx;

        //  ∂h/∂y
        double dvij_dy = (-v_rel_y * d_inv) - (v_rel_x * dx + v_rel_y * dy) * d_cub_inv * (-dy);
        grad_h(1) = -dy * d_inv + common_term * dvij_dy;

        //  ∂h/∂theta
        double dvij_dtheta = (x.v * sin(x.theta) * dx + -x.v * cos(x.theta) * dy) * d_inv;
        grad_h(2) = common_term * dvij_dtheta;

        //  ∂h/∂v
        double dvij_dv = (-cos(x.theta) * dx - sin(x.theta) * dy) * d_inv;
        grad_h(3) = common_term * dvij_dv;

        return grad_h;
    }


    // 求解MPC问题
    std::vector<Control> solve(const State& current_state) {
        const int M = obstacles.size();
        const int n_states = 4;
        const int n_controls = 2;
        
        // === 1. 定义优化变量和约束数量 ===
        // 优化变量 Z = [x_0, ..., x_N, u_0, ..., u_{N-1}, λ_0, ..., λ_{M*N-1}]
        const int n_var = n_states * (params.N + 1) + n_controls * params.N + M * params.N;
        
        // 约束: 初始状态 + 动力学 + 控制变化率 + SCBF
        const int n_con = n_states + n_states * params.N + n_controls * (params.N - 1) + M * params.N;
        
        // === 2. 生成用于线性化的参考轨迹 ===
        std::vector<State> linearization_states(params.N + 1);
        std::vector<Control> linearization_controls(params.N);
        linearization_states[0] = current_state;
        for (int k = 0; k < params.N; ++k) {
            // 使用简单的控制策略（例如，朝向参考点）进行前向仿真
            // 这里为了简化，我们假设一个简单的控制输入
            linearization_controls[k].tan_phi = 0;
            linearization_controls[k].a = 0.5;
            linearization_states[k+1] = dynamics(linearization_states[k], linearization_controls[k]);
        }

        // === 3. 构建QP问题 (H, g, A, lb, ub) ===
        // Eigen矩阵最好在堆上分配以避免栈溢出
        MatrixXd H = MatrixXd::Zero(n_var, n_var);
        VectorXd g = VectorXd::Zero(n_var);
        MatrixXd A = MatrixXd::Zero(n_con, n_var);
        VectorXd lb(n_con), ub(n_con);
        
        // --- 3.1 目标函数 H 和 g ---
        // 状态和控制变量的索引
        auto x_idx = [&](int k) { return k * n_states; };
        auto u_idx = [&](int k) { return (params.N + 1) * n_states + k * n_controls; };
        auto lambda_idx = [&](int k, int obs_i) { return (params.N + 1) * n_states + params.N * n_controls + obs_i * params.N + k; };

        // a. 状态跟踪成本
        for (int k = 0; k <= params.N; ++k) {
            Vector4d x_ref = reference_trajectory.row(k);
            H(x_idx(k), x_idx(k)) = w_x;                     // x
            H(x_idx(k) + 1, x_idx(k) + 1) = w_x;             // y
            H(x_idx(k) + 2, x_idx(k) + 2) = w_theta;         // theta
            H(x_idx(k) + 3, x_idx(k) + 3) = w_v;             // v
            g.segment(x_idx(k), n_states) = -H.block(x_idx(k), x_idx(k), n_states, n_states) * x_ref;
        }

        // b. 控制量大小成本
        for (int k = 0; k < params.N; ++k) {
            H(u_idx(k), u_idx(k)) = w_u_tan_phi;
            H(u_idx(k) + 1, u_idx(k) + 1) = w_u_a;
        }

        // c. 控制变化率成本
        for (int k = 0; k < params.N - 1; ++k) {
            H(u_idx(k), u_idx(k)) += w_du_tan_phi;
            H(u_idx(k+1), u_idx(k+1)) += w_du_tan_phi;
            H(u_idx(k), u_idx(k+1)) -= w_du_tan_phi;
            H(u_idx(k+1), u_idx(k)) -= w_du_tan_phi;
            
            H(u_idx(k)+1, u_idx(k)+1) += w_du_a;
            H(u_idx(k+1)+1, u_idx(k+1)+1) += w_du_a;
            H(u_idx(k)+1, u_idx(k+1)+1) -= w_du_a;
            H(u_idx(k+1)+1, u_idx(k)+1) -= w_du_a;
        }

        // d. 松弛变量成本
        for (int i = 0; i < M; ++i) {
            for (int k = 0; k < params.N; ++k) {
                H(lambda_idx(k, i), lambda_idx(k, i)) = w_lambda;
            }
        }

        // --- 3.2 约束矩阵 A, lb, ub ---
        int con_idx = 0;

        // a. 初始状态约束 x_0 = current_state (等式约束)
        A.block(con_idx, x_idx(0), n_states, n_states) = Matrix4d::Identity();
        Vector4d current_state_vec = current_state.toVector();
        lb.segment(con_idx, n_states) = current_state_vec;
        ub.segment(con_idx, n_states) = current_state_vec;
        con_idx += n_states;

        // b. 动力学约束 A_k*x_k + B_k*u_k - x_{k+1} = ... (等式约束)
        for (int k = 0; k < params.N; ++k) {
            State x_bar = linearization_states[k];
            Control u_bar = linearization_controls[k];
            Matrix4d Ak;
            // FIX: Replaced ambiguous 'Matrix' with 'Eigen::Matrix'
            Eigen::Matrix<double, 4, 2> Bk;
            compute_dynamics_jacobians(x_bar, u_bar, Ak, Bk);
            
            A.block(con_idx, x_idx(k), n_states, n_states) = Ak;
            A.block(con_idx, u_idx(k), n_states, n_controls) = Bk;
            A.block(con_idx, x_idx(k+1), n_states, n_states) = -Matrix4d::Identity();
            
            Vector4d rhs = Ak * x_bar.toVector() + Bk * u_bar.toVector() - dynamics(x_bar, u_bar).toVector();
            lb.segment(con_idx, n_states) = -rhs;
            ub.segment(con_idx, n_states) = -rhs;
            con_idx += n_states;
        }
        
        // c. 控制变化率约束
        for (int k = 0; k < params.N - 1; ++k) {
            // tan_phi
            A(con_idx, u_idx(k+1)) = 1.0;
            A(con_idx, u_idx(k)) = -1.0;
            lb(con_idx) = -params.max_delta_tan_phi;
            ub(con_idx) = params.max_delta_tan_phi;
            con_idx++;
            // a
            A(con_idx, u_idx(k+1)+1) = 1.0;
            A(con_idx, u_idx(k)+1) = -1.0;
            lb(con_idx) = -params.max_delta_a;
            ub(con_idx) = params.max_delta_a;
            con_idx++;
        }

        // d. SCBF 障碍物软约束 ∇h*x_k + λ_ik >= ∇h*x_bar - h_bar
        for (int i = 0; i < M; ++i) {
            for (int k = 0; k < params.N; ++k) {
                State x_bar = linearization_states[k];
                Vector4d grad_h = compute_scbf_gradient(x_bar, obstacles[i]);
                double h_val = computeSCBFValue(x_bar, obstacles[i]);
                
                A.block<1, 4>(con_idx, x_idx(k)) = grad_h.transpose();
                A(con_idx, lambda_idx(k, i)) = 1.0; // Slack variable
                
                double rhs = grad_h.transpose() * x_bar.toVector() - h_val;
                lb(con_idx) = rhs;
                ub(con_idx) = 1e10; // 无上界
                con_idx++;
            }
        }
        
        // --- 3.3 变量边界 var_lb, var_ub ---
        VectorXd var_lb(n_var), var_ub(n_var);
        
        // 状态边界 (可以设置，这里设为较宽松)
        for (int k = 0; k <= params.N; ++k) {
            var_lb.segment(x_idx(k), n_states) << -1e10, -1e10, -1e10, 0.0;
            var_ub.segment(x_idx(k), n_states) << 1e10, 1e10, 1e10, 30.0; // e.g., max speed 30m/s
        }
        
        // 控制边界
        for (int k = 0; k < params.N; ++k) {
            var_lb.segment(u_idx(k), n_controls) << -params.max_tan_phi, params.max_decel;
            var_ub.segment(u_idx(k), n_controls) << params.max_tan_phi, params.max_a;
        }

        // 松弛变量边界 (非负)
        for(int i=0; i < M * params.N; ++i){
            var_lb( (params.N+1)*n_states + params.N*n_controls + i ) = 0.0;
            var_ub( (params.N+1)*n_states + params.N*n_controls + i ) = 1e10;
        }


        // === 4. 求解QP问题 ===
        QProblem qp(n_var, n_con);
        Options options;
        options.setToDefault();
        options.printLevel = PL_NONE;
        qp.setOptions(options);
        
        int nWSR = 500; // 迭代次数可以根据问题复杂度调整
        returnValue status = qp.init(H.data(), g.data(), A.data(), 
                                    var_lb.data(), var_ub.data(), 
                                    lb.data(), ub.data(), nWSR);

        if (status != SUCCESSFUL_RETURN) {
            cerr << "QP求解失败，状态码: " << status << endl;
            // 返回默认的紧急制动控制
            std::vector<Control> emergency_controls(params.N, {0.0, params.max_decel});
            return emergency_controls;
        }

        // === 5. 提取控制序列 ===
        VectorXd result(n_var);
        qp.getPrimalSolution(result.data());
        
        std::vector<Control> controls;
        for (int k = 0; k < params.N; ++k) {
            controls.push_back({
                result(u_idx(k)),
                result(u_idx(k) + 1)
            });
        }
        
        return controls;
    }
};

int main() {
    VehicleParams params;
    VehicleMPC mpc(params);
    
    // 参考轨迹（缓慢加速）
    MatrixXd ref(params.N + 1, 4);
    for (int i = 0; i <= params.N; ++i) {
        double t = i * params.dt;
        double v_ref = 5.0;
        ref(i, 0) = v_ref*t;    // 位移
        ref(i, 1) = 0.0;
        ref(i, 2) = 0.0;
        ref(i, 3) = v_ref;
    }
    mpc.setReferenceTrajectory(ref);
    
    // 障碍物设置
    std::vector<Obstacle> obstacles;
    // obstacles.push_back({80.0, 0.5, 0, 0, 0.5}); // x, y, vx, vy, radius
    mpc.setObstacles(obstacles);
    
    // 关键修改：让当前状态与参考轨迹的起点完全一致
    State current_state;
    current_state.x = ref(0, 0);      // 0.0
    current_state.y = ref(0, 1);      // 0.0
    current_state.theta = ref(0, 2);  // 0.0
    current_state.v = ref(0, 3);      // 0.0 <--- 这里是关键

    
    // 在 main 函数中调用 solve 之前
    cout << "--- Initial State Check ---" << endl;
    for(size_t i = 0; i < obstacles.size(); ++i) {
        double h_val = mpc.computeSCBFValue(current_state, obstacles[i]);
        cout << "SCBF value for obstacle " << i << ": " << h_val << endl;
        if (h_val < 0) {
            cout << "WARNING: Initial state is already unsafe!" << endl;
        }
    }
    cout << "--------------------------" << endl;

    std::vector<Control> controls = mpc.solve(current_state);
    
    // 输出结果
    cout << "第一个控制指令:" << endl;
    cout << "tan_phi = " << controls[0].tan_phi << ", "
         << "a = " << controls[0].a << endl;
    
    State x = current_state;
    cout << "\n未来10步的预测轨迹:" << endl;
    for (int i = 0; i < 10; ++i) {
        x = mpc.dynamics(x, controls[i]);
        cout << "时刻 " << (i+1)*params.dt << "s: "
             << "x = " << x.x << "m, y = " << x.y << "m, v = " << x.v << "m/s" << endl;
    }
    
    return 0;
}