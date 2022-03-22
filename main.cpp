#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include <exception>
#include <algorithm>

using namespace std;

//设置一些全局变量
const double GAMMA = 1.4;

typedef struct{
    double rho_;
    double u_;
    double p_;
    double c_;  //auxiliary variables
}PrimitiveVar;
double PCenterSolver( PrimitiveVar &area2, PrimitiveVar &area1);
double ShockF(double p_center, PrimitiveVar &area);
double RarefactionF(double p_center, PrimitiveVar &area);
double SonicSpeed(const PrimitiveVar &ws);

int main(int argc, char *argv[]){

    PrimitiveVar area1 {1, 0, 1}, area2 {0.125, 0, 0.1};  //left and right status
    double  dt {0.01};  //delta t for step length of time
    int step_number {15};  //number of time steps

    //get input data from file
//    try
//    {
//        fstream openfile;
//        openfile.open("input.txt",ios::in);
//        cin >> area1.rho_ >> area1.u_ >> area1.p_;
//        cin >> area2.rho_ >> area2.u_ >> area2.p_;
//        cin >> dt >> step_number;
//        // cin >>
//    }
//    catch (...)
//    {
//        std::cerr << "invalid data!" << endl;
//        return -1;
//    }
//    area1.rho_ = 1;
//    area1.p_ = 1;
//    area1.u_ = 0;


    
    //Time steps
    vector<double> t_sample  {vector<double>(step_number)};
    for(int n =0; n < step_number; n++){
        t_sample[n] = n * dt;
    }
    

    //x-direction coordinates
    double xl {-50}, xr {50};
    int number_of_left_point {250}, number_of_right_point {250};
    double dxl {0}, dxr {0};
    dxl = -xl / number_of_left_point;
    dxr = xr / number_of_right_point;
    int num_of_point {number_of_left_point + number_of_right_point};
    vector<double> x_sample {vector<double>(num_of_point)};
    for (int i = 0; i < num_of_point; i++)
    {
        x_sample[i] = xl + i * dxl;
    }
    for (int i = 0; i < number_of_right_point; i++)
    {
        x_sample[num_of_point - i - 1] = xr - i * dxr;
    }



    //solve p* first, using Newtonian iteration method
    PrimitiveVar local_status {0,0,0,0}, area3 {0,0,0,0},
    area4 {0,0,0,0};
    area3.p_ = area4.p_ = PCenterSolver(area1,area2);

    //calculate rho and u of zone 3 and 4 using p*
    double right_wave_speed {0};
    area3.u_ = area4.u_ = 0.5 * (area2.u_ + area1.u_ + ShockF(area3.p_, area2) + RarefactionF(area3.p_, area1));
    //接触间断左侧密度----3区
    area3.rho_ = area1.rho_ * pow(area2.p_ / area1.p_, 1 / GAMMA );
    right_wave_speed = area2.u_ - (area4.p_ - area2.p_) / (area2.rho_ * (area2.u_ - area4.u_));  //右侧激波速度
    //接触间断右侧密度----4区
    area4.rho_ = area4.rho_ * (area4.u_ - right_wave_speed) / (area4.u_ - right_wave_speed);
    
    //左侧膨胀波波头速度和波尾速度---5区
    double left_wave_speed_head {0}, left_wave_speed_tail {0};
    left_wave_speed_head = area1.u_ - SonicSpeed(area1);
    left_wave_speed_tail = area3.u_ - sqrt(GAMMA * area3.p_ / area1.rho_);

     //location of rarefaction
    double x_rarefaction_tail {0}, x_rarefaction_head {0};
    double t_current {0};
    t_current = dt * step_number;  //当前的时间
    x_rarefaction_tail = - left_wave_speed_tail * t_current;  //当前膨胀波尾波的位置
    x_rarefaction_head = - left_wave_speed_head * t_current;  //当前膨胀波波头的位置
    
    //location of shockwave
    double x_shockwave {0};
    x_shockwave = right_wave_speed * t_current;  //当前激波的位置

    //接触间断的位置
    double x_discontinuity {0};
    x_discontinuity = area3.u_ * t_current;  

    //计算膨胀波内部密度、压强和速度----5区
    PrimitiveVar area5 {0,0,0,0};
    area1.c_ = SonicSpeed(area1);
    for (int i = 0; i < num_of_point; i++) {
        if (x_sample[i] < x_rarefaction_tail & x_sample[i] > x_rarefaction_head) {
            area5.c_ = (GAMMA - 1) / (GAMMA + 1) * (area1.u_ - x_sample[i] / t_current)
            + 2 / (GAMMA +1) * area1.c_;
            area5.u_ = area1.c_ + x_sample[i] / t_current;
            area5.p_ = area1.p_ * pow(area5.c_ / area1.c_, 2 * GAMMA / (GAMMA - 1) );
            area5.u_ = GAMMA * area5.p_ / pow(area5.c_, 2);
        }
            
    }

    // for (int i = 0; i < num_of_point; i++) {
    //     if (x_sample[i] < x_rarefaction_head) {
    //         local_status = area1;
    //     }
    //     else if (x_sample[i] < x_rarefaction_tail) {
    //         local_status = area5;
    //     }
    //     else if (x_sample[i] < x_discontinuity) {
    //         local_status = area3;
    //     }
    //     else if (x_sample[i] < x_shockwave) {
    //         local_status = area4;
    //     }
    //     else {
    //         local_status = area2;
    //     }
    // }

    //output
    ofstream output;
    output.open("output.txt",ios::app);
    //x_ coordinates
    for (int i = 0; i < num_of_point; i++)
    {
        if (x_sample[i] < x_rarefaction_head) {
            local_status = area1;
        }
        else if (x_sample[i] < x_rarefaction_tail) {
            local_status = area5;
        }
        else if (x_sample[i] < x_discontinuity) {
            local_status = area3;
        }
        else if (x_sample[i] < x_shockwave) {
            local_status = area4;
        }
        else {
            local_status = area2;
        }
        output << x_sample[i] << local_status.p_ << endl;
    }
    //pressure, density and velocity
    for (int i = 0; i < num_of_point; i++)
    {
        if (x_sample[i] < x_rarefaction_head) {
            local_status = area1;
        }
        else if (x_sample[i] < x_rarefaction_tail) {
            local_status = area5;
        }
        else if (x_sample[i] < x_discontinuity) {
            local_status = area3;
        }
        else if (x_sample[i] < x_shockwave) {
            local_status = area4;
        }
        else {
            local_status = area2;
        }
        output << x_sample[i] << local_status.p_ << endl;
    }

    for (int i = 0; i < num_of_point; i++)
    {
        if (x_sample[i] < x_rarefaction_head) {
            local_status = area1;
        }
        else if (x_sample[i] < x_rarefaction_tail) {
            local_status = area5;
        }
        else if (x_sample[i] < x_discontinuity) {
            local_status = area3;
        }
        else if (x_sample[i] < x_shockwave) {
            local_status = area4;
        }
        else {
            local_status = area2;
        }
        output << x_sample[i] << local_status.p_ << endl;
    }

    // double GAMMA = 1.4;  //比热比
    
    // double t_current = 0.14;  //求解的时间点

    // //initial state of SOD shock tube
    // //u_rhs and rho_rhs stands for status of right-hand side
    // //u_lhs for left-hand side
    // //c_>0
    // double u_rhs = 0;
    // double rho_rhs = 0.124;
    // double p_rhs = 0.1;
    // double c_rhs {0};  //声速

    // //c_<0
    // double u_lhs {0};
    // double rho_lhs {1};
    // double p_lhs {1};     
    // double c_lhs {0};  //声速

    // c_rhs = sqrt(GAMMA * p_rhs / rho_rhs);
    // c_lhs = sqrt(GAMMA * p_lhs / rho_lhs);

    // //计算函数F(p_*) = f(p_*, p_lhs, rho_lhs) + f(p_*, p_rhs, rho_rhs)
    // //并计算F(0), F(p_lhs)
    // double p_center;  //refer to p_* in the slides
    // p_center = PCenterSolver();
    
    // //计算激波和膨胀波中间的速度
    // double u_center;
    // u_center = 0.5 * (u_rhs + u_lhs + f(p_center,) + g(p_center));

    // //计算接触间断左侧的密度和右侧的密度,同时计算出激波和膨胀波的速度
    // double z_rhs {0};  //右侧膨胀波的速度
    // double right_shock_x {0};  //右侧激波的位置
    // double mid_discontinuity_x {0};  //接触间断的位置
    // z_rhs = u_rhs - (p_center - p_rhs) / (rho_rhs * (u_rhs - u_center));
    // double rho_right {0};  //接触间断右侧的密度
    // rho_right = rho_rhs * (u_rhs - z_rhs) / (u_center - z_rhs);

    // double z_lhs_head {0};  //左侧膨胀波波头速度
    // double z_lhs_tail {0};  //左侧膨胀波波尾速度
    // double rho_left {0};  //接触间断左侧的密度
    // z_lhs_head = u_lhs - c_lhs;
    // rho_left = rho_lhs * pow(p_center / p_lhs,(1 / GAMMA));
    // z_lhs_tail = u_center - sqrt(GAMMA * p_center / rho_left);



    // //计算左边膨胀波的内部状态值
    // double x_current {0};
    // double u_inner_expand_wave {0};
    // double rho_inner_expand_wave {0};
    // double p_inner_expand_wave {0};
    // double c_inner_expand_wave {0}; //声速
    // c_inner_expand_wave = (GAMMA - 1) / (GAMMA + 1) * (u_lhs - x_current / t_current) + 2 / (GAMMA + 1) * c_lhs;
    // u_inner_expand_wave = (x_current / t_current) + c_inner_expand_wave;
    // //根据等熵关系得到膨胀波内部的压力和密度
    // p_inner_expand_wave = p_lhs * pow((c_inner_expand_wave / c_lhs), 2 * GAMMA / (GAMMA -1));
    // rho_inner_expand_wave = GAMMA * p_inner_expand_wave / pow(c_inner_expand_wave, 2);


    // //数据的输出
    // //密度
    // fstream rho_output;
    // rho_output.open("rho_.txt", ios::app);
    // for (int i = 0; i < 50; i++){
    //     rho_output << rho_inner_expand_wave << "   " << std::endl;
    // }

    // //压力
    // fstream p_output;
    // p_output.open("p_.txt", ios::app);
    // for (int i = 0; i < 50; i++){
    //     p_output << rho_inner_expand_wave << "   " << std::endl;
    // }

    // //速度
    //  fstream u_output;
    // u_output.open("u_.txt", ios::app);
    // for (int i = 0; i < 50; i++){
    //     u_output << rho_inner_expand_wave << "   " << std::endl;
    // }



    return 0;
}

// //计算声速
// inline double sound_speed(double p_, double rho_){
//     return sqrt(GAMMA * p_ / rho_);
// }
//
//PCenterSolver is applied to calculate p*, using Newtonian Iteration method.
double PCenterSolver( PrimitiveVar &area2, PrimitiveVar &area1){
    double p_center {0}, p_newton_ {0.0001};  //a temperary var stores p*
    //p_newton stands for P(n+1)-P(n)
    double p_center_last {0};
    double deltap {0.0001};  //delta p to calculate difference of F(p)

    double const p_tol {1e-4};  //tolerance of Newtonian Iteration methods.
    while (abs(p_center - p_center_last) > p_tol){
        p_center_last = p_center;
        p_center = p_center - (ShockF(p_center, area2) + RarefactionF(p_center, area1)) / ((ShockF(p_center + deltap, area2) + RarefactionF(p_center + deltap, area1)) / deltap);
    }
    return  p_center;
}

//计算激波的f值
double ShockF(double p_center, PrimitiveVar &area){
    double shock_f_value {0};
    area.c_ = SonicSpeed(area);
    shock_f_value = (p_center - area.p_)/((area.rho_ * area.c_ * sqrt((GAMMA + 1)/2/GAMMA
                                                                      * (p_center / area.p_) + (GAMMA - 1)/2/GAMMA)));
    return shock_f_value;
}

//计算膨胀波的f值
double RarefactionF(double p_center, PrimitiveVar &area){
    area.c_ = SonicSpeed(area);
    double rare_f_value {0};
    rare_f_value = 2 * area.c_ / area.c_ * area.rho_ * (pow(p_center / area.p_, ((GAMMA - 1)/ 2 / GAMMA)) - 1);
    return rare_f_value;
}

// double DifferationofTotalF(){
//     return (ShockF(p_center + deltap_) + RarefactionF(p_center + deltap))
//      / deltap);
// }

//计算声速
double SonicSpeed(const PrimitiveVar &ws) {
    return sqrt(GAMMA * ws.p_ / ws.rho_ );
}


