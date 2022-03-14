#include <iostream>
#include <fstream>
#include <math.h>
using namespace std;

int main(){

    
    double gamma = 1.4;  //比热比
    double t_current = 0.14;  //求解的时间点

    //initial state of SOD shock tube
    //u_rhs and rho_rhs stands for status of right-hand side
    //u_lhs for left-hand side
    //x>0
    double u_rhs = 0;
    double rho_rhs = 0.124;
    double p_rhs = 0.1;
    double c_rhs {0};  //声速

    //x<0
    double u_lhs {0};
    double rho_lhs {1};
    double p_lhs {1};     
    double c_lhs {0};  //声速

    c_rhs = sqrt(gamma * p_rhs / rho_rhs);
    c_lhs = sqrt(gamma * p_lhs / rho_lhs);

    //计算函数F(p*) = f(p*, p_lhs, rho_lhs) + f(p*, p_rhs, rho_rhs)
    //并计算F(0), F(p_lhs)
    double p_center;  //refer to p* in the slides
    p_center = PCenterSolver();
    
    //计算激波和膨胀波中间的速度
    double u_center;
    u_center = 0.5 * (u_rhs + u_lhs + f(p_center,) + g(p_center));

    //计算接触间断左侧的密度和右侧的密度,同时计算出激波和膨胀波的速度
    double z_rhs {0};  //右侧膨胀波的速度
    double right_shock_x {0};  //右侧激波的位置
    double mid_discontinuity_x {0};  //接触间断的位置
    z_rhs = u_rhs - (p_center - p_rhs) / (rho_rhs * (u_rhs - u_center));
    double rho_right {0};  //接触间断右侧的密度
    rho_right = rho_rhs * (u_rhs - z_rhs) / (u_center - z_rhs);

    double z_lhs_head {0};  //左侧膨胀波波头速度
    double z_lhs_tail {0};  //左侧膨胀波波尾速度
    double rho_left {0};  //接触间断左侧的密度
    z_lhs_head = u_lhs - c_lhs;
    rho_left = rho_lhs * pow(p_center / p_lhs,(1 / gamma));
    z_lhs_tail = u_center - sqrt(gamma * p_center / rho_left);



    //计算左边膨胀波的内部状态值
    double x_current {0};
    double u_inner_expand_wave {0};
    double rho_inner_expand_wave {0};
    double p_inner_expand_wave {0};
    double c_inner_expand_wave {0}; //声速
    c_inner_expand_wave = (gamma - 1) / (gamma + 1) * (u_lhs - x_current / t_current) + 2 / (gamma + 1) * c_lhs;
    u_inner_expand_wave = (x_current / t_current) + c_inner_expand_wave;
    //根据等熵关系得到膨胀波内部的压力和密度
    p_inner_expand_wave = p_lhs * pow((c_inner_expand_wave / c_lhs), 2 * gamma / (gamma -1));
    rho_inner_expand_wave = gamma * p_inner_expand_wave / pow(c_inner_expand_wave, 2);


    //数据的输出
    //密度
    fstream rho_output;
    rho_output.open("rho.txt", ios::app);
    for (int i = 0; i < 50; i++){
        rho_output << rho_inner_expand_wave << "   " << std::endl;
    }

    //压力
    fstream p_output;
    p_output.open("p.txt", ios::app);
    for (int i = 0; i < 50; i++){
        p_output << rho_inner_expand_wave << "   " << std::endl;
    }

    //速度
     fstream u_output;
    u_output.open("u.txt", ios::app);
    for (int i = 0; i < 50; i++){
        u_output << rho_inner_expand_wave << "   " << std::endl;
    }



    return 0;
}

//p_center is applied to get the U_1 - U_2 equation.
//使用二分法计算中间的压力p*，这里变量设置为p_center
double PCenterSolver(u_lhs_, u_rhs_){  
    double p_center_ {2};
    for (int i; i < 1000; i++){
        epsilon_ = 1.0e-8;
        p_center_ = p_center_ / 2;
        if(u_lhs_ - u_rhs_ - (p_center_ - p_rhs_)/((rho_rhs_ * c_rhs_ * sqrt((gamma_rhs_ + 1)/2/gamma_rhs_ 
        * (p_center_ / p_rhs_) + (gamma_rhs_ - 1)/2/gamma_rhs_))) - 2 * c_lhs_/(gamma_lhs_ - 1) * 
        (pow(p_center_ / p_lhs_, ((gamma_lhs_ - 1)/ 2 / gamma_rhs_)) - 1) <= epsilon_){
            break;
        }
                          
    } 

    return p_center_;
}

//计算激波的f值
double f(){
    double f_;
    f_ = (p_center_ - p_rhs_)/((rho_rhs_ * c_rhs_ * sqrt((gamma_rhs_ + 1)/2/gamma_rhs_ 
        * (p_center_ / p_rhs_) + (gamma_rhs_ - 1)/2/gamma_rhs_)));
    return f_;
}

//计算膨胀波的f值
double g(){
    double g_;
    g_ = 2 * c_lhs_/(gamma_lhs_ - 1) * 
        (pow(p_center_ / p_lhs_, ((gamma_lhs_ - 1)/ 2 / gamma_rhs_)) - 1);
    return g_;
}