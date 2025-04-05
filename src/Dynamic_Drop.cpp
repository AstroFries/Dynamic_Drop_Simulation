#define MATPLOTLIBCPP_ENABLE_QT5
#include <iostream>
#include <cmath>
#include <cstdio>
#include <vector>
#include <chrono>
#include <thread>
#include <fstream>
#include <string>
#include <utility>
#include <Eigen/Dense>
#include <Eigen/SparseLU>
#include <ctime>
#include <random>
//#include "matplotlibcpp.h"
#include <SFML/Graphics.hpp>
#include <SFML/Audio.hpp>
//using namespace std;
double R = 1.6e-3;
const double pi = 3.1416;
const double g = 9.8;
const double sigma = 7.3e-2;
const double rho = 1e3;
const double k = -g * R * 0.05;
const double deltax = R*0.0555;
const double deltat_0 = 2*1e-7;
double deltat = deltat_0;
double v0 = 16*1e-2;
const double eta = 9.0e-4;
const double l_max = 5e-4;//if l > l_max ,devide
const double V_min = 25*1e-13;
const int l = 3;
double t;
class slice{
  public:
    double V;
    double r;
    double v;
    double x;
    double l;
};
std::vector<slice> slices;
void slice_init(){
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0, v0 * 0.05);
    for (int i = 0; double(i) * deltax < R; ++i){
        slice current_slice;
        current_slice.V = pi * (R * R - double(i) * deltax * double(i) * deltax) * deltax;
        double randomDouble = dis(gen);
        current_slice.v = v0*(1+double(i) * deltax/R) + randomDouble;
        current_slice.r = sqrt(R * R - double(i) * deltax * double(i) * deltax);
        current_slice.x = deltax * (double)i;
        slices.push_back(current_slice);
    }
    
    slices[slices.size() - 1].V = 0;
    slices[slices.size() - 2].V *= 0.33;
}
std::string dataFilePath = std::string(DATA_DIR) + "/workspace.txt";
std::string initFilePath = std::string(DATA_DIR) + "/init.txt";
std::ofstream outFile;
bool read(){
    //data struct : r,delta_t0,v0
    bool match = 0;
    double r_, delta_t_0_, v0_;
    std::ifstream file(dataFilePath);
    if (!file.is_open()) {
        std::cout << "无法打开文件: " << dataFilePath <<" ,已经新建。" <<std::endl;
        //file.close();
        outFile.open(dataFilePath);
        //file.open(dataFilePath);
        outFile << R <<' ' <<  deltat_0 << ' '<< v0 << std::endl;
        outFile.close();
        //slice_init();
        //return 1;
        //if (!file.is_open()) return 1;
    }
    else{
        
        if (!(file >> r_ ) || !(file >> delta_t_0_) || !(file >> v0_) ||
            abs((r_ - R)/R) > 0.001 || abs((delta_t_0_ - deltat_0)/R) > 0.001 || abs((v0_ - v0)/v0) > 0.01)
        {
            std::cout << "与文件: " << dataFilePath <<" 参数不匹配，将新建。" <<std::endl;
            outFile.open(dataFilePath);
            outFile << R <<' ' <<  deltat_0 << ' '<< v0 << std::endl;
            outFile.close();
            file.close();
            return 1;
        }
        std::cout << "成功检查文件: " << dataFilePath <<" ,参数相匹配。" <<std::endl;
        match = 1;
        file.close();
    }
    match = 1;
    if (match){
        file.open(initFilePath);
        if (!file.is_open()) {
            std::cerr << "无法打开文件: " << initFilePath <<" ,将按程序内初始化模拟。" <<std::endl;
            slice_init();
        }else {
            int n;
            if ((file >> n)){
                
            
                for (int i = 0; i< n; ++i){
                    slice c_slice;
                    file >> c_slice.V >> c_slice.v >> c_slice.x;
                    slices.push_back(c_slice);
                }
            }else {
                slice_init();
            }
            file.close();
        }
    }
    outFile.open(dataFilePath, std::ios::app);
    return 0;
}
void flow(){
    if (slices.empty())return;
    if (slices[0].x > deltax){
        slice current_slice;
        current_slice.V = pi * R * R * deltax;
        current_slice.v = v0;
        current_slice.r = R;
        current_slice.x = slices[0].x - deltax;
        slices.insert(slices.begin(), current_slice);
    }
}
double h_c(double r_1, double r_2){
    return (0.75* r_1 * r_1 + 0.5 * r_1 * r_2 + 0.25 * r_2 * r_2)/(r_1 * r_1 + r_1 * r_2 + r_2 * r_2);
}
void compute(){
    t += deltat;
    flow();
    //slices[1].V = slices[2].V;
    //slices[1].V = 0.5*(slices[0].V + slices[2].V);
    for (int i = 0; i < slices.size() - 1; ++i){
        slices[i].r = sqrt(slices[i].V / pi / (slices[i+1].x - slices[i].x));
        slices[i].l = sqrt(pow(slices[i+1].x - slices[i].x,2) + pow(slices[i+1].r - slices[i].r,2));
    }
    for (int i = 0; i < slices.size() - 1; ++i){
        if (slices[i].V < V_min){
            double nV = 0;
            for (int j = i; j < slices.size();j++){
                nV += slices[j].V;
            }
            if (nV * 1e6 > 1e-4)std::cout << t <<":" << nV * 1e6 <<"mL" << std::endl;
            slices.erase(slices.begin() + i + 1, slices.end());
            deltat = deltat_0;
            break;
        }
    }
    for (int i = 0; i + 2 < slices.size();) {
        if (slices[i+1].x - slices[i].x > l_max) {
            //std::cout << slices[i].l << std::endl;
            slice add_slice;
            add_slice.x = 0.5 * slices[i].x + 0.5 * slices[i+1].x;
            //add_slice.r =  slices[i].r;
            add_slice.r =  0.5 * slices[i].r * 0.5 * slices[i].r; 
            add_slice.v = 0.5 * slices[i].v + 0.5 * slices[i+1].v;
            add_slice.V = slices[i].V * 0.5;
            //add_slice.V = pi * (slices[i+1].x - slices[i].x) / 3. *(add_slice[i].)
            //std::cout << i << std::endl;
            slices.insert(slices.begin() + i + 1, add_slice);
            if (deltat > 1e-8)deltat *= 0.9;
            slices[i].V *= 0.5;
            
            slices[i].l = sqrt(pow(slices[i].x - slices[i+1].x, 2) + pow(slices[i].r - slices[i+1].r, 2));
            slices[i+1].l = sqrt(pow(slices[i+1].x - slices[i+2].x, 2) + pow(slices[i+1].r - slices[i+2].r, 2));
            
        } else {
            ++i; 
        }
    }//*/
    //boundary
    slices[0].r = R;
    
    slices[slices.size()-1].r = 0;

    Eigen::SparseMatrix<double> A_;
    Eigen::VectorXd B_ = Eigen::VectorXd::Zero(slices.size() - 2);
    std::vector<Eigen::Triplet<double>> tripletList;
    Eigen::VectorXd solution_;
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver_;
    for (int i = 2; i < slices.size() ; ++i){
        double pEp_px = 0, pr0_px = 0, prl_px = 0;
        pEp_px -=  rho * g * slices[i-1].V * h_c(slices[i].r, slices[i-1].r);
        pEp_px += 3 * eta * (slices[i].v - slices[i-1].v)/pow(slices[i].x - slices[i-1].x,2)*slices[i-1].V;
        if (i < slices.size() - 1){
            pEp_px -=  rho * g * slices[i].V * h_c(slices[i].r, slices[i+1].r);
            pEp_px += 3 * eta * (slices[i].v - slices[i+1].v)/pow(slices[i+1].x - slices[i].x,2)*slices[i].V;
        }
        //r_i-1 to i-1,i
        prl_px = -0.5 * sqrt(slices[i-1].V / pi / pow(slices[i].x - slices[i-1].x,3));
        pEp_px += sigma * pi *(slices[i-1].l + (pow(slices[i-1].r,2) - pow(slices[i].r,2))/slices[i-1].l)*prl_px
            +sigma * pi *(slices[i].r + slices[i-1].r) * (slices[i].x - slices[i-1].x)/slices[i-1].l;
        
        if (i < slices.size() - 1){
            pr0_px = 0.5 * sqrt(slices[i].V / pi / pow(slices[i+1].x - slices[i].x,3));
            //r_i to i,i+1
            pEp_px += sigma * pi *(slices[i].l + (pow(slices[i].r,2) - pow(slices[i+1].r,2))/slices[i].l)*pr0_px
                -sigma * pi *(slices[i+1].r + slices[i].r) * (slices[i+1].x - slices[i].x)/slices[i].l;
            //r_i to i-1,i
            pEp_px += sigma * pi *(slices[i-1].l + (pow(slices[i].r,2) - pow(slices[i-1].r,2))/slices[i-1].l)*pr0_px;
        }
        if (i>1){
            //r_i-1 to i-2,i-1
            pEp_px += sigma * pi *(slices[i-2].l + (pow(slices[i-1].r,2) - pow(slices[i-2].r,2))/slices[i-2].l)*prl_px;
        }
        double xi_tt = 0;
        //xi_tt = -1 / rho / slices[i-1].V * pEp_px;
        B_(i-2) = -pEp_px;
        tripletList.emplace_back(i-2, i-2,  rho * slices[i-1].V);
        /*
        if (i > 2){
            tripletList.emplace_back(i-2, i-2, 0.25* rho * slices[i-1].V);
            tripletList.emplace_back(i-2, i-3, 0.25* rho * slices[i-1].V);
            tripletList.emplace_back(i-2, i-2, 1.0/8.0* rho * pow(slices[i-1].V, 2) /pi /pow(slices[i].x-slices[i-1].x, 3));
            tripletList.emplace_back(i-2, i-3, -1.0/8.0* rho * pow(slices[i-1].V, 2) /pi /pow(slices[i].x-slices[i-1].x, 3));
            B_(i-2) += 3.0/16.0 *rho * pow(slices[i-1].V,2)/pi * pow(slices[i].v - slices[i-1].v,2) / pow(slices[i].x - slices[i-1].x, 4);
        }
        if (i < slices.size() - 1){
            tripletList.emplace_back(i-2, i-2, 0.25* rho * slices[i].V);
            tripletList.emplace_back(i-2, i-1, 0.25* rho * slices[i].V);
            tripletList.emplace_back(i-2, i-2, 1.0/8.0* rho * pow(slices[i].V, 2) /pi /pow(slices[i+1].x-slices[i].x, 3));
            tripletList.emplace_back(i-2, i-1, -1.0/8.0* rho * pow(slices[i].V, 2) /pi /pow(slices[i+1].x-slices[i].x, 3));
            B_(i-2) -= 3.0/16.0 *rho * pow(slices[i].V,2)/pi * pow(slices[i+1].v - slices[i].v,2) / pow(slices[i+1].x - slices[i].x, 4);
        }//*/
    }

    A_.resize(slices.size()-2, slices.size()-2);
    A_.setFromTriplets(tripletList.begin(), tripletList.end());
    //std::cout << A_ << std::endl;
    solver_.compute(A_);
    solution_ = solver_.solve(B_);
    if (solver_.info() != Eigen::Success) {
        std::cerr << "fail" << std::endl;
    }
    //std::cout << solution_ << std::endl;
    double E_k_x = 0,E_k_r = 0;
    for (int i = 0; i < slices.size() ; ++i){
        if (i >= 2)slices[i].v += solution_(i-2) * deltat;
        slices[i].x += slices[i].v * deltat ;
        E_k_x += 1.0/2.0 * pow(slices[i].v ,2) * rho * slices[i].V ;
        if (i < slices.size() - 1)E_k_r += 1.0/4.0 * rho * slices[i].V * pow((slices[i+1].v - slices[i].v)/2/(slices[i+1].x - slices[i].x) * slices[i].r,2) ;
        //std::cout << i << ':' <<slices[i].v << std::endl;
    }
    //std::cout << E_k_r / E_k_x << std::endl;
    
}



//namespace plt = matplotlibcpp;
int main(){
    std::cout << "begin" << std::endl;
    //slice_init();
    read();
    flow();
    
    sf::RenderWindow window(sf::VideoMode({800, 600}), "My window");
    sf::CircleShape circle(1);
    circle.setFillColor({255,255,255,255});
    
    // run the program as long as the window is open

    int t_ = 0;
    while (window.isOpen())
    {
        compute();
        ++t_;
        //if (t > 2000)break;
        while (const std::optional event = window.pollEvent())
        {
            // "close requested" event: we close the window
            if (event->is<sf::Event::Closed>())
                window.close();
        }

        if (t_%1000 == 0){
            window.clear(sf::Color::Black);
            sf::CircleShape shape(50.f);
            //std::cout << t << std::endl;
            shape.setFillColor(sf::Color(255, 255, 255));
            if (sf::Keyboard::isKeyPressed(sf::Keyboard::Key::Left))
            {
                // left key is pressed: move our character
                deltat *=1.01;
            }
            outFile << slices.size() << ' ' <<t << std::endl;
            for (int j = 0; j < slices.size(); ++j){
                circle.setPosition(sf::Vector2f(4e4 * slices[j].x, 300 - 4e4 * slices[j].r));
                window.draw(circle);
                circle.setPosition(sf::Vector2f(4e4 * slices[j].x, 300 + 4e4 * slices[j].r));
                outFile << slices[j].x << ' ' <<slices[j].r << std::endl;
                //circle.setPosition(sf::Vector2f(400, 600));
                //std::cout << 1e6 * slices[j].x << std::endl;
                window.draw(circle);
            }
            //circle.setPosition(sf::Vector2f(400, 200));
            window.draw(circle);
            window.display();
            
        }
        //_sleep(5000);
    }
    outFile.close();
    outFile.open(initFilePath);
    outFile << slices.size() << std::endl;
    for (int i = 0; i < slices.size(); ++i){
        outFile << slices[i].V << ' ' << slices[i].v << ' ' << slices[i].x << std::endl;
    }
    return 0;
}

/*
    Py_Initialize();
    PyRun_SimpleString("import matplotlib; matplotlib.use('Qt5Agg')");
    PyRun_SimpleString("import matplotlib; print('Current backend:', matplotlib.get_backend())");
    //*
    plt::ion(); // 启用交互模式
    //plt::figure();
    std::vector<double> x, y, x2, y2;
    plt::clf();
    for (int i = 0; i < 2000; ++i) {
        // 生成新数据点
        compute();
        //std::cout << i << std::endl;
        plt::clf();
        if (i%100 == 0){
            std::cout << i << std::endl;
            x.clear();
            y.clear();
            x2.clear();
            y2.clear();
            for (int j = 0; j < slices.size(); ++j){
                x.push_back(slices[j].x);
                y.push_back(slices[j].r);
            }
            plt::clf();

            
            plt::subplot(1, 2, 2);
            
            plt::plot(x2, y2, "r-");
            plt::title("Dynamic Plot 2");
            plt::xlabel("X");
            plt::ylabel("Value");
            ///
            // 短暂暂停以更新图表
            plt::pause(0.001);
            plt::show();
            std::this_thread::sleep_for(std::chrono::milliseconds(500));
        }
        // 清除当前图形
    }

    plt::show(); // 保持窗口显示
    //
    Py_Finalize();
    //*/
