#define MATPLOTLIBCPP_ENABLE_QT5
#include <iostream>
#include <algorithm>
#include <cmath>
#include <vector>
#include <chrono>
#include <thread>
#include <fstream>
#include <string>
#include <ctime>
#include <memory>
#include <SFML/Graphics.hpp>
#include <SFML/Audio.hpp>

const double R = 1.6e-3;
const double pi = 3.1416;
const double g = 9.8;
double slower = 10.0;
double delta_t = 0.002;
class slice{
  public:
    double V;
    double r;
    double v;
    double x;
    double l;
};
std::vector<slice> slices;
struct st{
    std::vector<slice> slices_;
    double t;
};
std::shared_ptr<std::vector<st>> st_1 = nullptr,st_2 = nullptr;
//namespace plt = matplotlibcpp;
bool init(std::string name, std::shared_ptr<std::vector<st>> sl_t){
    //std::string dataFilePath = std::string(DATA_DIR) + "/example.txt";
    std::string dataFilePath = std::string(DATA_DIR) + "/" + name;
    std::ifstream file(dataFilePath);
    if (!file.is_open()) {
        std::cerr << "无法打开文件: " << dataFilePath << std::endl;
        return 1;
    }
    std::string line;
    if (!std::getline(file, line))  return 1;
    int n;
    while (file >> n) {
        double t ;
        file >> t;
        //std::cout << t << std::endl;
        slices.clear();
        slice c_slice;
        st c_st;
        for (int i = 0; i < n; ++i){
            file >> c_slice.x;
            file >> c_slice.r;
            slices.push_back(c_slice);
        }
        c_st.slices_ = slices;
        c_st.t = t;
        sl_t->push_back(c_st);
        //std::cout << sl_t->size() << std::endl;
    }
    file.close();
    return 0;
}
int main(){
    std::cout << "begin" << std::endl;
    st_1 = std::make_shared<std::vector<st>>();
    st_2 = std::make_shared<std::vector<st>>();
    init("16cms2.txt",st_1);
    init("workspace.txt",st_2);
    std::cout << st_1->size() << std::endl;
    sf::RenderWindow window(sf::VideoMode({800, 600}), "My window");
    sf::CircleShape circle(2);
    circle.setFillColor({0,0,255,255});
    sf::CircleShape circle_2(2);
    circle_2.setFillColor({255,0,0,255});
    // run the program as long as the window is open

    int t = 0;
    double last_t = 0;
    auto it_1 = st_1->begin(), it_2 = st_2->begin();
    std::cout << st_1->size() << std::endl;
    std::cout << st_2->size() << std::endl;
    while (window.isOpen())
    {
        last_t += delta_t;
        //std::cout << last_t << std::endl;
        while (it_1->t < last_t  && it_1 != st_1->end()) {
            ++it_1;
            //std::cout << it_1->t << std::endl;
        }
        while (it_2->t < last_t  && it_2 != st_2->end()) ++it_2;
        if (it_1 == st_1->end() || it_2 == st_2->end()){
            it_1 = st_1->begin(), it_2 = st_2->begin();
            last_t = 0;
        }
        std::cout << last_t << std::endl;
        //if (t > 2000)break;
        while (const std::optional event = window.pollEvent())
        {
            // "close requested" event: we close the window
            if (event->is<sf::Event::Closed>())
                window.close();
        }
        //std::cout << it->size();
        
        window.clear(sf::Color::White);
        sf::CircleShape shape(50.f);
        //std::cout << t << std::endl;
        shape.setFillColor(sf::Color(255, 255, 255));
        for (int j = 0; j < it_1->slices_.size(); ++j){
            circle.setPosition(sf::Vector2f(4e4 * (it_1->slices_)[j].x, 300 - 4e4 * (it_1->slices_)[j].r));
            window.draw(circle);
            circle.setPosition(sf::Vector2f(4e4 * (it_1->slices_)[j].x, 300 + 4e4 * (it_1->slices_)[j].r));
            //circle.setPosition(sf::Vector2f(400, 600));
            //std::cout << 1e6 * slices[j].x << std::endl;
            window.draw(circle);
        }
        for (int j = 0; j < it_2->slices_.size(); ++j){
            circle_2.setPosition(sf::Vector2f(4e4 * (it_2->slices_)[j].x, 300 - 4e4 * (it_2->slices_)[j].r));
            window.draw(circle_2);
            circle_2.setPosition(sf::Vector2f(4e4 * (it_2->slices_)[j].x, 300 + 4e4 * (it_2->slices_)[j].r));
            //circle.setPosition(sf::Vector2f(400, 600));
            //std::cout << 1e6 * slices[j].x << std::endl;
            window.draw(circle_2);
        }
        //circle.setPosition(sf::Vector2f(400, 200));
        window.draw(circle);
        window.display();
        std::this_thread::sleep_for(std::chrono::milliseconds(int(delta_t * 1e3 * slower)));
        //_sleep(5000);
    }
    return 0;
}


