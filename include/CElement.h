#ifndef CELEMENT_H
#define CELEMENT_H

#include <map>
#include <include/CONTACT.h>



class CONTACT;

class CElement
{
public:
    static double Kn;
    static double Ks;
    static double Cp;
    static double Ft;
    static double Fc;
    static double m;
    static double I;


//main函数用到,先设置成public
public:
    static const double r = 0.2; //半径


private:
    static const double g = 9.80665;
    static const double deltaTime = 1e-6;

    /**
     * 材料参数
     */
    static const double fc = 14.3; //抗压强度
    static const double ft = 1.43; //抗拉强度
    static const double Ec = 3.0 * 1e4; //弹性模量
    static const double vc = 0.2 ;//泊松比
    static const double Gc = 1.2 * 1e4; //剪切模量
    static const double W = 25; //重度
    static const double rho = 2500; //密度

    /**
     * 单元参数
     */
    static const double delta = 0.01; //厚度
    static const double alpha = 0.28; //质量阻尼系数

    static const double C = 3.18 * 1e6; //粘性系数
    static const double miu = 0.6; //摩擦系数
    static const double beta = 0.5; //刚强阻尼系数

    /**
     * 纵筋
     */
    static const double zEs = 2.00 * 1e5; //弹性模量
    static const double zfy = 300; //抗拉强度
    static const double zfyp = 300; //抗压强度 一撇
    static const double zA = 1017.876002; //钢筋界面

    /**
     * 箍筋
     */
    static const double gEs = 2.1 * 1e5; //弹性模量
    static const double gfy = 210; //抗拉强度
    static const double gfyp = 210; //抗压强度 一撇

    static bool IsBreak(double nowFn, double nowFs, double tau);
public:
    static void staticInit();


public:
    //        CElement();
    //        virtual ~CElement();
    double m_mass; //单元的质量
    int m_number; //单元的编号
    double m_moment; //单元的转动惯量
    double m_velx, m_vely, m_rotate; //x、y平动速度, 单元转动速度
    double m_x, m_y, m_r; //圆心x坐标, y坐标, 半径
    bool m_restx, m_resty, m_resto; //x、y、转动方向是否受约束
    double m_loadx, m_loady, m_loado; //x、y所受合力, 单元所受合力矩
    double m_disx, m_disy, m_diso; //x、y、转角方向的位移增量
    std::map<int, CONTACT> contactMap; //接触

    double vt_t2[4], vtt2[4], mtt2, mt_t2; //vt_t2和mt_t2分别表示t-△t/2时刻的速度和转动速度，vtt2和mtt2分别表示t+△t/2时刻的速度和转动速度
    double m_disn, m_diss; //分别表示单元法向方向位移增量，切向方向位移增量
    double tau; //极限切向力
    double m_ag; //地震加速度记录，从地震数据中读取
    double m_Feqx, m_Feqy; //x方向和y方向地震力
    double m_Fxsum, m_Fysum, m_Msum; //x和y方向和合力及合力矩

    //材料参数
    double m_young;//弹性模量
    double m_possion;//泊松比
    double m_thickness;//单元厚度
    double m_weight;//单元容重

    double moment;
    double m_o;

public:
    CElement(int number, int coordX, int coordY);
    double AnsysElement();
//    void Initial(int number, int coordX, int coordY);
    void calForce(CElement *p1, CElement *p2);
    void union_lisan();
    void calContactForce(CElement* p2, CONTACT* cont1);
    void calCollisionForce(CElement* p2, CONTACT* cont1);

private:
    void addSumF(CElement *p2, CONTACT *cont1);
    void cal_vtt2();
    void cal_vt();
    void cal_dis();
    void cal_utt();
    void change_of_data();
    void calSumF();
};

#endif // CELEMENT_H
