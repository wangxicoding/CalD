#include <iostream>
#include <stdio.h>
#include <math.h>
#include <vector>
#include <include/CElement.h>
#include <include/CONTACT.h>

using namespace std;

const int maxn = 5;

/**
 * 材料参数
 */
const double fc = 14.3; //抗压强度
const double ft = 1.43; //抗拉强度
const double Ec = 3.0 * 1e4; //弹性模量
const double vc = 0.2 ;//泊松比
const double Gc = 1.2 * 1e4; //剪切模量
const double W = 25; //重度
const double rho = 2500; //密度

/**
 * 单元参数
 */
const double r = 0.2; //半径
const double delta = 0.01; //厚度
const double alpha = 0.28; //质量阻尼系数

const double C = 3.18 * 1e6; //粘性系数
const double miu = 0.6; //摩擦系数
const double beta = 0.5; //刚强阻尼系数

/**
 * 纵筋
 */
const double zEs = 2.00 * 1e5; //弹性模量
const double zfy = 300; //抗拉强度
const double zfyp = 300; //抗压强度 一撇
const double zA = 1017.876002; //钢筋界面

/**
 * 箍筋
 */
const double gEs = 2.1 * 1e5; //弹性模量
const double gfy = 210; //抗拉强度
const double gfyp = 210; //抗压强度 一撇

double Kn;
double Ks;
double Cp;
double Ft;
double Fc;
double m;
double I;

extern double slqh(double x, double y);

/**
 * 前处理——参数计算
 */
void init()
{
    /** 1 **/

    /** 2 **/
    Kn = (sqrt(3) * Ec * delta) / (3 * (1 - vc));
    Ks = Kn * (1 - 3 * vc);


    /** 3 **/

    /** 4 **/
    Cp = C * 2 * r * delta;

    /** 5 **/
    Ft = ft * 2 *r * delta;
    Fc = fc * 2 *r * delta;


    /** 6 **/
    /** 7 **/
    m = rho * 0.5 * M_PI * r * r * 0.01;

    /** 8 **/
    I = 0.5 * m * r * r;



}


vector<CElement> unit;

void unitInit()
{
    /** 初始化单元 **/
    for (int i = 0; i < maxn; i++)
    {
        CElement cElement;

        cElement.m_mass = m;
        cElement.m_number = i;
        cElement.m_moment = I;

        cElement.m_velx = 0;
        cElement.m_vely = 0;
        cElement.m_rotate = 0;

        cElement.m_r = r;
        cElement.m_x = 0;
        cElement.m_y = r + 2 * r * i;

        cElement.m_disx = 0;
        cElement.m_disy = 0;
        cElement.m_diso = 0;

        unit.push_back(cElement);
    }

    /** 初始化接触对象 **/
    for (int i = 0; i < maxn -1; i++)
    {
        CONTACT contact(&unit[i + 1]);
        CONTACT contact1(&unit[i]);

        unit[i].contactMap.insert(pair<int, CONTACT>(i + 1, contact));
        unit[i + 1].contactMap.insert(pair<int, CONTACT>(i, contact1));

        cout<< unit[i].m_mass <<endl;

    }

}

void calculate()
{
    for (int i = 0; i < maxn; i++)
    {
        CElement &iUnit = unit[i];
        for (int j = i + 1; j < maxn; j++)
        {
            CElement &jUnit = unit[j];
            map<int, CONTACT>::iterator itMap;

            itMap = iUnit.contactMap.find(j);

            double D;
            D = slqh(iUnit.m_x - jUnit.m_x, iUnit.m_y - jUnit.m_y);

            if (D > iUnit.m_r + jUnit.m_r) //如果大于
            {
                if (itMap == iUnit.contactMap.end()) //没弹簧、没有碰撞
                {
                    continue; //跳过
                }
                else // 有弹簧或者碰撞
                {
                    CONTACT &contact = itMap->second;
                    if (contact.isSpring == false) //不是弹簧,是碰撞,删除contact
                    {
                        iUnit.contactMap.erase(itMap);
                        jUnit.contactMap.erase(i);
                    }
                    else //是弹簧,拉力
                    {
                        //计算
                    }
                } //end else
            } //end if
            else //如果小于, 那一定会有力
            {
                if (itMap == iUnit.contactMap.end()) //没弹簧、没碰撞, 那么添加碰撞
                {
                    CONTACT contact(&unit[j], false); //添加碰撞
                    CONTACT contact1(&unit[i], false); //添加碰撞
                    iUnit.contactMap.insert(pair<int, CONTACT>(j, contact));
                    jUnit.contactMap.insert(pair<int, CONTACT>(i, contact1));

                    //计算碰撞力

                }
                else //有弹簧或者碰撞
                {
                    CONTACT &contact = itMap->second;
                    if (contact.isSpring == false) //不是弹簧,是碰撞
                    {
                        //计算碰撞力
                    }
                    else //是弹簧
                    {
                        //计算
                    }
                } //end else
            } //end else

        } //end for j
    } //end for i
}

int main()
{
    init();
    unitInit();
    calculate();
    cout << "Hello world!" << endl;
    return 0;
}


















