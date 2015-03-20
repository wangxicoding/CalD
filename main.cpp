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
    for (int i = 0; i < maxn; i++)
    {
        CONTACT contact;
        CONTACT contact1;

        if( i != 0)
        {
            contact.p_element = &unit[i - 1];
            contact.m_fn = 0;
            contact.m_fs = 0;
            contact.fenli = false;
            contact.fencon = false;
            unit[i].contactList.push_back(contact);
        }

        if( i != maxn - 1 )
        {
            contact1.p_element = &unit[i + 1];
            contact1.m_fn = 0;
            contact1.m_fs = 0;
            contact1.fenli = false;
            contact1.fencon = false;
            unit[i].contactList.push_back(contact1);
        }

        cout<< unit[i].m_mass <<endl;

    }

}

double dis(double ui, double vi, double uj, double vj)
{
    double ans = (ui - uj) * (ui - uj) + (vi - vj) * (vi - vj);
    return sqrt(ans);
}

void calculate()
{
    for (int i = 0; i < maxn; i++)
    {
        CElement *iUnit = &unit[i];
        CElement *jUnit;

        list<CONTACT>::iterator spring; //连接的弹簧
        for (spring = iUnit->contactList.begin(); spring != iUnit->contactList.end(); ++spring)
        {
            jUnit = spring->p_element;

            /** 第二步 **/
            if (spring->fenli | spring->fencon) //如果有一个被破坏了
            {
                /**
                *跳到 第D步
                */
                continue;
            }
            else;

            /** 第三步 **/
            double ui,vi,uj,vj;

            ui = iUnit->m_x;
            vi = iUnit->m_y;

            uj = jUnit->m_x;
            vj = jUnit->m_y;

            double D, cosTheta, sinTheta;

            D = dis(ui, vi, uj, vj);
            cosTheta = (uj - ui) / D;
            sinTheta = sqrt(1 - cosTheta * cosTheta);

            /** 第四步 **/
            double deltaUi, deltaVi, deltaPi, deltaUj, deltaVj, deltaPj;
            double deltaUn, deltaUs;
            double ri, rj;

            ri = iUnit->m_r;
            rj = jUnit->m_r;

            deltaUi = iUnit->m_disx, deltaVi = iUnit->m_disy, deltaPi = iUnit->m_diso;
            deltaUj = jUnit->m_disx, deltaVj = jUnit->m_disy, deltaPj = jUnit->m_diso;

            deltaUn = (deltaUj - deltaUi) * cosTheta + (deltaVj - deltaVi) * sinTheta;
            deltaUs = -(deltaUj - deltaUi) * sinTheta + (deltaVj - deltaVi) * cosTheta - (ri * deltaPi + rj * deltaPj);

            /** 第五步 **/
            double fnc0, fsc0;
            double fnc, fsc;


            fnc0 = spring->m_fn, fsc0 = spring->m_fs;

            fnc = fnc0 - Kn * deltaUn;
            fsc = fsc0 + Ks * deltaUs;

            /** 第六步 **/
            double Dn, Ds;

            Dn = -beta * Kn * deltaUn;
            Ds = -beta * Ks * deltaUs;

            /** 第七步 **/
            double tao;

            tao = Cp + miu * fnc;

            /** 第八步 **/
            if (fnc > 0)
            {
                if (fnc > Fc)
                {
                    /**弹簧断裂,进入阶段2**/
                }
                else
                {
                    /**弹簧未断裂,仍是阶段1**/
                }

            }
            else
            {
                if (-fnc > Ft)
                {
                    /**弹簧断裂,进入阶段2**/
                }
                else
                {
                    /**弹簧未断裂,仍是阶段1**/
                }
            }

            if (fsc > tao)
            {
                /**弹簧断裂,进入阶段2,切向力 u*fnc **/
            }
            else
            {
                /**弹簧未断裂,仍是阶段1**/
            }

/******* 阶段1: ********/

            /** 第九步 **/
            double Fx, Fy;

            Fx = -(fsc + Ds) * sinTheta - (fnc + Dn) * cosTheta;
            Fy = (fsc + Ds) * cosTheta - (fnc + Dn) * sinTheta;

            /** 第十步 **/
            double Fxload, ag;
            ag = 0;  //这里到时候读取文件
            Fxload = -m * ag;

            /** 第十一步 **/




        }
    }
}

int main()
{
    init();
    unitInit();
    calculate();
    cout << "Hello world!" << endl;
    return 0;
}


















