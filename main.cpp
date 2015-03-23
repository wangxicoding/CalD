#include <iostream>
#include <stdio.h>
#include <math.h>
#include <vector>
#include <include/CElement.h>
#include <include/CONTACT.h>

using namespace std;

const int maxn = 5;

extern double slqh(double x, double y);

vector<CElement> unit;

void unitInit()
{
    int r = CElement::r;
    /** 初始化单元 **/
    for (int i = 0; i < maxn; i++)
    {
        CElement cElement;

        cElement.m_mass = CElement::m;
        cElement.m_number = i;
        cElement.m_moment = CElement::I;

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
    CElement::staticInit();
    unitInit();
    calculate();
    cout << "Hello world!" << endl;
    return 0;
}


















