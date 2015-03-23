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
        CElement cElement(i, 0, r + 2*r*i);

        unit.push_back(cElement);
    }

    /** 初始化接触对象 **/
    for (int i = 0; i < maxn -1; i++)
    {
        CONTACT contact(&unit[i + 1]);
        CONTACT contact1(&unit[i]);

        /** 接触,成对的出现  ╭(●｀∀´●)╯╰(●'◡'●)╮ 好基友,一辈子 **/
//        contact.p_partner = &contact1;
//        contact1.p_partner = &contact;

        unit[i].contactMap.insert(pair<int, CONTACT>(i + 1, contact));
        unit[i + 1].contactMap.insert(pair<int, CONTACT>(i, contact1));

        // 因为c++构造函数的问题, 上面那个会出现指针错误
        (unit[i].contactMap[i + 1]).p_partner = &(unit[i + 1].contactMap[i]);
        (unit[i + 1].contactMap[i]).p_partner = &(unit[i].contactMap[i + 1]);

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

                        continue; //跳过
                    }
                    else //是弹簧,拉力
                    {
                        //计算弹簧力
                        iUnit.calContactForce(&jUnit, &contact);

                    }
                } //end else
            } //end if
            else //如果小于, 那一定会有力
            {
                if (itMap == iUnit.contactMap.end()) //没弹簧、没碰撞, 那么添加碰撞
                {
                    CONTACT contact(&unit[j], false); //添加碰撞
                    CONTACT contact1(&unit[i], false); //添加碰撞

                    /** 接触,成对的出现  ╭(●｀∀´●)╯╰(●'◡'●)╮ 好基友,一辈子 **/
//                    contact.p_partner = &contact1;
//                    contact1.p_partner = &contact;

                    iUnit.contactMap.insert(pair<int, CONTACT>(j, contact));
                    jUnit.contactMap.insert(pair<int, CONTACT>(i, contact1));

                    // 因为c++构造函数的问题, 上面那个会出现指针错误
                    iUnit.contactMap[j].p_partner = &(jUnit.contactMap[i]);
                    jUnit.contactMap[i].p_partner = &(iUnit.contactMap[j]);


                    //计算碰撞力
                    iUnit.calCollisionForce(&jUnit, &contact);

                }
                else //有弹簧或者碰撞
                {
                    CONTACT &contact = itMap->second;
                    if (contact.isSpring == false) //不是弹簧,是碰撞
                    {
                        //计算碰撞力
                        iUnit.calCollisionForce(&jUnit, &contact);
                    }
                    else //是弹簧
                    {
                        //计算弹簧力
                        iUnit.calContactForce(&jUnit, &contact);
                    }
                } //end else
            } //end else

        } //end for j

        // 计算完与所有单元接触的力
        /** 十二步以后的事儿 **/
        iUnit.union_lisan();

    } //end for i
}

int main()
{
    int maxStep = 1e2; /**总步长,每一步长相当于每一时间间隔
                         *如果deltaTime = 1e-6,步长为 1e6
                         *那总时间为1s.
                         */

    CElement::staticInit();
    unitInit();
    for (int i = 0; i < maxStep; i++) //计算相应步长的单元离散过程
    {
        calculate();
    }
    cout << "Hello world!" << endl;
    return 0;
}


















