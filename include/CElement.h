#ifndef CELEMENT_H
#define CELEMENT_H

#include <list>
#include <include/CONTACT.h>

const double g = 9.80665;
const double dt = 1e-6;

class CONTACT;

class CElement
{
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
        std::list<CONTACT> contactList; //接触链表

        double vt_t2[4], vtt2[4], mtt2, mt_t2; //vt_t2和mt_t2分别表示t-△t/2时刻的速度和转动速度，vtt2和mtt2分别表示t+△t/2时刻的速度和转动速度
        double m_disn, m_diss; //分别表示单元法向方向位移增量，切向方向位移增量
        double m_ag; //地震加速度记录，从地震数据中读取
        double m_Feqx, m_Feqy; //x方向和y方向地震力
        double m_Fxsum, m_Fysum, m_Msum; //x和y方向和合力及合力矩

        //材料参数
        double m_young;//弹性模量
        double m_possion;//泊松比
        double m_thickness;//单元厚度
        double m_weight;//单元容重

    public:
        double AnsysElement();
        double Initial();
        void calForce(CElement *p1, CElement *p2);
        void cal_vtt2(CElement *p);





};

#endif // CELEMENT_H
