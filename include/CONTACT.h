#ifndef CONTACT_H
#define CONTACT_H

#include <include/CElement.h>

class CElement;

class CONTACT
{
    public:

        CElement *p_element; //接触单元的指针
        double m_fn; //接触法向力
        double m_fs; //接触切向力
        double m_dn;//法向阻尼力Dn
        double m_ds;//法向阻尼力Ds
        bool fenli; //判断单元间切向弹簧是否被破坏
        bool fencon; //判断单元间法向弹簧是否被破坏
        bool isSpring; //是否是弹簧

        double m_Fx;
        double m_Fy;
        double m_M;

        CONTACT(CElement *cElement, bool spring = true);

    protected:
    private:
};

#endif // CONTACT_H
