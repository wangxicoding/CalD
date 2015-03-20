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
        bool fenli; //判断单元间切向弹簧是否被破坏
        bool fencon; //判断单元间法向弹簧是否被破坏

    protected:
    private:
};

#endif // CONTACT_H
