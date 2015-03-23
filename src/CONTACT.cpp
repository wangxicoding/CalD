#include "include/CONTACT.h"

CONTACT::CONTACT(CElement *cElement, bool spring)
{
     p_element = cElement;
     m_fn = 0; //接触法向力
     m_fs = 0; //接触切向力
     m_dn = 0;//法向阻尼力Dn
     m_ds = 0;//法向阻尼力Ds
     fenli = false; //切向弹簧没有被破坏
     fencon = false; //法向弹簧没有被破坏
     isSpring = spring; //是弹簧

     m_Fx = 0;
     m_Fy = 0;
     m_M = 0;
}

//
//CONTACT::~CONTACT()
//{
//    //dtor
//}
