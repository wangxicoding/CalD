#include "include/CElement.h"
#include <math.h>


double CElement::Kn;
double CElement::Ks;
double CElement::Cp;
double CElement::Ft;
double CElement::Fc;
double CElement::m;
double CElement::I;


/**
 * @brief slqh 矢量求和函数
 * @param x
 * @param y
 * @return
 */
double slqh(double x, double y)
{
    double xy;
    xy=sqrt(x*x+y*y);
    return xy;
}

/**
 * @brief sl_angel 根据向量求角度，x，y分别表示XY方向的分量
 * @param x
 * @param y
 * @return
 */
double sl_angel(double x, double y)
{
    double angel;
    if(x==0)
    {
        if(y>0)
        {
            angel=0.5*M_PI;//Y轴正方向
        }
        else
        {
            angel=1.5*M_PI;//Y轴负方向
        }
    }    //过滤分母为零的情况，以免出现不必要的错误
    else
    {
        if(x>0)
        {
            if(y>0)
            {
                angel=fabs(atan(y/x));//第一象限角
            }
            else
            {
                angel=2*M_PI-fabs(atan(fabs(y)/x));//第四象限角
            }
        }
        else
        {
            if(y>0)
            {
                angel=M_PI-fabs(atan(y/fabs(x)));//第二象限角
            }
            else
            {
                angel=M_PI+fabs(atan(fabs(y)/fabs(x)));//第三象限角
            }

        }

    }
    if(fabs(angel-2*M_PI)<0.00001)//过滤360度情况
        {
            angel=0;
        }
    return angel;//返回向量角度，单位为弧度制
}


/**
 * 前处理——参数计算
 */
void CElement::staticInit()
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

double CElement::AnsysElement()
{
    double m_area=0.5*M_PI*m_r*m_r;
    m_moment=0.5*m_mass*m_r*m_r;
    return m_area;
}

double CElement::Initial()
{
    //单元性质
    m_number=0; //单元编号初始化
    m_x=m_y=m_r=0.0; //单元x方向的坐标，单元y方向的坐标，单元的半径初始化
    m_velx=m_vely=m_rotate=0.0; //单元x，y方向的速度及单元转动速度初始化
    vt_t2[4]=vtt2[4]=mtt2=mt_t2=0.0; //vt_t2和mt_t2分别表示t-△t/2时刻的速度和转动速度，vtt2和mtt2分别表示t+△t/2时刻的速度和转动速度初始化
    m_mass=m_moment=0.0; //单元的质量和转动惯量初始化
    m_disx=m_disy=m_diso=0.0; //分别表示单元x，y方向的位移增量及转角增量初始化
    m_disn=m_diss=0.0; //分别表示单元法向方向位移增量，切向方向位移增量初始化
    m_ag=0.0; //地震加速度记录，从地震数据中读取初始化
    m_Feqx=m_Feqy=0.0; //x方向和y方向地震力初始化
    m_Fxsum=m_Fysum=m_Msum=0.0; //x和y方向和合力及合力矩初始化
    //材料性质参数
    m_young=0.0;
    m_possion=0.0;
    m_thickness=0.0;
    m_weight=0.0;

    return 0;
}


/**
 * @brief CElement::calForce 计算接触力，p2单元作用到p1单元上的力，
 * 同时也储存到p2的接触力容器中，便于计算p2的接触力
 * @param p1
 * @param p2
 */
void CElement::calContactForce(CElement* p2, CONTACT* cont1)//计算接触力有弹簧，p2单元作用到p1单元上的力，同时也储存到p2的接触力容器中，便于计算p2的接触力
{
   double Kn;
double Ks;
double Cp;

    double l_xy; //代表两颗粒之间绝对距离
    double l_xc; //l_xc表示接触点到p1单元圆心的距离
    double lx, ly, c_x, c_y; //
    double sl_u[3]; //定义接触力矢量数组
    double sumR; //两个单元半径和
    double cos_theta,sin_theta; //定义两个单元的圆心连线与x轴正方向的余弦值和正弦值
    double delta_n,delta_s; //法向位移增量，切向位移增量
    double b1, b2;//法向力增量，切向力增量，法向刚度增量，切向刚度增量

    // 15963 下面几个我定义的
    //double beta = 0;
    //double deltaTime = 0;
    //double C = 0;
    //double mu = 0;
    //double Kn = 0;
    //double Ks = 0;

    sl_u[0]=(p2->m_x) - (this->m_x);
    sl_u[1]=(p2->m_y) - (this->m_y);

    l_xy=slqh(sl_u[0], sl_u[1]);//调用求和程序，给定两个单元的x，y坐标值，计算两个单元的绝对距离
    cos_theta=(p2->m_x - this->m_x)/l_xy;//计算两单元型心连线与x轴的夹角余弦值
    sin_theta=(p2->m_y - this->m_y)/l_xy;//计算两单元型心连线与x轴的夹角正弦值

    l_xc=0.5*(this->m_r - p2->m_r + l_xy); //

    lx=l_xc * cos_theta; // 接触点到 p1 x坐标距离
    ly=l_xc * sin_theta; // 接触点到 p1 y坐标距离

    c_x = lx + this->m_x; //接触点c的x坐标
    c_y = ly + this->m_y; //接触点c的y坐标

    b1 = p2->m_disx - this->m_disx; //
    b2 = p2->m_disy - this->m_disy; //
    //对应第2页第4步
    this->m_disn = b1*cos_theta + b2*sin_theta;//计算法向方向位移增量
    this->m_diss = -b1*sin_theta + b2*cos_theta - ((this->m_r)*(this->m_diso) + (p2->m_r)*(p2->m_diso));//切向方向位移增量
    //对应第2页第5步
    cont1->m_fn = cont1->m_fn - Kn*this->m_disn;//求法向力fn
    cont1->m_fs = cont1->m_fs + Ks*this->m_diss;//求切向力fs
    //对应第2页第6步
    cont1->m_dn=-beta*Kn*(this->m_disn);//求法向阻尼力Dn，beta为参数，设置为0.6
    cont1->m_ds=-beta*Ks*(this->m_diss);//求切向阻尼力Ds
    //对应第3页第7步
    this->tau=Cp+miu*(cont1->m_fn);//求极限剪力值tau，C为常数,这个tau用于后面的判断，在DiscreteElement.cpp中



    //求合力及合力矩
    this->m_Fxsum = -(cont1->m_fs + cont1->m_ds)*sin_theta - (cont1->m_fn + cont1->m_dn)*cos_theta;
    this->m_Fysum = (cont1->m_fs + cont1->m_ds)*cos_theta-(cont1->m_fn + cont1->m_dn)*sin_theta;
    this->m_Msum = this->m_Fysum * (c_x - this->m_x) - this->m_Fxsum * (c_y - this->m_y);
    //作用到p2单元上的力为
    p2->m_Fxsum=-this->m_Fxsum;
    p2->m_Fysum=-this->m_Fysum;
    p2->m_Msum=-this->m_Msum;
    // 15963 所有的力都放在contact里，放的地方不一致，用的时候怎么办， 这地方你看一下我这样对不对
    cont1->m_Fx=this->m_Fxsum;
    cont1->m_Fy=this->m_Fysum;
    cont1->m_M=this->m_Msum;

}
void CElement::calCollisionForce(CElement* p2, CONTACT* cont1) //结算碰撞力，无弹簧 15963 p2 是另一个单元 cont1 是他们两个的关系 原来的p1现在是当前单元this
{double Kn;
double Ks;
double Cp;

    int i; //循环变量
    double l_xy; //代表两颗粒之间绝对距离
    double sl_u[3]; //定义接触力矢量数组
    double cos_theta,sin_theta; //定义两个单元的圆心连线与x轴正方向的余弦值和正弦值
//double delta_n,delta_s; //法向位移增量，切向位移增量
    double sumR; //两个单元半径和
    double e[2],t[2];
    double f_n, f_s, d_n, d_s;//法向力增量，切向力增量，法向刚度增量，切向刚度增量
    double delta_Fx, delta_Fy;

    sumR = this->m_r + p2->m_r; //求两个单元半径和

    sl_u[0]=(p2->m_x)-(this->m_x);
    sl_u[1]=(p2->m_y)-(this->m_y);

    l_xy=slqh(sl_u[0],sl_u[1]);//调用求和程序，给定两个单元的x，y坐标值，计算两个单元的绝对距离
    cos_theta=(p2->m_x-this->m_x)/l_xy;//计算两单元型心连线与x轴的夹角余弦值
    sin_theta=(p2->m_y-this->m_y)/l_xy;//计算两单元型心连线与x轴的夹角正弦值


    if(l_xy>sumR) // 15963 这个什么意思，是到这就不算了吧
    {
        return;
    }
    else
    {
        e[0]=((p2->m_x)-(this->m_x))/l_xy; //cos_theta
        e[1]=((p2->m_y)-(this->m_y))/l_xy; //sin_theta
        t[0]=e[1]; //sin_theta
        t[1]=-e[0]; //-cos_theta
        // 15963 下面几个我定义的
        double vn_XD[10];
        double vs_XD[10];
        double drt_n[10];
        double drt_s[10];


        for(i=0;i<2;i++)
        {
            vn_XD[0]=((this->m_velx)-(p2->m_velx))*e[0];//单元相对速度x方向分量
            vn_XD[1]=((this->m_vely)-(p2->m_vely))*e[1];//单元相对速度y方向分量
            vs_XD[0]=(((this->m_velx)-(p2->m_velx))*t[0])-(this->m_rotate*this->m_r+p2->m_rotate*p2->m_r);//相对速度
            vs_XD[1]=(((this->m_vely)-(p2->m_vely))*t[1])-(this->m_rotate*this->m_r+p2->m_rotate*p2->m_r);//相对速度
            drt_n[i]=vn_XD[i]*deltaTime;//相对位移增量
            drt_s[i]=vs_XD[i]*deltaTime;//相对位移增量
        }
            vn_XD[3]=slqh(vn_XD[0],vn_XD[1]);//法向相对速度大小
            vs_XD[3]=slqh(vs_XD[0],vs_XD[1]);//切向相对速度大小
            drt_n[3]=slqh(drt_n[0],drt_n[1]);//法向的位移增量
            drt_s[3]=slqh(drt_s[0],drt_s[1]);//切向位移增量值

            f_n=Kn*drt_n[3];//法向力增量
            f_s=Ks*drt_s[3];//切向力增量
            d_n=-beta*Kn*vn_XD[3];//法向刚度增量
            d_s=-beta*Ks*vs_XD[3];//切向刚度增量

            cont1->m_fn=cont1->m_fn+f_n+d_n;
            cont1->m_fs=cont1->m_fs+f_s+d_s;
            this->tau=(cont1->m_fn)*miu+Cp;//库伦摩擦力

            if(abs(cont1->m_fs) > abs(this->tau))
            {
                cont1->m_fs=abs(this->tau)*(cont1->m_fs/abs(cont1->m_fs));//方向与原摩阻力方向一致
            }
            else
            {
                delta_Fx=cont1->m_fn*e[0]+cont1->m_fs*e[1];
                delta_Fy=cont1->m_fn*e[1]-cont1->m_fs*e[0];
                cont1->m_Fx=cont1->m_Fx+delta_Fx;//求x方向合力，未考虑外力和重力
                cont1->m_Fy=cont1->m_Fy+delta_Fy;//求x方向合力，未考虑外力和重力
                cont1->m_M=cont1->m_M+cont1->m_fs*this->m_r;//求合力矩
            }
    }
}



void CElement::cal_vtt2() //这个单元就是他自己，全部使用this
{
    double c1, c2;

    c1=1-alpha*deltaTime/2;
    c2=1+alpha*deltaTime/2;
    this->vtt2[0]=(this->vt_t2[0] * c1 + this->m_Fxsum * deltaTime / this->m_mass) / c2;
    this->vtt2[1]=(this->vt_t2[1]*c1+(this->m_Fysum/this->m_mass+g)*deltaTime)/c2;
    this->mtt2=(this->mt_t2*c1+(this->m_Msum/this->moment)*deltaTime)/c2;
    this->vtt2[2]=sl_angel(this->vtt2[0],this->vtt2[1]);
    this->vtt2[3]=slqh(this->vtt2[0],this->vtt2[1]);
}

void CElement::cal_vt()
{
    this->m_velx=(this->vtt2[0]+this->vt_t2[0])/2;
    this->m_vely=(this->vtt2[1]+this->vt_t2[1])/2;
    this->m_rotate=(this->mtt2+this->mt_t2)/2;
}

void CElement::cal_dis()
{

    this->m_disx=this->vtt2[0]*deltaTime;
    this->m_disy=this->vtt2[1]*deltaTime;
    this->m_diso=this->mtt2*deltaTime;
}

void CElement::cal_utt()
{
    this->m_x = this->m_x + this->m_disx;
    this->m_y = this->m_y + this->m_disy;
    this->m_o = this->m_o + this->m_diso;
}

void CElement::change_of_data()
{
    int i;
    for(i=0; i<4; i++)
    {
        this->vt_t2[i]=this->vtt2[i];
    }
}

void CElement::union_lisan()
{
    cal_vtt2();
    cal_vt();
    cal_dis();
    cal_utt();
    change_of_data();
}
