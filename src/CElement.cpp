#include "include/CElement.h"
#include <math.h>

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

double distance(double ui, double vi, double uj, double vj)
{
    double ans = (ui - uj) * (ui - uj) + (vi - vj) * (vi - vj);
    return sqrt(ans);
}

double distance(double deltaU, double deltaV)
{
    double ans = deltaU * deltaU + deltaV * deltaV;
    return sqrt(ans);
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
}


/**
 * @brief CElement::calForce 计算接触力，p2单元作用到p1单元上的力，
 * 同时也储存到p2的接触力容器中，便于计算p2的接触力
 * @param p1
 * @param p2
 */
//void CElement::calForce(CElement* iUnit, CElement* jUnit)//
//{
//    /** 第三步 **/
//    double ui,vi,uj,vj;

//    ui = iUnit->m_circle.m_x;
//    vi = iUnit->m_circle.m_y;

//    uj = jUnit->m_circle.m_x;
//    vj = jUnit->m_circle.m_y;

//    double dis; //代表两颗粒之间绝对距离
//    double cosTheta,sinTheta; //定义两个单元的圆心连线与x轴正方向的余弦值和正弦值

//    dis = distance(ui, vi, uj, vj);
//    cosTheta = (uj - ui) / dis;
//    sinTheta = sqrt(1 - cosTheta * cosTheta);

//    bool mark_of_contract = true;
//    int i; //循环变量
//    double sl_u[3]; //定义接触力矢量数组
//    double delta_n,delta_s; //法向位移增量，切向位移增量
//    double sumR;
//    double b1, b2, b3, b4;//法向力增量，切向力增量，法向刚度增量，切向刚度增量
//    double delta_Fx, delta_Fy;

//    sl_u[0]=(jUnit->m_x)-(iUnit->m_x);
//    sl_u[1]=(jUnit->m_y)-(iUnit->m_y);
//    if(mark_of_contract==1)//弹簧未断开的计算
//    {
//        dis=slqh(sl_u[0],sl_u[1]);//调用求和程序，给定两个单元的x，y坐标值，计算两个单元的绝对距离
//        cosTheta=(jUnit->m_x-iUnit->m_x)/dis;//计算两单元型心连线与x轴的夹角余弦值
//        sinTheta=(jUnit->m_y-iUnit->m_y)/dis;//正弦值

//        b1=jUnit->m_disx-iUnit->m_disx;
//        b2=jUnit->m_disy-iUnit->m_disy;
//        //对应第2页第4步
//        iUnit->m_disn=b1*cosTheta+b2*sinTheta;//计算法向方向位移增量
//        iUnit->m_diss=-b1*sinTheta+b2*cosTheta-((iUnit->m_r)*(iUnit->m_diso)+(jUnit->m_r)*(jUnit->m_diso));//切向方向位移增量
//        //对应第2页第5步
//        iUnit->p_contract->m_fn=iUnit->p_contract->m_fn-kn*iUnit->m_disn;//求法向力fn
//        iUnit->p_contract->m_fs=iUnit->p_contract->m_fn+ks*iUnit->m_diss;//求切向力fs
//        //对应第2页第6步
//        iUnit->p_contract->m_dn=-beta*kn*(iUnit->m_disn);//求法向阻尼力Dn，beta为参数，设置为0.6
//        iUnit->p_contract->m_ds=-beta*ks*(iUnit->m_diss);//求切向阻尼力Ds
//        //对应第3页第7步
//        iUnit->tao=C+0.6*(iUnit->p_contract->m_fn);//求极限剪力值tao，C为常数,这个tao用于后面的判断，在DiscreteElement.cpp中
//        //求合力
//        iUnit->p_contract->m_Fx=-(iUnit->p_contract->m_fs+iUnit->p_contract->m_ds)*sinST-(iUnit->p_contract->m_fn+iUnit->p_contract->m_dn)*cosST;
//        iUnit->p_contract->m_Fy=(iUnit->p_contract->m_fs+iUnit->p_contract->m_ds)*cosST-(iUnit->p_contract->m_fn+iUnit->p_contract->m_dn)*sinST;
//    }
//    if(mark_of_contract==0)//弹簧断开的计算
//    {
//        if(dis>(iUnit->m_r+jUnit->m_r))
//        {
//            continue;
//        }
//        else
//        {
//            sl_u[2]=PI+sl_u[2];//如果受压状态，旋转180度
//            if(sl_u[2]>2*PI)
//            {
//                sl_u[2]=sl_u[2]-2*PI;//如果大于360度，予以修正
//            }
//            kn1=1.00*kn;
//            ks1=1.00*ks;
//            cn1=0.6*kn;
//            cs1=0.6*ks;

//            for(i=0;i<2;i++)//求方向向量
//            {
//                e[i]=((jUnit->ut[i])-(iUnit->ut[i]))/(sqrt(dis));
//                t[0]=e[1];
//                t[1]=-e[0];
//            }
//            for(i=0;i<2;i++)
//            {
//                vn_XD[0]=((iUnit->m_velx)-(jUnit->m_velx))*e[0];//单元相对速度x方向分量
//                vn_XD[1]=((iUnit->m_vely)-(jUnit->m_vely))*e[1];//单元相对速度y方向分量
//                vs_XD[0]=(((iUnit->m_velx)-(jUnit->m_velx))*t[0])-(iUnit->m_rotate*iUnit->m_r+jUnit->m_rotate*jUnit->m_r));//相对速度
//                vs_XD[1]=(((iUnit->m_vely)-(jUnit->m_vely))*t[1])-(iUnit->m_rotate*iUnit->m_r+jUnit->m_rotate*jUnit->m_r));//相对速度
//                drt_n[i]=vn_XD[i]*dt;//相对位移增量
//                drt_s[i]=vs_XD[i]*dt;//相对位移增量
//            }
//                vn_XD[2]=sl_angel(vn_XD[0],vn_XD[1]);//法向相对速度分量
//                vs_XD[2]=sl_angel(vs_XD[0],vs_XD[1]);//切向相对速度分量
//                vn_XD[3]=slqh(vn_XD[0],vn_XD[1]);//法向相对速度大小
//                vs_XD[3]=slqh(vs_XD[0],vs_XD[1]);//切向相对速度大小
//                drt_n[2]=sl_angel(drt_n[0],drt_n[1]);//法向相对位移分量
//                drt_n[3]=slqh(drt_n[0],drt_n[1]);//法向的位移增量值
//                drt_s[2]=sl_angel(drt_s[0],drt_s[1]);//切向相对位移分量
//                drt_s[3]=slqh(drt_s[0],drt_s[1]);//切向位移增量值

//                iUnit->f_n[3]=kn*drt_n[3];//法向力增量
//                iUnit->f_s[3]=ks*drt_s[3];//切向力增量
//                iUnit->d_n[3]=cn*vn_XD[3];//法向刚度增量
//                iUnit->d_s[3]=cs*vs_XD[3];//切向刚度增量
//                iUnit->p_contract->m_fn=iUnit->p_contract->m_fn+iUnit->f_n[3]+iUnit->d_n[3];
//                iUnit->p_contract->m_fs=iUnit->p_contract->m_fs+iUnit->f_s[3]+iUnit->d_s[3];
//                iUnit->tao=(iUnit->p_contract->m_fn)*0.6+C;//库伦摩擦力
//                if(abs(iUnit->p_contract->m_fs)>abs(iUnit->tao))
//                {
//                    iUnit->p_contract->m_fs=abs(iUnit->tao)*(iUnit->p_contract->m_fs/abs(iUnit->p_contract->m_fs));//方向与原摩阻力方向一致
//                }
//                else
//                delta_Fx=iUnit->p_contract->m_fn*e[0]+iUnit->p_contract->m_fs*e[1];
//                delta_Fy=iUnit->p_contract->m_fn*e[1]-iUnit->p_contract->m_fs*e[0];
//                iUnit->p_contract->m_Fx=iUnit->p_contract->m_Fx+delta_Fx;//求x方向合力，未考虑外力和重力
//                iUnit->p_contract->m_Fy=iUnit->p_contract->m_Fy+delta_Fy;//求x方向合力，未考虑外力和重力
//                iUnit->p_contract->m_M=iUnit->p_contract->m_M+iUnit->p_contract->m_fs*iUnit->m_r;//求合力矩
//        }
//    }


//}


void CElement::cal_vtt2(CElement* p)
{
	double c1, c2;
	c1=1-alpha*drt/2;
	c2=1+alpha*drt/2;
	p->vtt2[0]=(p->vt_t2[0]*c1+p->m_Fxsum*drt/p->m_mass)/c2;
	p->vtt2[1]=(p->vt_t2[1]*c1+(p->m_Fysum/p->m_mass+g)*drt)/c2;
	p->mtt2=(p->mt_t2*c1+(p->m_Msum/p->moment)*drt)/c2;
    p->vtt2[2]=sl_angel(p->vtt2[0],p->vtt2[1]);
	p->vtt2[3]=slqh(p->vtt2[0],p->vtt2[1]);
}

void CElement::cal_vt(CElement* p)
{
	p->m_velx=(p->vtt2[0]+p->vt_t2[0])/2;
	p->m_vely=(p->vtt2[1]+p->vt_t2[1])/2;
	p->m_rotate=(p->mtt2+p->mt_t2)/2;
}

void CElement::cal_dis(CElement* p)
{
	p->m_disx=p->vtt2[0]*drt;
	p->m_disy=p->vtt2[1]*drt;
	p->m_diso=p->mtt2*drt;
}

void CElement::cal_utt(CElement* p)
{
	p->m_x=p->m_x+p->m_disx;
	p->m_y=p->m_y+p->m_disy;
	p->m_o=p->m_o+p->m_diso;
}

void CElement::change_of_data(CElement* p)
{
	int i;
	for(i=0; i<4; i++)
	{
		p->vt_t2[i]=p->vtt2[i];
	}
}

void CElement::union_lisan(CElement* p)
{
	cal_vtt2(p);
	cal_vt(p);
	cal_dis(p);
	cal_utt(p);
	change_of_data(p);
}

